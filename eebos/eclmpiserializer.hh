/*
  This file is part of the eWoms project.

  eWoms is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  eWoms is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with eWoms.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef ECL_MPI_SERIALIZER_HH
#define ECL_MPI_SERIALIZER_HH

#include <ewoms/eclsimulators/utils/parallelrestart.hh>

namespace Ewoms {

/*! \brief Class for (de-)serializing and broadcasting data in parallel.
 *!  \details Can be called on any class with a serializeOp member. Such classes
 *!           are referred to as 'complex types' in the documentation.
*/

class EclMpiSerializer {
public:
    //! \brief Constructor.
    //! \param comm The global communicator to broadcast using
    explicit EclMpiSerializer(Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm) :
        comm_(comm)
    {}

    //! \brief (De-)serialization for simple types.
    //! \details The data handled by this depends on the underlying serialization used.
    //!          Currently you can call this for scalars, and stl containers with scalars.
    template<class T>
    void operator()(const T& data)
    {
        if constexpr (is_ptr<T>::value) {
            ptr(data);
        } else if constexpr (is_pair<T>::value) {
            pair(data);
        } else {
          if (m_op == Operation::PACKSIZE)
              m_packSize += Mpi::packSize(data, comm_);
          else if (m_op == Operation::PACK)
              Mpi::pack(data, m_buffer, m_position, comm_);
          else if (m_op == Operation::UNPACK)
              Mpi::unpack(const_cast<T&>(data), m_buffer, m_position, comm_);
        }
    }

    //! \brief Handler for vectors.
    //! \tparam T Type for vector elements
    //! \tparam complexType Whether or not T is a complex type
    //! \param data The vector to (de-)serialize
    template<class T, bool complexType = true>
    void vector(std::vector<T>& data)
    {
        auto handle = [&](auto& d)
        {
            for (auto& it : d) {
              if constexpr (is_pair<T>::value)
                  pair(it);
              else if constexpr (is_ptr<T>::value)
                  ptr(it);
              else if constexpr (!complexType)
                  (*this)(it);
              else
                  it.serializeOp(*this);
            }
        };

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), comm_);
            handle(data);
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, comm_);
            handle(data);
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, comm_);
            data.resize(size);
            handle(data);
        }
    }

    //! \brief Handler for maps.
    //! \tparam Map map type
    //! \tparam complexType Whether or not Data in map is a complex type
    //! \param map The map to (de-)serialize
    template<class Map, bool complexType = true>
    void map(Map& data)
    {
        using Key = typename Map::key_type;
        using Data = typename Map::mapped_type;

        auto handle = [&](auto& d)
        {
            if constexpr (is_vector<Data>::value)
                this->template vector<typename Data::value_type,complexType>(d);
            else if constexpr (is_ptr<Data>::value)
                ptr(d);
            else if constexpr (is_dynamic_state<Data>::value)
                d.template serializeOp<EclMpiSerializer, complexType>(*this);
            else if constexpr (complexType)
                d.serializeOp(*this);
            else
                (*this)(d);
        };

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), comm_);
            for (auto& it : data) {
                m_packSize += Mpi::packSize(it.first, comm_);
                handle(it.second);
            }
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, comm_);
            for (auto& it : data) {
                Mpi::pack(it.first, m_buffer, m_position, comm_);
                handle(it.second);
            }
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, comm_);
            for (size_t i = 0; i < size; ++i) {
                Key key;
                Mpi::unpack(key, m_buffer, m_position, comm_);
                Data entry;
                handle(entry);
                data.insert(std::make_pair(key, entry));
            }
        }
    }

    //! \brief Call this to serialize data.
    //! \tparam T Type of class to serialize
    //! \param data Class to serialize
    template<class T>
    void pack(T& data)
    {
        m_op = Operation::PACKSIZE;
        m_packSize = 0;
        data.serializeOp(*this);
        m_position = 0;
        m_buffer.resize(m_packSize);
        m_op = Operation::PACK;
        data.serializeOp(*this);
    }

    //! \brief Call this to de-serialize data.
    //! \tparam T Type of class to de-serialize
    //! \param data Class to de-serialize
    template<class T>
    void unpack(T& data)
    {
        m_position = 0;
        m_op = Operation::UNPACK;
        data.serializeOp(*this);
    }

    //! \brief Serialize and broadcast on root process, de-serialize on others.
    //! \tparam T Type of class to broadcast
    //! \param data Class to broadcast
    template<class T>
    void broadcast(T& data)
    {
        if (comm_.size() == 1)
            return;

        if (comm_.rank() == 0) {
            pack(data);
            comm_.broadcast(&m_position, 1, 0);
            comm_.broadcast(m_buffer.data(), m_position, 0);
        } else {
            comm_.broadcast(&m_packSize, 1, 0);
            m_buffer.resize(m_packSize);
            comm_.broadcast(m_buffer.data(), m_packSize, 0);
            unpack(data);
        }
    }

    //! \brief Returns current position in buffer.
    size_t position() const
    {
        return m_position;
    }

    //! \brief Returns true if we are currently doing a serialization operation.
    bool isSerializing() const
    {
        return m_op != Operation::UNPACK;
    }

protected:
    //! \brief Enumeration of operations.
    enum class Operation {
        PACKSIZE, //!< Calculating serialization buffer size
        PACK,     //!< Performing serialization
        UNPACK    //!< Performing de-serialization
    };

    //! \brief Predicate for detecting pairs.
    template<class T>
    struct is_pair {
        constexpr static bool value = false;
    };

    template<class T1, class T2>
    struct is_pair<std::pair<T1,T2>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for detecting vectors.
    template<class T>
    struct is_vector {
        constexpr static bool value = false;
    };

    template<class T1>
    struct is_vector<std::vector<T1>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for smart pointers.
    template<class T>
    struct is_ptr {
        constexpr static bool value = false;
    };

    template<class T1>
    struct is_ptr<std::shared_ptr<T1>> {
        constexpr static bool value = true;
    };

    template<class T1>
    struct is_ptr<std::unique_ptr<T1>> {
        constexpr static bool value = true;
    };

    //! \brief Predicate for DynamicState.
    template<class T>
    struct is_dynamic_state {
        constexpr static bool value = false;
    };

    template<class T1>
    struct is_dynamic_state<DynamicState<T1>> {
        constexpr static bool value = true;
    };

    //! \brief Handler for pairs.
    //! \details If data is POD or a string, we pass it to the underlying serializer,
    //!          if not we assume a complex type.
    template<class T1, class T2>
    void pair(const std::pair<T1,T2>& data)
    {
        if constexpr (std::is_pod<T1>::value || std::is_same<T1,std::string>::value)
            (*this)(data.first);
        else
            data.first.serializeOp(*this);

        if constexpr (std::is_pod<T2>::value || std::is_same<T2,std::string>::value)
            (*this)(data.second);
        else
            const_cast<T2&>(data.second).serializeOp(*this);
    }

    //! \brief Handler for smart pointers.
    //! \details If data is POD or a string, we pass it to the underlying serializer,
    //!          if not we assume a complex type.
    template<template<class T> class PtrType, class T1>
    void ptr(const PtrType<T1>& data)
    {
        bool value = data ? true : false;
        (*this)(value);
        if (m_op == Operation::UNPACK && value) {
            const_cast<PtrType<T1>&>(data).reset(new T1);
        }
        if (data)
            data->serializeOp(*this);
    }

    Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm_; //!< Communicator to broadcast using

    Operation m_op = Operation::PACKSIZE; //!< Current operation
    size_t m_packSize = 0; //!< Required buffer size after PACKSIZE has been done
    int m_position = 0; //!< Current position in buffer
    std::vector<char> m_buffer; //!< Buffer for serialized data
};

}

#endif
