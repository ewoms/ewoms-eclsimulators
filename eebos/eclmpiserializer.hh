// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
*/
#ifndef ECL_MPI_SERIALIZER_HH
#define ECL_MPI_SERIALIZER_HH

#include <ewoms/eclsimulators/utils/parallelrestart.hh>

#include <ewoms/eclio/parser/eclipsestate/schedule/dynamicstate.hh>

#include <ewoms/common/optional.hh>
#include <ewoms/common/variant.hh>

namespace Ewoms {

/*! \brief Class for (de-)serializing and broadcasting data for MPI runs.
 *!  \details Can be called on any class with a serializeOp member. Such classes
 *!           are referred to as 'complex types' in the documentation.
*/

class EclMpiSerializer
{
    //! \brief Predicate for detecting pairs.
    template<class T>
    struct is_pair : public std::false_type {};

    template<class T1, class T2>
    struct is_pair<std::pair<T1,T2>> : public std::true_type {};

    //! \brief Predicate for Ewoms::optional.
    template<class T>
    struct is_optional : public std::false_type {};

    template<class T1>
    struct is_optional<Ewoms::optional<T1>> : public std::true_type {};

    //! \brief Predicate for detecting variants.
    template<class T>
    struct is_variant : public std::false_type {};

    template<class... Ts>
    struct is_variant<Ewoms::variant<Ts...>> : public std::true_type {};

    //! \brief Predicate for detecting vectors.
    template<class T>
    struct is_vector : public std::false_type {};

    template<class T1>
    struct is_vector<std::vector<T1>> : public std::true_type {};

    //! \brief Predicate for smart pointers.
    template<class T>
    struct is_ptr : public std::false_type {};

    template<class T1>
    struct is_ptr<std::shared_ptr<T1>> : public std::true_type {};

    template<class T1, class Deleter>
    struct is_ptr<std::unique_ptr<T1, Deleter>> : public std::true_type {};

    //! \brief Predicate for DynamicState.
    template<class T>
    struct is_dynamic_state : public std::false_type { };

    template<class T1>
    struct is_dynamic_state<DynamicState<T1>>  : public std::true_type { };

    //! \brief Predicate for determining wheter an onbject has a template method called serializeOp().
    template <typename... Ts>
    using void_t = void;

    template <typename T, typename = void>
    struct hasSerializeOp : public std::false_type {};

    template <typename T>
    struct hasSerializeOp<T, void_t<decltype(&T::template serializeOp<EclMpiSerializer>)>> : public std::true_type {};

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
        bool isHandled = false;

        isHandled = isHandled || ptr_(data, is_ptr<T>());
        isHandled = isHandled || pair_(data, is_pair<T>());
        isHandled = isHandled || variant_(data, is_variant<T>());
        isHandled = isHandled || optional_(data, is_optional<T>());
        isHandled = isHandled || vector_(data, is_vector<T>());
        isHandled = isHandled || serializeOp_(data, hasSerializeOp<T>());

        if (!isHandled)
            primitive(data);
    }

    template<class T>
    bool primitive(const T& data)
    {
#if HAVE_MPI
        if (m_op == Operation::PACKSIZE)
            m_packSize += Mpi::packSize(data, comm_);
        else if (m_op == Operation::PACK)
            Mpi::pack(data, m_buffer, m_position, comm_);
        else if (m_op == Operation::UNPACK)
            Mpi::unpack(const_cast<T&>(data), m_buffer, m_position, comm_);
#endif

        return true;
    }

    //! \brief Handler for vectors.
    //! \tparam T Type for vector elements
    //! \param data The vector to (de-)serialize
    template<class T, bool dummy = true>
    void vector(const std::vector<T>& data)
    { vector_(data, std::true_type()); }

private:
    template<class T>
    bool vector_(T& data, std::true_type)
    {
#if HAVE_MPI
        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), comm_);
            for (const auto& it : data)
                (*this)(it);
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, comm_);
            for (const auto& it : data)
                (*this)(it);
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, comm_);
            const_cast<std::remove_const_t<T>&>(data).resize(size);
            for (const auto& it : data)
                (*this)(it);
        }
#endif

        return true;
    }

    template<class T>
    bool vector_(T& data, std::false_type)
    { return false; }

private:
    //! \brief Handler for Ewoms::optional.
    //! \tparam T Type for data
    //! \param data The optional to (de-)serialize
    template<class T>
    bool optional_(const Ewoms::optional<T>& data, std::true_type)
    {
#if HAVE_MPI
        if (m_op == Operation::PACKSIZE) {
            bool hasValue = data.operator bool();
            if (hasValue)
                this->operator()(*data);

            m_packSize += Mpi::packSize(hasValue, comm_);
        } else if (m_op == Operation::PACK) {
            bool hasValue = data.operator bool();
            Mpi::pack(hasValue, m_buffer, m_position, comm_);
            if (hasValue) {
                (*this)(*data);
            }
        } else if (m_op == Operation::UNPACK) {
            bool hasValue;
            Mpi::unpack(hasValue,
                        m_buffer,
                        m_position,
                        comm_);

            if (!hasValue)
                const_cast<Ewoms::optional<T>&>(data) = Ewoms::nullopt;
            else {
                T val;
                (*this)(val);
                const_cast<Ewoms::optional<T>&>(data) = Ewoms::optional<T>(val);
            }
        }
#endif

        return true;
    }

    template<class T>
    bool optional_(const T&, std::false_type)
    { return false; }

    //! \brief Handler for Ewoms::variant<> with three types
    template<class... T>
    bool variant_(const Ewoms::variant<T...>& data, std::true_type)
    {
#if HAVE_MPI
        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data, comm_);
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data, m_buffer, m_position, comm_);
        } else if (m_op == Operation::UNPACK) {
            Mpi::unpack(const_cast<Ewoms::variant<T...>&>(data),
                        m_buffer,
                        m_position,
                        comm_);
        }
#endif

        return true;
    }

    template<class T>
    bool variant_(const T&, std::false_type)
    { return false; }

public:
    //! \brief Handler for maps.
    //! \tparam Map map type
    //! \tparam complexType Whether or not Data in map is a complex type
    //! \param map The map to (de-)serialize
    template<class Map, bool complexType = true>
    bool map(Map& data)
    {
#if HAVE_MPI
        using Key = typename Map::key_type;
        using Data = typename Map::mapped_type;

        if (m_op == Operation::PACKSIZE) {
            m_packSize += Mpi::packSize(data.size(), comm_);
            for (auto& it : data) {
                m_packSize += Mpi::packSize(it.first, comm_);
                (*this)(it.second);
            }
        } else if (m_op == Operation::PACK) {
            Mpi::pack(data.size(), m_buffer, m_position, comm_);
            for (auto& it : data) {
                Mpi::pack(it.first, m_buffer, m_position, comm_);
                (*this)(it.second);
            }
        } else if (m_op == Operation::UNPACK) {
            size_t size;
            Mpi::unpack(size, m_buffer, m_position, comm_);
            for (size_t i = 0; i < size; ++i) {
                Key key;
                Mpi::unpack(key, m_buffer, m_position, comm_);
                Data entry;
                (*this)(entry);
                data.insert(std::make_pair(key, entry));
            }
        }
#endif

        return true;
    }

    //! \brief Handler for pairs.
    //! \details If data is POD or a string, we pass it to the underlying serializer,
    //!          if not we assume a complex type.
    template<class T1, class T2>
    void pair(const std::pair<T1,T2>& data)
    { pair_(data, std::true_type()); }

private:
    template<class T1, class T2>
    bool pair_(const std::pair<T1,T2>& data, std::true_type)
    {
        (*this)(data.first);
        (*this)(data.second);

        return true;
    }

    template<class T>
    bool pair_(const T& data, std::false_type)
    { return false; }

public:
    //! \brief Handler for smart pointers.
    //! \details If data is POD or a string, we pass it to the underlying serializer,
    //!          if not we assume a complex type.
    template<class PtrType>
    void ptr(const PtrType& data)
    { ptr_(data, std::true_type()); }

private:
    template<class PtrType>
    bool ptr_(const PtrType& data, std::true_type)
    {
        using T1 = typename PtrType::element_type;
        bool value = data ? true : false;
        (*this)(value);
        if (m_op == Operation::UNPACK && value)
            const_cast<PtrType&>(data).reset(new T1);
        if (data)
            data->serializeOp(*this);

        return true;
    }

    template<class PtrType>
    bool ptr_(const PtrType& data, std::false_type)
    { return false; }


    template<class T>
    bool serializeOp_(const T& data, std::true_type)
    {
        const_cast<T&>(data).serializeOp(*this);

        return true;
    }

    template<class PtrType>
    bool serializeOp_(const PtrType& data, std::false_type)
    { return false; }

public:
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
            m_packSize = m_position;
            try {
                pack(data);
                m_packSize = m_position;
                comm_.broadcast(&m_packSize, 1, 0);
                comm_.broadcast(m_buffer.data(), m_position, 0);
            } catch (...) {
                m_packSize = std::numeric_limits<size_t>::max();
                comm_.broadcast(&m_packSize, 1, 0);
                throw;
            }
        } else {
            comm_.broadcast(&m_packSize, 1, 0);
            if (m_packSize == std::numeric_limits<size_t>::max()) {
                throw std::runtime_error("Error detected in parallel serialization");
            }
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

    Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm_; //!< Communicator to broadcast using

    Operation m_op = Operation::PACKSIZE; //!< Current operation
    size_t m_packSize = 0; //!< Required buffer size after PACKSIZE has been done
    int m_position = 0; //!< Current position in buffer
    std::vector<char> m_buffer; //!< Buffer for serialized data
};

}

#endif
