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

class EclMpiSerializer {
public:
    enum class Operation {
        PACKSIZE,
        PACK,
        UNPACK
    };

    explicit EclMpiSerializer(Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm) :
        comm_(comm)
    {}

    template<class T>
    void operator()(const T& data)
    {
        if (m_op == Operation::PACKSIZE)
            m_packSize += Mpi::packSize(data, comm_);
        else if (m_op == Operation::PACK)
            Mpi::pack(data, m_buffer, m_position, comm_);
        else if (m_op == Operation::UNPACK)
            Mpi::unpack(const_cast<T&>(data), m_buffer, m_position, comm_);
    }

    template<class T>
    void vector(std::vector<T>& data)
    {
        static_assert(!std::is_pod<T>::value, "Do not call this for POD vectors");
        auto handle = [&](auto& d)
        {
            for (auto& it : d) {
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

    template<class T>
    void unpack(T& data)
    {
        m_position = 0;
        m_op = Operation::UNPACK;
        data.serializeOp(*this);
    }

    template<class T>
    void broadcast(T& data)
    {
        if (comm_.size() == 1)
            return;

#if HAVE_MPI
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
#endif
    }

    size_t position() const
    {
        return m_position;
    }

protected:
    Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm_;

    Operation m_op = Operation::PACKSIZE;
    size_t m_packSize = 0;
    int m_position = 0;
    std::vector<char> m_buffer;
};

}

#endif
