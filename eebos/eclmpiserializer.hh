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
    explicit EclMpiSerializer(Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm) :
        comm_(comm)
    {}

    template<class T>
    std::size_t packSize(const T& data) {
        return Mpi::packSize(data, comm_);
    }

    template<class T>
    void pack(const T& data, std::vector<char>& buffer, int& pos) {
        Mpi::pack(data, buffer, pos, comm_);
    }

    template<class T>
    void unpack(T& data, std::vector<char>& buffer, int& pos) {
        Mpi::unpack(data, buffer, pos, comm_);
    }

    template<class T>
    void staticBroadcast()
    {
        if (comm_.size() == 1)
            return;

#if HAVE_MPI
        if (comm_.rank() == 0) {
            size_t size = T::packSize(*this);
            std::vector<char> buffer(size);
            int position = 0;
            T::pack(buffer, position, *this);
            comm_.broadcast(&position, 1, 0);
            comm_.broadcast(buffer.data(), position, 0);
        } else {
            int size;
            comm_.broadcast(&size, 1, 0);
            std::vector<char> buffer(size);
            comm_.broadcast(buffer.data(), size, 0);
            int position = 0;
            T::unpack(buffer, position, *this);
        }
#endif
    }

    template<class T>
    void broadcast(T& data)
    {
        if (comm_.size() == 1)
            return;

#if HAVE_MPI
        if (comm_.rank() == 0) {
            size_t size = data.packSize(*this);
            std::vector<char> buffer(size);
            int position = 0;
            data.pack(buffer, position, *this);
            comm_.broadcast(&position, 1, 0);
            comm_.broadcast(buffer.data(), position, 0);
        } else {
            int size;
            comm_.broadcast(&size, 1, 0);
            std::vector<char> buffer(size);
            comm_.broadcast(buffer.data(), size, 0);
            int position = 0;
            data.unpack(buffer, position, *this);
        }
#endif
    }

protected:
    Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm_;
};

}

#endif
