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
#ifndef PARALLEL_RESTART_HH
#define PARALLEL_RESTART_HH

#if HAVE_MPI
#include <mpi.h>
#endif

#include <ewoms/eclio/output/restartvalue.hh>
#include <ewoms/eclio/output/eclipseio.hh>
#include <ewoms/eclio/output/summary.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <set>
#include <tuple>
#include <vector>
#include <map>
#include <unordered_map>

namespace Ewoms
{

namespace Mpi
{
template<class T>
std::size_t packSize(const T*, std::size_t, Dune::MPIHelper::MPICommunicator,
                     std::integral_constant<bool, false>);

template<class T>
std::size_t packSize(const T*, std::size_t l, Dune::MPIHelper::MPICommunicator comm,
                     std::integral_constant<bool, true>);

template<class T>
std::size_t packSize(const T* data, std::size_t l, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const T&, Dune::MPIHelper::MPICommunicator,
                     std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
std::size_t packSize(const T&, Dune::MPIHelper::MPICommunicator comm,
                     std::integral_constant<bool, true>)
{
#if HAVE_MPI
    int size{};
    MPI_Pack_size(1, Dune::MPITraits<T>::getType(), comm, &size);
    return size;
#else
    (void) comm;
    return 0;
#endif
}

template<class T>
std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.first, comm) + packSize(data.second, comm);
}

template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class K, class C, class A>
std::size_t packSize(const std::set<K,C,A>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class T, class H, class KE, class A>
std::size_t packSize(const std::unordered_set<T,H,KE,A>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class A>
std::size_t packSize(const std::vector<bool,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class... Ts>
std::size_t packSize(const std::tuple<Ts...>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T, std::size_t N>
std::size_t packSize(const std::array<T,N>& data, Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const char* str, Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const std::string& str, Dune::MPIHelper::MPICommunicator comm);

inline void pack(const std::string& str, std::vector<char>& buffer, int& position,
                 Dune::MPIHelper::MPICommunicator comm)
{
    pack(str.c_str(), buffer, position, comm);
}

template<class T1, class T2, class C, class A>
std::size_t packSize(const std::map<T1,T2,C,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
std::size_t packSize(const std::unordered_map<T1,T2,H,P,A>& data, Dune::MPIHelper::MPICommunicator comm);

////// pack routines

template<class T>
void pack(const T*, std::size_t, std::vector<char>&, int&,
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const T&, std::vector<char>&, int&,
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>)
{
#if HAVE_MPI
    MPI_Pack(&data, 1, Dune::MPITraits<T>::getType(), buffer.data(),
             buffer.size(), &position, comm);
#else
    (void) data;
    (void) comm;
    (void) buffer;
    (void) position;
#endif
}

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.first, buffer, position, comm);
    pack(data.second, buffer, position, comm);
}

template<class T, class A>
void pack(const std::vector<T,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class A>
void pack(const std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class... Ts>
void pack(const std::tuple<Ts...>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class K, class C, class A>
void pack(const std::set<K,C,A>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T, class H, class KE, class A>
void pack(const std::unordered_set<T,H,KE,A>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T, size_t N>
void pack(const std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
void pack(const std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
void pack(const std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

void pack(const char* str, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

/// unpack routines

template<class T>
void unpack(T*, const std::size_t&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm,
            std::integral_constant<bool, true>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(T&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>)
{
#if HAVE_MPI
    MPI_Unpack(buffer.data(), buffer.size(), &position, &data, 1,
               Dune::MPITraits<T>::getType(), comm);
#else
    (void) data;
    (void) comm;
    (void) buffer;
    (void) position;
#endif
}

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
void unpack(std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.first, buffer, position, comm);
    unpack(data.second, buffer, position, comm);
}

template<class T, class A>
void unpack(std::vector<T,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class A>
void unpack(std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class... Ts>
void unpack(std::tuple<Ts...>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class K, class C, class A>
void unpack(std::set<K,C,A>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T, class H, class KE, class A>
void unpack(std::unordered_set<T,H,KE,A>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T, size_t N>
void unpack(std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
void unpack(std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
void unpack(std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

/// prototypes for complex types

#define ADD_PACK_PROTOTYPES(T) \
  std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm); \
  void pack(const T& data, std::vector<char>& buffer, int& position, \
          Dune::MPIHelper::MPICommunicator comm); \
  void unpack(T& data, std::vector<char>& buffer, int& position, \
              Dune::MPIHelper::MPICommunicator comm);

ADD_PACK_PROTOTYPES(data::CellData)
ADD_PACK_PROTOTYPES(data::Connection)
ADD_PACK_PROTOTYPES(data::CurrentControl)
ADD_PACK_PROTOTYPES(data::Rates)
ADD_PACK_PROTOTYPES(data::Segment)
ADD_PACK_PROTOTYPES(data::Solution)
ADD_PACK_PROTOTYPES(data::Well)
ADD_PACK_PROTOTYPES(data::WellRates)
ADD_PACK_PROTOTYPES(RestartKey)
ADD_PACK_PROTOTYPES(RestartValue)
ADD_PACK_PROTOTYPES(std::string)

} // end namespace Mpi

RestartValue loadParallelRestart(const EclipseIO* eclIO, SummaryState& summaryState,
                                 const std::vector<Ewoms::RestartKey>& solutionKeys,
                                 const std::vector<Ewoms::RestartKey>& extraKeys,
                                 Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm);

} // end namespace Ewoms
#endif // PARALLEL_RESTART_HH
