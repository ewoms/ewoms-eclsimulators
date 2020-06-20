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
#include <ewoms/eclio/parser/eclipsestate/aquancon.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/summarystate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/state.hh>

#include <ewoms/common/variant.hh>

#include <dune/common/parallel/mpitraits.hh>
#include <dune/common/classname.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <utility>
#include <set>
#include <tuple>
#include <vector>
#include <map>
#include <unordered_map>

namespace Ewoms {
namespace Mpi {

// helper class to extract the type at a given position from a parameter pack
template <int targetIdx, int curIdx, class T, class ...Ts>
struct ExtractType
{
    using type = typename ExtractType<targetIdx, curIdx + 1, Ts...>::type;
};

template <int targetIdx, class T, class ...Ts>
struct ExtractType<targetIdx, targetIdx, T, Ts...>
{
    typedef T type;
};

///////////////////
// pack size, pack and unpack routines for integral constants
///////////////////

// "real" integral constant
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

// complex object that can be memcopied
template<class T>
std::size_t packSize(const T*, std::size_t l, Dune::MPIHelper::MPICommunicator comm,
                     std::integral_constant<bool, true>)
{
#if HAVE_MPI
    int size;
    MPI_Pack_size(1, Dune::MPITraits<std::size_t>::getType(), comm, &size);
    std::size_t totalSize = size;
    MPI_Pack_size(l, Dune::MPITraits<T>::getType(), comm, &size);
    return totalSize + size;
#else
    (void) comm;
    return l-l;
#endif
}

// oopsies: we don't know how we could handle non-raw objects generically
template<class T>
std::size_t packSize(const T&, Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Determination of pack size not (yet) supported for objects of non-trivial type '" + Dune::className<T>() + ".");
}

template<class T>
std::size_t packSize(const T*, std::size_t, Dune::MPIHelper::MPICommunicator,
                     std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Determination of pack size not (yet) supported for arrays of non-trivial type '" + Dune::className<T>() + ".");
}

// MPI-pack a single data object
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

// MPI-pack an array of 'l' data objects
template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm,
          std::integral_constant<bool, true>)
{
#if HAVE_MPI
    MPI_Pack(&l, 1, Dune::MPITraits<std::size_t>::getType(), buffer.data(),
             buffer.size(), &position, comm);
    MPI_Pack(data, l, Dune::MPITraits<T>::getType(), buffer.data(),
             buffer.size(), &position, comm);
#else
    (void) data;
    (void) comm;
    (void) l;
    (void) buffer;
    (void) position;
#endif
}

// oopsies: we don't know how we could handle non-raw objects generically
template<class T>
void pack(const T&,
          std::vector<char>&,
          int&,
          Dune::MPIHelper::MPICommunicator,
          std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Packing not (yet) supported for objects of non-trivial type '" + Dune::className<T>() + ".");
}

template<class T>
void pack(const T*,
          std::size_t, std::vector<char>&,
          int&,
          Dune::MPIHelper::MPICommunicator,
          std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Array packing not (yet) supported for arrays of non-trivial type '" + Dune::className<T>() + ".");
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
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm,
            std::integral_constant<bool, true>)
{
#if HAVE_MPI
    MPI_Unpack(buffer.data(), buffer.size(), &position, data, l,
               Dune::MPITraits<T>::getType(), comm);
#else
    (void) data;
    (void) comm;
    (void) l;
    (void) buffer;
    (void) position;
#endif
}

// oopsies: we don't know how we could handle non-raw objects generically
template<class T>
void unpack(T&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Unpacking not (yet) supported for objects of non-trivial type '" + Dune::className<T>() + ".");
}

template<class T>
void unpack(T*, const std::size_t&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Unpacking not (yet) supported for arrays of non-trivial type '" + Dune::className<T>() + ".");
}

///////////////////
// pack size, pack and unpack routines for data types that can be copyied using memcpy()
///////////////////

#define HANDLE_AS_TRIVIAL(T) \
  inline std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm) \
  { \
      return packSize(data, comm, std::integral_constant<bool,true>()); \
  } \
  inline void pack(const T& data, std::vector<char>& buffer, int& position, \
            Dune::MPIHelper::MPICommunicator comm) \
  { \
      pack(data, buffer, position, comm, std::integral_constant<bool,true>()); \
  } \
  inline void unpack(T& data, std::vector<char>& buffer, int& position, \
              Dune::MPIHelper::MPICommunicator comm) \
  { \
      unpack(data, buffer, position, comm, std::integral_constant<bool,true>()); \
  }


HANDLE_AS_TRIVIAL(data::Connection)
HANDLE_AS_TRIVIAL(data::CurrentControl)
HANDLE_AS_TRIVIAL(data::Rates)
HANDLE_AS_TRIVIAL(data::Segment)

#undef HANDLE_AS_TRIVIAL

///////////////////
// pack size routines for complex data types
///////////////////

// prototypes
std::size_t packSize(const char* str, Dune::MPIHelper::MPICommunicator comm);
std::size_t packSize(const std::string& str, Dune::MPIHelper::MPICommunicator comm);
template<class T>std::size_t packSize(const T* data, std::size_t l, Dune::MPIHelper::MPICommunicator comm);
template<class T> std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm);
template<class T1, class T2> std::size_t packSize(const std::pair<T1,T2>& data, Dune::MPIHelper::MPICommunicator comm);
template<class T, class A> std::size_t packSize(const std::vector<T,A>& data, Dune::MPIHelper::MPICommunicator comm);
template<class A> std::size_t packSize(const std::vector<bool,A>& data, Dune::MPIHelper::MPICommunicator comm);
template<class... Ts> std::size_t packSize(const std::tuple<Ts...>& data, Dune::MPIHelper::MPICommunicator comm);
template<class... Ts> std::size_t packSize(const Ewoms::variant<Ts...>& data, Dune::MPIHelper::MPICommunicator comm);
template<class T, class H, class KE, class A> std::size_t packSize(const std::unordered_set<T,H,KE,A>& data, Dune::MPIHelper::MPICommunicator comm);
template<class K, class C, class A> std::size_t packSize(const std::set<K,C,A>& data, Dune::MPIHelper::MPICommunicator comm);
template<class T1, class T2, class C, class A> std::size_t packSize(const std::map<T1,T2,C,A>& data, Dune::MPIHelper::MPICommunicator comm);
template<class T1, class T2, class H, class P, class A> std::size_t packSize(const std::unordered_map<T1,T2,H,P,A>& data, Dune::MPIHelper::MPICommunicator comm);
template<class T, std::size_t N>
std::size_t packSize(const std::array<T,N>& data, Dune::MPIHelper::MPICommunicator comm);
std::size_t packSize(const data::Well& data, Dune::MPIHelper::MPICommunicator comm);
std::size_t packSize(const data::CellData& data, Dune::MPIHelper::MPICommunicator comm);
std::size_t packSize(const RestartKey& data, Dune::MPIHelper::MPICommunicator comm);
std::size_t packSize(const data::Solution& data, Dune::MPIHelper::MPICommunicator comm);
std::size_t packSize(const data::WellRates& data, Dune::MPIHelper::MPICommunicator comm);
std::size_t packSize(const RestartValue& data, Dune::MPIHelper::MPICommunicator comm);

// C-Strings
inline std::size_t packSize(const char* str, Dune::MPIHelper::MPICommunicator comm)
{
#if HAVE_MPI
    int size;
    MPI_Pack_size(1, Dune::MPITraits<std::size_t>::getType(), comm, &size);
    int totalSize = size;
    MPI_Pack_size(strlen(str)+1, MPI_CHAR, comm, &size);
    return totalSize + size;
#else
    (void) str;
    (void) comm;
    return 0;
#endif
}

// C++-Strings
inline std::size_t packSize(const std::string& str, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(str.c_str(), comm);
}

// trivial and unknown objects
template<class T>
std::size_t packSize(const T* data, std::size_t l, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data, l, comm, typename std::is_trivial<T>::type());
}

template<class T>
std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data, comm, typename std::is_trivial<T>::type());
}

// std::pair
template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.first, comm) + packSize(data.second, comm);
}

// generic Ewoms::optional objects
template<class T>
std::size_t packSize(const Ewoms::optional<T>& data, Dune::MPIHelper::MPICommunicator comm)
{
    bool hasValue = data.operator bool();
    std::size_t pkSize = packSize(hasValue, comm);
    if (hasValue)
        pkSize += packSize(*data, comm);
    return pkSize;
}

// generic std::vector objects
template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    if (std::is_trivial<T>::value)
        // size written automatically
        return packSize(data.data(), data.size(), comm);

    std::size_t size = packSize(data.size(), comm);

    for (const auto& entry: data)
        size += packSize(entry, comm);

    return size;
}

// generic std::vector<bool> objects with a special allocator
template<class A>
std::size_t packSize(const std::vector<bool,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    bool entry;
    return packSize(data.size(), comm) + data.size()*packSize(entry,comm);
}

// std::tuple objects
template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, std::size_t>::type
packSizeTupleEntry(const Tuple&, Dune::MPIHelper::MPICommunicator)
{
    return 0;
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, std::size_t>::type
packSizeTupleEntry(const Tuple& tuple, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(std::get<I>(tuple), comm) + packSizeTupleEntry<I+1>(tuple, comm);
}

template<class... Ts>
std::size_t packSize(const std::tuple<Ts...>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSizeTupleEntry(data, comm);
}

// Ewoms::variant objects
template<int staticIdx, class... Ts>
std::enable_if_t<staticIdx == sizeof...(Ts), std::size_t>
packSizeVariantContent(const Ewoms::variant<Ts...>& data,
                       Dune::MPIHelper::MPICommunicator comm)
{ return 0; }

template<int staticIdx, class... Ts>
std::enable_if_t<(staticIdx < sizeof...(Ts)), std::size_t>
packSizeVariantContent(const Ewoms::variant<Ts...>& data,
                       Dune::MPIHelper::MPICommunicator comm)
{
    if (Ewoms::variantIndex(data) < staticIdx)
        return packSizeVariantContent<staticIdx + 1, Ts...>(data, comm);
    else
        return packSize(Ewoms::get<typename ExtractType<staticIdx, 0, Ts...>::type>(data), comm);
}

template<class... Ts>
std::size_t packSize(const Ewoms::variant<Ts...>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    int dummy;
    return packSize(dummy, comm) + packSizeVariantContent<0, Ts...>(data, comm);
}

// std::unordered_set
template<class T, class H, class KE, class A>
std::size_t packSize(const std::unordered_set<T,H,KE,A>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry : data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

// std::set
template<class K, class C, class A>
std::size_t packSize(const std::set<K,C,A>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry : data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

// std::map
template<class T1, class T2, class C, class A>
std::size_t packSize(const std::map<T1,T2,C,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry: data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

// std::unordered_map
template<class T1, class T2, class H, class P, class A>
std::size_t packSize(const std::unordered_map<T1,T2,H,P,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry: data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

// std::array
template<class T, std::size_t N>
std::size_t packSize(const std::array<T,N>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return N*packSize(data[0], comm);
}

// Ewoms::data::Well
inline std::size_t packSize(const data::Well& data, Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.rates, comm);
    size += packSize(data.bhp, comm) + packSize(data.thp, comm);
    size += packSize(data.temperature, comm);
    size += packSize(data.control, comm);
    size += packSize(data.connections, comm);
    size += packSize(data.segments, comm);
    size += packSize(data.current_control, comm);
    return size;
}

// Ewoms::data::CellData
inline std::size_t packSize(const data::CellData& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.dim, comm) + packSize(data.data, comm) + packSize(data.target, comm);
}

// Ewoms::RestartKey
inline std::size_t packSize(const RestartKey& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.key, comm) + packSize(data.dim, comm) + packSize(data.required, comm);
}

// Ewoms::data::Solution
inline std::size_t packSize(const data::Solution& data, Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    return packSize(static_cast<const std::map< std::string, data::CellData>&>(data), comm);
}

// Ewoms::data::WellRates
inline std::size_t packSize(const data::WellRates& data, Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    return packSize(static_cast<const std::map< std::string, data::Well>&>(data), comm);
}

// Ewoms::RestartValue
inline std::size_t packSize(const RestartValue& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.solution, comm) + packSize(data.wells, comm) + packSize(data.extra, comm);
}

///////////////////
// pack routines
///////////////////

// prototypes
void pack(const char* str, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void pack(const std::string& str, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T> void pack(const T& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T> void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T1, class T2> void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T, class A> void pack(const std::vector<T, A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class K, class C, class A> void pack(const std::set<K,C,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T, class H, class KE, class A> void pack(const std::unordered_set<T,H,KE,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T, size_t N> void pack(const std::array<T,N>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class A> void pack(const std::vector<bool,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class... Ts> void pack(const std::tuple<Ts...>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class... Ts> void pack(const Ewoms::variant<Ts...>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T1, class T2, class C, class A> void pack(const std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T1, class T2, class H, class P, class A> void pack(const std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void pack(const data::Well& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void pack(const RestartKey& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void pack(const data::CellData& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void pack(const data::Solution& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void pack(const data::WellRates& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void pack(const RestartValue& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);

// character array
inline void pack(const char* str,
                 std::vector<char>& buffer,
                 int& position,
                 Dune::MPIHelper::MPICommunicator comm)
{
#if HAVE_MPI
    std::size_t length = strlen(str)+1;
    MPI_Pack(&length, 1, Dune::MPITraits<std::size_t>::getType(), buffer.data(),
        buffer.size(), &position, comm);
    MPI_Pack(str, strlen(str)+1, MPI_CHAR, buffer.data(), buffer.size(),
         &position, comm);
#else
    (void) str;
    (void) comm;
    (void) buffer;
    (void) position;
#endif
}

// std::string
inline void pack(const std::string& str, std::vector<char>& buffer, int& position,
                 Dune::MPIHelper::MPICommunicator comm)
{
    pack(str.c_str(), buffer, position, comm);
}


// trivial or unknown objects
template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data, buffer, position, comm, typename std::is_trivial<T>::type());
}

// array of trivial or unknown objects
template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data, l, buffer, position, comm, typename std::is_trivial<T>::type());
}

// std::pair
template<class T1, class T2>
void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.first, buffer, position, comm);
    pack(data.second, buffer, position, comm);
}

// Ewoms::optional
template<class T>
void pack(const Ewoms::optional<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    bool hasValue = data.operator bool();
    pack(hasValue, buffer, position, comm);
    if (hasValue)
        pack(*data, buffer, position, comm);
}

// std::vector with custom allocator
template<class T, class A>
void pack(const std::vector<T, A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    if (std::is_trivial<T>::value)
    {
        // size written automatically
        pack(data.data(), data.size(), buffer, position, comm);
        return;
    }

    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
        pack(entry, buffer, position, comm);
}

// std::set
template<class K, class C, class A>
void pack(const std::set<K,C,A>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry : data)
    {
        pack(entry, buffer, position, comm);
    }
}

// std::unordered_set
template<class T, class H, class KE, class A>
void pack(const std::unordered_set<T,H,KE,A>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry : data)
    {
        pack(entry, buffer, position, comm);
    }
}

// std::array
template<class T, size_t N>
void pack(const std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    for (const T& entry : data)
        pack(entry, buffer, position, comm);
}

// boolean std::vector with custom allocator
template<class A>
void pack(const std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);
    for (const auto& entry : data) {
        bool b = entry;
        pack(b, buffer, position, comm);
    }
}

//  std::tuple
template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, void>::type
packTupleEntry(const Tuple&, std::vector<char>&, int&,
                      Dune::MPIHelper::MPICommunicator)
{
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, void>::type
packTupleEntry(const Tuple& tuple, std::vector<char>& buffer,
                 int& position, Dune::MPIHelper::MPICommunicator comm)
{
    pack(std::get<I>(tuple), buffer, position, comm);
    packTupleEntry<I+1>(tuple, buffer, position, comm);
}

template<class... Ts>
void pack(const std::tuple<Ts...>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm)
{
    packTupleEntry(data, buffer, position, comm);
}

// Ewoms::variant objects
template<int staticIdx, class... Ts>
std::enable_if_t<(staticIdx >= sizeof...(Ts)), void>
packVariantContent(const Ewoms::variant<Ts...>& data,
                   std::vector<char>& buffer,
                   int& position,
                   Dune::MPIHelper::MPICommunicator comm)
{
}

template<int staticIdx, class... Ts>
std::enable_if_t<(staticIdx < sizeof...(Ts)), void>
packVariantContent(const Ewoms::variant<Ts...>& data,
                   std::vector<char>& buffer,
                   int& position,
                   Dune::MPIHelper::MPICommunicator comm)
{
    if (Ewoms::variantIndex(data) != staticIdx)
        packVariantContent<staticIdx + 1, Ts...>(data, buffer, position, comm);
    else
        pack(Ewoms::get<typename ExtractType<staticIdx, 0, Ts...>::type>(data), buffer, position, comm);
}

template<class... Ts>
void pack(const Ewoms::variant<Ts...>& data,
          std::vector<char>& buffer,
          int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(Ewoms::variantIndex(data), buffer, position, comm);
    packVariantContent<0>(data, buffer, position, comm);
}

// std::map
template<class T1, class T2, class C, class A>
void pack(const std::map<T1,T2,C,A>& data,
          std::vector<char>& buffer,
          int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
    {
        pack(entry, buffer, position, comm);
    }
}

// std::unordered_map
template<class T1, class T2, class H, class P, class A>
void pack(const std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
    {
        pack(entry, buffer, position, comm);
    }
}

// Ewoms::data::Well
inline void pack(const data::Well& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.rates, buffer, position, comm);
    pack(data.bhp, buffer, position, comm);
    pack(data.thp, buffer, position, comm);
    pack(data.temperature, buffer, position, comm);
    pack(data.control, buffer, position, comm);
    pack(data.connections, buffer, position, comm);
    pack(data.segments, buffer, position, comm);
    pack(data.current_control, buffer, position, comm);
}

// Ewoms::RestartKey
inline void pack(const RestartKey& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.key, buffer, position, comm);
    pack(data.dim, buffer, position, comm);
    pack(data.required, buffer, position, comm);
}

// Ewoms::data::CellData
inline void pack(const data::CellData& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.dim, buffer, position, comm);
    pack(data.data, buffer, position, comm);
    pack(data.target, buffer, position, comm);
}

// Ewoms::data::Solution
inline void pack(const data::Solution& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    pack(static_cast<const std::map< std::string, data::CellData>&>(data),
         buffer, position, comm);
}

// Ewoms::data::WellRates
inline void pack(const data::WellRates& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    pack(static_cast<const std::map< std::string, data::Well>&>(data),
         buffer, position, comm);
}

// Ewoms::data::RestartValue
inline void pack(const RestartValue& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.solution, buffer, position, comm);
    pack(data.wells, buffer, position, comm);
    pack(data.extra, buffer, position, comm);
}

///////////////////
// unpack routines
///////////////////

// prototypes
template<class T> void unpack(T& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T> void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void unpack(std::string& str, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T, class A> void unpack(std::vector<T,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class A> void unpack(std::vector<bool,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class... Ts> void unpack(std::tuple<Ts...>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class... Ts> void unpack(Ewoms::variant<Ts...>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class K, class C, class A> void unpack(std::set<K,C,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T, class H, class KE, class A> void unpack(std::unordered_set<T,H,KE,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T, size_t N> void unpack(std::array<T,N>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T1, class T2, class C, class A> void unpack(std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T1, class T2, class H, class P, class A> void unpack(std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void unpack(data::Well& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void unpack(RestartKey& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void unpack(data::CellData& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void unpack(data::Solution& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void unpack(data::WellRates& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
void unpack(RestartValue& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);
template<class T1, class T2> void unpack(std::pair<T1,T2>& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm);

// trivial or unknown objects
template<class T>
void unpack(T& data,
            std::vector<char>& buffer,
            int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data, buffer, position, comm, typename std::is_trivial<T>::type());
}

// arrays of trivial or unknown objects
template<class T>
void unpack(T* data,
            const std::size_t& l,
            std::vector<char>& buffer,
            int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data, l, buffer, position, comm, typename std::is_trivial<T>::type());
}


// C-Strings
inline void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
#if HAVE_MPI
    MPI_Unpack(buffer.data(), buffer.size(), &position, const_cast<char*>(str), length, MPI_CHAR, comm);
#else
    (void) str;
    (void) comm;
    (void) length;
    (void) buffer;
    (void) position;
#endif
}

// C++ strings
inline void unpack(std::string& str, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t length=0;
    unpack(length, buffer, position, comm);
    std::vector<char> cStr(length, '\0');
    unpack(cStr.data(), length, buffer, position, comm);
    str.clear();
    str.append(cStr.data());
}

// Ewoms::optional
template<class T>
void unpack(Ewoms::optional<T>&data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    bool hasValue;
    unpack(hasValue, buffer, position, comm);
    if (hasValue)
        unpack(*data, buffer, position, comm);
    else
        data = Ewoms::nullopt;
}

// dynamic arrays with arbitrary allocators
template<class T, class A>
void unpack(std::vector<T,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t length = 0;
    unpack(length, buffer, position, comm);
    data.resize(length);

    if (std::is_trivial<T>::value)
    {
        unpack(data.data(), data.size(), buffer, position, comm);
        return;
    }

    for (auto& entry: data)
        unpack(entry, buffer, position, comm);
}

// dynamic boolean arrays with custom allocators
template<class A>
void unpack(std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    size_t size;
    unpack(size, buffer, position, comm);
    data.clear();
    data.reserve(size);
    for (size_t i = 0; i < size; ++i) {
        bool entry;
        unpack(entry, buffer, position, comm);
        data.push_back(entry);
    }
}

// std::tuple
template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, void>::type
unpackTupleEntry(Tuple&, std::vector<char>&, int&,
                   Dune::MPIHelper::MPICommunicator)
{
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, void>::type
unpackTupleEntry(Tuple& tuple, std::vector<char>& buffer,
                   int& position, Dune::MPIHelper::MPICommunicator comm)
{
    unpack(std::get<I>(tuple), buffer, position, comm);
    unpackTupleEntry<I+1>(tuple, buffer, position, comm);
}

template<class... Ts>
void unpack(std::tuple<Ts...>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm)
{
    unpackTupleEntry(data, buffer, position, comm);
}

// Ewoms::variant objects
template<int staticIdx, class... Ts>
typename std::enable_if<staticIdx == sizeof...(Ts), void>::type
unpackVariantContent(int dynamicIdx,
                     Ewoms::variant<Ts...>& data,
                     std::vector<char>& buffer,
                     int& position,
                     Dune::MPIHelper::MPICommunicator comm)
{
    throw std::runtime_error("Cannot unpack element "+std::to_string(dynamicIdx)+
                             " for a"+std::to_string(sizeof...(Ts))+"-sized variant");
}

template<int staticIdx, class... Ts>
typename std::enable_if<(staticIdx < sizeof...(Ts)), void>::type
unpackVariantContent(int dynamicIdx,
                     Ewoms::variant<Ts...>& data,
                     std::vector<char>& buffer,
                     int& position,
                     Dune::MPIHelper::MPICommunicator comm)
{
    if (dynamicIdx != staticIdx)
        unpackVariantContent<staticIdx + 1, Ts...>(dynamicIdx, data, buffer, position, comm);
    else {
        using VT = typename ExtractType<staticIdx, 0, Ts...>::type;
        VT val;
        unpack(val, buffer, position, comm);
        data = val;
    }
}

template<class... Ts>
void unpack(Ewoms::variant<Ts...>& data,
            std::vector<char>& buffer,
            int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int which;
    unpack(which, buffer, position, comm);
    unpackVariantContent<0, Ts...>(which, data, buffer, position, comm);
}

// std::set
template<class K, class C, class A>
void unpack(std::set<K,C,A>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = 0;
    unpack(size, buffer, position, comm);

    for (;size>0; size--)
    {
        K entry;
        unpack(entry, buffer, position, comm);
        data.insert(entry);
    }
}

// std::unordered_set
template<class T, class H, class KE, class A>
void unpack(std::unordered_set<T,H,KE,A>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size=0;
    unpack(size, buffer, position, comm);

    for (;size>0; size--)
    {
        T entry;
        unpack(entry, buffer, position, comm);
        data.insert(entry);
    }
}

// std::array
template<class T, size_t N>
void unpack(std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    for (T& entry : data)
        unpack(entry, buffer, position, comm);
}

// std::map
template<class T1, class T2, class C, class A>
void unpack(std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size=0;
    unpack(size, buffer, position, comm);

    for (;size>0; size--)
    {
        std::pair<T1,T2> entry;
        unpack(entry, buffer, position, comm);
        data.insert(entry);
    }
}

// std::unordered_map
template<class T1, class T2, class H, class P, class A>
void unpack(std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size=0;
    unpack(size, buffer, position, comm);

    for (; size > 0; --size)
    {
        std::pair<T1,T2> entry;
        unpack(entry, buffer, position, comm);
        data.insert(entry);
    }
}

// Ewoms::data::Well
inline void unpack(data::Well& data, std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.rates, buffer, position, comm);
    unpack(data.bhp, buffer, position, comm);
    unpack(data.thp, buffer, position, comm);
    unpack(data.temperature, buffer, position, comm);
    unpack(data.control, buffer, position, comm);
    unpack(data.connections, buffer, position, comm);
    unpack(data.segments, buffer, position, comm);
    unpack(data.current_control, buffer, position, comm);
}

// Ewoms::RestartKey
inline void unpack(RestartKey& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.key, buffer, position, comm);
    unpack(data.dim, buffer, position, comm);
    unpack(data.required, buffer, position, comm);
}

// Ewoms::CellData
inline void unpack(data::CellData& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.dim, buffer, position, comm);
    unpack(data.data, buffer, position, comm);
    unpack(data.target, buffer, position, comm);
}

// Ewoms::Solution
inline void unpack(data::Solution& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    unpack(static_cast<std::map< std::string, data::CellData>&>(data),
           buffer, position, comm);
}

// Ewoms::WellRates
inline void unpack(data::WellRates& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    unpack(static_cast<std::map< std::string, data::Well>&>(data),
           buffer, position, comm);
}

// Ewoms::RestartValue
inline void unpack(RestartValue& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.solution, buffer, position, comm);
    unpack(data.wells, buffer, position, comm);
    unpack(data.extra, buffer, position, comm);
}

// std::pair
template<class T1, class T2>
void unpack(std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.first, buffer, position, comm);
    unpack(data.second, buffer, position, comm);
}

} // end namespace Mpi

inline RestartValue loadParallelRestart(const EclipseIO* eclIO,
                                        SummaryState& summaryState,
                                        Action::State& actionState,
                                        const std::vector<Ewoms::RestartKey>& solutionKeys,
                                        const std::vector<Ewoms::RestartKey>& extraKeys,
                                        Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm)
{
#if HAVE_MPI
    data::Solution sol;
    data::Wells wells;
    RestartValue restartValues(sol, wells);

    if (eclIO)
    {
        assert(comm.rank() == 0);
        restartValues = eclIO->loadRestart(actionState, summaryState, solutionKeys, extraKeys);
        int packedSize = Mpi::packSize(restartValues, comm);
        std::vector<char> buffer(packedSize);
        int position=0;
        Mpi::pack(restartValues, buffer, position, comm);
        comm.broadcast(&position, 1, 0);
        comm.broadcast(buffer.data(), position, 0);
    }
    else
    {
        int bufferSize{};
        comm.broadcast(&bufferSize, 1, 0);
        std::vector<char> buffer(bufferSize);
        comm.broadcast(buffer.data(), bufferSize, 0);
        int position{};
        Mpi::unpack(restartValues, buffer, position, comm);
    }
    return restartValues;
#else
    (void) comm;
    return eclIO->loadRestart(summaryState, solutionKeys, extraKeys);
#endif
}

} // end namespace Ewoms


#endif // PARALLEL_RESTART_HH
