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
#include <config.h>
#if HAVE_MPI
#include <mpi.h>
#endif

#include "parallelrestart.hh"
#include <ewoms/eclio/opmlog/location.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipseconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/runspec.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/fault.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/faultcollection.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/faultface.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/multregtscanner.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/nnc.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/transmult.hh>
#include <ewoms/eclio/parser/eclipsestate/initconfig/equil.hh>
#include <ewoms/eclio/parser/eclipsestate/initconfig/foamconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/initconfig/initconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/ioconfig/ioconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/ioconfig/restartconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/edit/editnnc.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actionast.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actions.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actionx.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/astnode.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/condition.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/events.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/guiderateconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/messagelimits.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/icd.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/spiralicd.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/valve.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/oilvaporizationproperties.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/rftconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/scheduletypes.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/timemap.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/tuning.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqastnode.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqdefine.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqfunction.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqinput.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqfunctiontable.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/vfpinjtable.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/vfpprodtable.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/connection.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wellconnections.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wellfoamproperties.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wellpolymerproperties.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/welltracerproperties.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wlist.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wlistmanager.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/bcconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/rockconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/simulationconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/thresholdpressure.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/aqudims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/columnschema.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/eqldims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/flattable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/jfunc.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/plymwinjtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/polyinjtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/pvtgtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/pvtotable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/regdims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rock2dtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rock2dtrtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rocktabtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/plyshlogtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/simpletable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/skprpolytable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/skprwattable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/standardcond.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablecolumn.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablecontainer.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablemanager.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tableschema.hh>
#include <ewoms/eclio/parser/eclipsestate/util/iorderset.hh>
#include <dune/common/parallel/mpitraits.hh>

#define HANDLE_AS_POD(T) \
  std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm) \
  { \
      return packSize(data, comm, std::integral_constant<bool,true>()); \
  } \
  void pack(const T& data, std::vector<char>& buffer, int& position, \
            Dune::MPIHelper::MPICommunicator comm) \
  { \
      pack(data, buffer, position, comm, std::integral_constant<bool,true>()); \
  } \
  void unpack(T& data, std::vector<char>& buffer, int& position, \
              Dune::MPIHelper::MPICommunicator comm) \
  { \
      unpack(data, buffer, position, comm, std::integral_constant<bool,true>()); \
  }

namespace
{
template<template<class, class> class Map, class Type, class Key>
std::pair<std::vector<Type>, std::vector<std::pair<Key, std::vector<int>>>>
splitDynMap(const Map<Key, Ewoms::DynamicState<Type>>& map)
{
    // we have to pack the unique ptrs separately, and use an index map
    // to allow reconstructing the appropriate structures.
    std::vector<std::pair<Key, std::vector<int>>> asMap;
    std::vector<Type> unique;
    for (const auto& it : map) {
        for (const auto& w : it.second.data()) {
            if (std::find(unique.begin(), unique.end(), w) == unique.end())
                unique.push_back(w);
        }
    }
    for (const auto& it : map) {
        std::vector<int> idxVec;
        idxVec.reserve(it.second.size()+1);
        for (const auto& w : it.second.data()) {
            auto uIt = std::find(unique.begin(), unique.end(), w);
            idxVec.push_back(uIt-unique.begin());
        }
        idxVec.push_back(it.second.initialRange());
        asMap.push_back(std::make_pair(it.first, idxVec));
    }

    return std::make_pair(unique, asMap);
}

template<class Type>
std::pair<std::vector<Type>, std::vector<int>>
splitDynState(const Ewoms::DynamicState<Type>& state)
{
    std::vector<Type> unique;
    for (const auto& w : state.data()) {
        if (std::find(unique.begin(), unique.end(), w) == unique.end())
            unique.push_back(w);
    }
    std::vector<int> idxVec;
    idxVec.reserve(state.data().size()+1);
    for (const auto& w : state.data()) {
        auto uIt = std::find(unique.begin(), unique.end(), w);
        idxVec.push_back(uIt-unique.begin());
    }
    idxVec.push_back(state.initialRange());

    return std::make_pair(unique, idxVec);
}

template<class Type>
void reconstructDynState(const std::vector<Type>& unique,
                         const std::vector<int>& idxVec,
                         Ewoms::DynamicState<Type>& result)
{
    std::vector<Type> ptrData;
    for (size_t i = 0; i < idxVec.size()-1; ++i) {
        ptrData.push_back(unique[idxVec[i]]);
    }
    result = Ewoms::DynamicState<Type>(ptrData, idxVec.back());
}

template<template<class, class> class Map, class Type, class Key>
void reconstructDynMap(const std::vector<Type>& unique,
                       const std::vector<std::pair<Key, std::vector<int>>>& asMap,
                       Map<Key, Ewoms::DynamicState<Type>>& result)
{
    for (const auto& it : asMap) {
        reconstructDynState(unique, it.second, result[it.first]);
    }
}

template<template<class, class> class Map, class Type, class Key>
std::size_t packSizeDynMap(const Map<Key, Ewoms::DynamicState<Type>>& data,
                           Dune::MPIHelper::MPICommunicator comm)
{
    auto split = splitDynMap<Map,Type,Key>(data);
    return Ewoms::Mpi::packSize(split.first, comm) +
           Ewoms::Mpi::packSize(split.second, comm);
}

template<template<class, class> class Map, class Type, class Key>
void packDynMap(const Map<Key, Ewoms::DynamicState<Type>>& data,
                std::vector<char>& buffer,
                int& position,
                Dune::MPIHelper::MPICommunicator comm)
{
    auto split = splitDynMap<Map,Type,Key>(data);
    Ewoms::Mpi::pack(split.first, buffer, position, comm);
    Ewoms::Mpi::pack(split.second, buffer, position, comm);
}

template<template<class, class> class Map, class Type, class Key>
void unpackDynMap(Map<Key, Ewoms::DynamicState<Type>>& data,
                  std::vector<char>& buffer,
                  int& position,
                  Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Type> unique;
    std::vector<std::pair<Key, std::vector<int>>> indices;
    Ewoms::Mpi::unpack(unique, buffer, position, comm);
    Ewoms::Mpi::unpack(indices, buffer, position, comm);
    reconstructDynMap<Map,Type,Key>(unique, indices, data);
}

}

namespace Ewoms
{
namespace Mpi
{
template<class T>
std::size_t packSize(const T*, std::size_t, Dune::MPIHelper::MPICommunicator,
                     std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

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

template<class T>
std::size_t packSize(const T* data, std::size_t l, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data, l, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.first, comm) + packSize(data.second, comm);
}

template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    if (std::is_pod<T>::value)
        // size written automatically
        return packSize(data.data(), data.size(), comm);

    std::size_t size = packSize(data.size(), comm);

    for (const auto& entry: data)
        size += packSize(entry, comm);

    return size;
}

template<class A>
std::size_t packSize(const std::vector<bool,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    bool entry;
    return packSize(data.size(), comm) + data.size()*packSize(entry,comm);
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, std::size_t>::type
packSize_tuple_entry(const Tuple&, Dune::MPIHelper::MPICommunicator)
{
    return 0;
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, std::size_t>::type
packSize_tuple_entry(const Tuple& tuple, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(std::get<I>(tuple), comm) + packSize_tuple_entry<I+1>(tuple, comm);
}

template<class... Ts>
std::size_t packSize(const std::tuple<Ts...>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize_tuple_entry(data, comm);
}

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

template<class Key, class Value>
std::size_t packSize(const OrderedMap<Key,Value>& data, Dune::MPIHelper::MPICommunicator comm)
{
  return packSize(data.getIndex(), comm) + packSize(data.getStorage(), comm);
}

template<class T>
std::size_t packSize(const DynamicState<T>& data, Dune::MPIHelper::MPICommunicator comm)
{

    auto split = splitDynState(data);
    return packSize(split.first, comm) + packSize(split.second, comm);
}

template<class T>
std::size_t packSize(const DynamicVector<T>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.data(), comm);
}

std::size_t packSize(const char* str, Dune::MPIHelper::MPICommunicator comm)
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

std::size_t packSize(const std::string& str, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(str.c_str(), comm);
}

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

template<class T, std::size_t N>
std::size_t packSize(const std::array<T,N>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return N*packSize(data[0], comm);
}

HANDLE_AS_POD(Actdims)
HANDLE_AS_POD(Aqudims)
HANDLE_AS_POD(BCConfig::BCFace)
HANDLE_AS_POD(data::Connection)
HANDLE_AS_POD(data::CurrentControl)
HANDLE_AS_POD(data::Rates)
HANDLE_AS_POD(data::Segment)
HANDLE_AS_POD(DENSITYRecord)
HANDLE_AS_POD(DenT::entry)
HANDLE_AS_POD(Eqldims)
HANDLE_AS_POD(MLimits)
HANDLE_AS_POD(PlmixparRecord)
HANDLE_AS_POD(PlyvmhRecord)
HANDLE_AS_POD(PVTWRecord)
HANDLE_AS_POD(PVCDORecord)
HANDLE_AS_POD(Regdims)
HANDLE_AS_POD(RockConfig::RockComp)
HANDLE_AS_POD(ROCKRecord)
HANDLE_AS_POD(SatFuncControls)
HANDLE_AS_POD(ShrateRecord)
HANDLE_AS_POD(StandardCond)
HANDLE_AS_POD(Stone1exRecord)
HANDLE_AS_POD(Tabdims)
HANDLE_AS_POD(TimeStampUTC::YMD)
HANDLE_AS_POD(TlmixparRecord)
HANDLE_AS_POD(Tuning)
HANDLE_AS_POD(VISCREFRecord)
HANDLE_AS_POD(WATDENTRecord)
HANDLE_AS_POD(WellBrineProperties)
HANDLE_AS_POD(Welldims)
HANDLE_AS_POD(WellFoamProperties)
HANDLE_AS_POD(WellSegmentDims)

std::size_t packSize(const data::Well& data, Dune::MPIHelper::MPICommunicator comm)
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

std::size_t packSize(const data::CellData& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.dim, comm) + packSize(data.data, comm) + packSize(data.target, comm);
}

std::size_t packSize(const RestartKey& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.key, comm) + packSize(data.dim, comm) + packSize(data.required, comm);
}

std::size_t packSize(const data::Solution& data, Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    return packSize(static_cast<const std::map< std::string, data::CellData>&>(data), comm);
}

std::size_t packSize(const data::WellRates& data, Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    return packSize(static_cast<const std::map< std::string, data::Well>&>(data), comm);
}

std::size_t packSize(const RestartValue& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.solution, comm) + packSize(data.wells, comm) + packSize(data.extra, comm);
}

std::size_t packSize(const ThresholdPressure& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.active(), comm) +
          packSize(data.restart(), comm) +
          packSize(data.thresholdPressureTable(), comm) +
          packSize(data.pressureTable(), comm);
}

std::size_t packSize(const WellType& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.producer(), comm) +
           packSize(data.preferred_phase(), comm);
}

std::size_t packSize(const DenT& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.records(), comm);
}

std::size_t packSize(const Aquifetp& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.data(), comm);
}

std::size_t packSize(const Aquifetp::AQUFETP_data& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.aquiferID, comm) +
           packSize(data.pvttableID, comm) +
           packSize(data.J, comm) +
           packSize(data.C_t, comm) +
           packSize(data.V0, comm) +
           packSize(data.d0, comm) +
           packSize(data.p0, comm);
}

std::size_t packSize(const AquiferCT::AQUCT_data& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.aquiferID, comm) +
           packSize(data.inftableID, comm) +
           packSize(data.pvttableID, comm) +
           packSize(data.phi_aq, comm) +
           packSize(data.d0, comm) +
           packSize(data.C_t, comm) +
           packSize(data.r_o, comm) +
           packSize(data.k_a, comm) +
           packSize(data.c1, comm) +
           packSize(data.h, comm) +
           packSize(data.theta, comm) +
           packSize(data.c2, comm) +
           packSize(data.p0, comm) +
           packSize(data.td, comm) +
           packSize(data.pi, comm) +
           packSize(data.cell_id, comm);
}

std::size_t packSize(const AquiferCT& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.data(), comm);
}

std::size_t packSize(const AquiferConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.fetp(), comm) +
           packSize(data.ct(), comm) +
           packSize(data.connections(), comm);
}

std::size_t packSize(const Aquancon::AquancCell& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.aquiferID, comm) +
           packSize(data.global_index, comm) +
           packSize(data.influx_coeff, comm) +
           packSize(data.influx_mult, comm) +
           packSize(data.face_dir, comm);
}

std::size_t packSize(const Aquancon& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.data(), comm);
}

std::size_t packSize(const BCConfig& bc, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(bc.faces(), comm);
}

std::size_t packSize(const RockConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.active(), comm) +
           packSize(data.rocknum_property(), comm) +
           packSize(data.comp(), comm) +
           packSize(data.num_rock_tables(), comm) +
           packSize(data.water_compaction(), comm) +
           packSize(data.hysteresis_mode(), comm);
}

std::size_t packSize(const NNC& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.data(), comm);
}

std::size_t packSize(const EDITNNC& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.data(), comm);
}

std::size_t packSize(const Rock2dTable& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.pvmultValues(), comm) +
          packSize(data.pressureValues(), comm);
}

std::size_t packSize(const Rock2dtrTable& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.transMultValues(), comm) +
          packSize(data.pressureValues(), comm);
}

std::size_t packSize(const ColumnSchema& data, Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t res = packSize(data.name(), comm) +
                      packSize(data.order(), comm) +
                      packSize(data.getDefaultMode(), comm);
    if (data.getDefaultMode() == Table::DEFAULT_CONST) {
        res += packSize(data.getDefaultValue(), comm);
    }

    return res;
}

std::size_t packSize(const TableSchema& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.getColumns(), comm);
}

std::size_t packSize(const TableColumn& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.schema(), comm) +
          packSize(data.name(), comm) +
          packSize(data.values(), comm) +
          packSize(data.defaults(), comm) +
          packSize(data.defaultCount(), comm);
}

std::size_t packSize(const SimpleTable& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.schema(), comm) +
          packSize(data.columns(), comm) +
          packSize(data.jfunc(), comm);
}

std::size_t packSize(const TableContainer& data, Dune::MPIHelper::MPICommunicator comm)
{
    size_t res = 2*packSize(data.max(), comm);
    for (const auto& it : data.tables()) {
        if (it.second) {
            res += packSize(it.first, comm) + packSize(*it.second, comm);
        }
    }

    return res;
}

std::size_t packSize(const Equil& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.records(), comm);
}

std::size_t packSize(const FoamConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.records(), comm) +
           packSize(data.getTransportPhase(), comm) +
           packSize(data.getMobilityModel(), comm);
}

std::size_t packSize(const InitConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getEquil(), comm) +
           packSize(data.getFoamConfig(), comm) +
           packSize(data.filleps(), comm) +
           packSize(data.hasGravity(), comm) +
           packSize(data.restartRequested(), comm) +
           packSize(data.getRestartStep(), comm) +
           packSize(data.getRestartRootName(), comm);
}

std::size_t packSize(const SimulationConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getThresholdPressure(), comm) +
           packSize(data.bcconfig(), comm) +
           packSize(data.rock_config(), comm) +
           packSize(data.useCPR(), comm) +
           packSize(data.hasDISGAS(), comm) +
           packSize(data.hasVAPOIL(), comm) +
           packSize(data.isThermal(), comm);
}

std::size_t packSize(const TimeMap& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timeList(), comm);
}

std::size_t packSize(const RestartConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timeMap(), comm) +
           packSize(data.getFirstRestartStep(), comm) +
           packSize(data.writeInitialRst(), comm) +
           packSize(data.restartSchedule(), comm) +
           packSize(data.restartKeywords(), comm) +
           packSize(data.saveKeywords(), comm);
}

std::size_t packSize(const IOConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getWriteINITFile(), comm) +
           packSize(data.getWriteEGRIDFile(), comm) +
           packSize(data.getUNIFIN(), comm) +
           packSize(data.getUNIFOUT(), comm) +
           packSize(data.getFMTIN(), comm) +
           packSize(data.getFMTOUT(), comm) +
           packSize(data.getDeckFileName(), comm) +
           packSize(data.getOutputEnabled(), comm) +
           packSize(data.getOutputDir(), comm) +
           packSize(data.getNoSim(), comm) +
           packSize(data.getBaseName(), comm) +
           packSize(data.getEclCompatibleRST(), comm);
}

std::size_t packSize(const Phases& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getBits(), comm);
}

std::size_t packSize(const EndpointScaling& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getBits(), comm);
}

std::size_t packSize(const UDQParams& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.reseed(), comm) +
           packSize(data.rand_seed(), comm) +
           packSize(data.range(), comm) +
           packSize(data.undefinedValue(), comm) +
           packSize(data.cmpEpsilon(), comm);
}

std::size_t packSize(const Runspec& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.phases(), comm) +
           packSize(data.tabdims(), comm) +
           packSize(data.endpointScaling(), comm) +
           packSize(data.wellDimensions(), comm) +
           packSize(data.wellSegmentDimensions(), comm) +
           packSize(data.udqParams(), comm) +
           packSize(data.hysterPar(), comm) +
           packSize(data.actdims(), comm) +
           packSize(data.saturationFunctionControls(), comm);
}

std::size_t packSize(const PvtxTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getOuterColumnSchema(), comm) +
           packSize(data.getOuterColumn(), comm) +
           packSize(data.getUnderSaturatedSchema(), comm) +
           packSize(data.getSaturatedSchema(), comm) +
           packSize(data.getUnderSaturatedTables(), comm) +
           packSize(data.getSaturatedTable(), comm);
}

std::size_t packSize(const PvtgTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PvtxTable&>(data), comm);
}

std::size_t packSize(const PvtoTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PvtxTable&>(data), comm);
}

std::size_t packSize(const PvtwTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<PVTWRecord>&>(data), comm);
}

std::size_t packSize(const PvcdoTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<PVCDORecord>&>(data), comm);
}

std::size_t packSize(const DensityTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<DENSITYRecord>&>(data), comm);
}

std::size_t packSize(const ViscrefTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<VISCREFRecord>&>(data), comm);
}

std::size_t packSize(const WatdentTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<WATDENTRecord>&>(data), comm);
}

std::size_t packSize(const PolyInjTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getThroughputs(), comm) +
           packSize(data.getVelocities(), comm) +
           packSize(data.getTableNumber(), comm) +
           packSize(data.getTableData(), comm);
}

std::size_t packSize(const PlymwinjTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PolyInjTable&>(data), comm);
}

std::size_t packSize(const SkprpolyTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PolyInjTable&>(data), comm) +
           packSize(data.referenceConcentration(), comm);
}

std::size_t packSize(const SkprwatTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PolyInjTable&>(data), comm);
}

std::size_t packSize(const RockTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<ROCKRecord>&>(data), comm);
}

namespace {

struct SplitSimpleTables {
    size_t plyshMax = 0;
    size_t rockMax = 0;
    std::map<size_t, std::shared_ptr<PlyshlogTable>> plyshMap;
    std::map<size_t, std::shared_ptr<RocktabTable>> rockMap;
};

SplitSimpleTables
splitSimpleTable(std::map<std::string, TableContainer>& simpleTables)
{
    SplitSimpleTables result;

    // PlyshlogTable need special treatment
    auto it = simpleTables.find("PLYSHLOG");
    if (it != simpleTables.end()) {
        result.plyshMax = it->second.max();
        for (const auto& mapIt : it->second.tables()) {
            auto ptr = std::static_pointer_cast<PlyshlogTable>(mapIt.second);
            result.plyshMap.insert(std::make_pair(mapIt.first, ptr));
        }
        simpleTables.erase(it);
    }

    // RocktabTable need special treatment
    it = simpleTables.find("ROCKMAP");
    if (it != simpleTables.end()) {
        result.rockMax = it->second.max();
        for (const auto& mapIt : it->second.tables()) {
            auto ptr = std::static_pointer_cast<RocktabTable>(mapIt.second);
            result.rockMap.insert(std::make_pair(mapIt.first,  ptr));
        }
        simpleTables.erase(it);
    }

    return result;
}

}

std::size_t packSize(const TableManager& data, Dune::MPIHelper::MPICommunicator comm)
{
    auto simpleTables = data.getSimpleTables();
    auto splitTab = splitSimpleTable(simpleTables);

    return packSize(simpleTables, comm) +
           packSize(splitTab.plyshMax, comm) +
           packSize(splitTab.plyshMap, comm) +
           packSize(splitTab.rockMax, comm) +
           packSize(splitTab.rockMap, comm) +
           packSize(data.getPvtgTables(), comm) +
           packSize(data.getPvtoTables(), comm) +
           packSize(data.getRock2dTables(), comm) +
           packSize(data.getRock2dtrTables(), comm) +
           packSize(data.getPvtwTable(), comm) +
           packSize(data.getPvcdoTable(), comm) +
           packSize(data.getDensityTable(), comm) +
           packSize(data.getPlyvmhTable(), comm) +
           packSize(data.getRockTable(), comm) +
           packSize(data.getPlmixparTable(), comm) +
           packSize(data.getShrateTable(), comm) +
           packSize(data.getStone1exTable(), comm) +
           packSize(data.getTlmixparTable(), comm) +
           packSize(data.getViscrefTable(), comm) +
           packSize(data.getWatdentTable(), comm) +
           packSize(data.getPvtwSaltTables(), comm) +
           packSize(data.getBrineDensityTables(), comm) +
           packSize(data.getSolventDensityTables(), comm) +
           packSize(data.getPlymwinjTables(), comm) +
           packSize(data.getSkprwatTables(), comm) +
           packSize(data.getSkprpolyTables(), comm) +
           packSize(data.getTabdims(), comm) +
           packSize(data.getRegdims(), comm) +
           packSize(data.getEqldims(), comm) +
           packSize(data.getAqudims(), comm) +
           packSize(data.useImptvd(), comm) +
           packSize(data.useEnptvd(), comm) +
           packSize(data.useEqlnum(), comm) +
           packSize(data.useShrate(), comm) +
           packSize(data.useJFunc(), comm) +
          (data.useJFunc() ? packSize(data.getJFunc(), comm) : 0) +
           packSize(data.OilDenT(), comm) +
           packSize(data.GasDenT(), comm) +
           packSize(data.WatDenT(), comm) +
           packSize(data.stCond(), comm) +
           packSize(data.gas_comp_index(), comm) +
           packSize(data.rtemp(), comm);
}

template
std::size_t packSize(const std::map<Phase,Group::GroupInjectionProperties>& data,
                     Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const OilVaporizationProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getType(), comm) +
           packSize(data.vap1(), comm) +
           packSize(data.vap2(), comm) +
           packSize(data.maxDRSDT(), comm) +
           packSize(data.maxDRSDT_allCells(), comm) +
           packSize(data.maxDRVDT(), comm);
}

std::size_t packSize(const Events& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.events(), comm);
}

std::size_t packSize(const MessageLimits& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getLimits(), comm);
}

std::size_t packSize(const VFPInjTable& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getTableNum(), comm) +
           packSize(data.getDatumDepth(), comm) +
           packSize(data.getFloType(), comm) +
           packSize(data.getFloAxis(), comm) +
           packSize(data.getTHPAxis(), comm) +
           packSize(data.getTable(), comm);
}

std::size_t packSize(const VFPProdTable& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getTableNum(), comm) +
           packSize(data.getDatumDepth(), comm) +
           packSize(data.getFloType(), comm) +
           packSize(data.getWFRType(), comm) +
           packSize(data.getGFRType(), comm) +
           packSize(data.getALQType(), comm) +
           packSize(data.getFloAxis(), comm) +
           packSize(data.getTHPAxis(), comm) +
           packSize(data.getWFRAxis(), comm) +
           packSize(data.getGFRAxis(), comm) +
           packSize(data.getALQAxis(), comm) +
           packSize(data.getTable(), comm);
}

std::size_t packSize(const WellTestConfig::WTESTWell& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name, comm) +
           packSize(data.shut_reason, comm) +
           packSize(data.test_interval, comm) +
           packSize(data.num_test, comm) +
           packSize(data.startup_time, comm) +
           packSize(data.begin_report_step, comm);
}

std::size_t packSize(const WellTestConfig& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getWells(), comm);
}

std::size_t packSize(const WellTracerProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getConcentrations(), comm);
}

std::size_t packSize(const UDAValue& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.get_dim(), comm) +
           packSize(data.is<double>(), comm) +
           (data.is<double>() ? packSize(data.get<double>(), comm) :
                                packSize(data.get<std::string>(), comm));
}

std::size_t packSize(const Connection& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.dir(), comm) +
           packSize(data.depth(), comm) +
           packSize(data.state(), comm) +
           packSize(data.satTableId(), comm) +
           packSize(data.complnum(), comm) +
           packSize(data.CF(), comm) +
           packSize(data.Kh(), comm) +
           packSize(data.rw(), comm) +
           packSize(data.r0(), comm) +
           packSize(data.skinFactor(), comm) +
           packSize(data.getI(), comm) +
           packSize(data.getJ(), comm) +
           packSize(data.getK(), comm) +
           packSize(data.kind(), comm) +
           packSize(data.getSeqIndex(), comm) +
           packSize(data.getSegDistStart(), comm) +
           packSize(data.getSegDistEnd(), comm) +
           packSize(data.getDefaultSatTabId(), comm) +
           packSize(data.getCompSegSeqIndex(), comm) +
           packSize(data.segment(), comm) +
           packSize(data.wellPi(), comm);
}

std::size_t packSize(const Well::WellInjectionProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name, comm) +
           packSize(data.surfaceInjectionRate, comm) +
           packSize(data.reservoirInjectionRate, comm) +
           packSize(data.BHPTarget, comm) +
           packSize(data.THPTarget, comm) +
           packSize(data.bhp_hist_limit, comm) +
           packSize(data.thp_hist_limit, comm) +
           packSize(data.temperature, comm) +
           packSize(data.BHPH, comm) +
           packSize(data.THPH, comm) +
           packSize(data.VFPTableNumber, comm) +
           packSize(data.predictionMode, comm) +
           packSize(data.injectionControls, comm) +
           packSize(data.injectorType, comm) +
           packSize(data.controlMode, comm);
}

std::size_t packSize(const WellEconProductionLimits& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.minOilRate(), comm) +
           packSize(data.minGasRate(), comm) +
           packSize(data.maxWaterCut(), comm) +
           packSize(data.maxGasOilRatio(), comm) +
           packSize(data.maxWaterGasRatio(), comm) +
           packSize(data.workover(), comm) +
           packSize(data.endRun(), comm) +
           packSize(data.followonWell(), comm) +
           packSize(data.quantityLimit(), comm) +
           packSize(data.maxSecondaryMaxWaterCut(), comm) +
           packSize(data.workoverSecondary(), comm) +
           packSize(data.maxGasLiquidRatio(), comm) +
           packSize(data.minLiquidRate(), comm) +
           packSize(data.maxTemperature(), comm) +
           packSize(data.minReservoirFluidRate(), comm);
}

std::size_t packSize(const WellConnections& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getHeadI(), comm) +
           packSize(data.getHeadJ(), comm) +
           packSize(data.getNumRemoved(), comm) +
           packSize(data.getConnections(), comm);
}

std::size_t packSize(const Well::WellProductionProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name, comm) +
           packSize(data.OilRate, comm) +
           packSize(data.WaterRate, comm) +
           packSize(data.GasRate, comm) +
           packSize(data.LiquidRate, comm) +
           packSize(data.ResVRate, comm) +
           packSize(data.BHPTarget, comm) +
           packSize(data.THPTarget, comm) +
           packSize(data.bhp_hist_limit, comm) +
           packSize(data.thp_hist_limit, comm) +
           packSize(data.BHPH, comm) +
           packSize(data.THPH, comm) +
           packSize(data.VFPTableNumber, comm) +
           packSize(data.ALQValue, comm) +
           packSize(data.predictionMode, comm) +
           packSize(data.controlMode, comm) +
           packSize(data.whistctl_cmode, comm) +
           packSize(data.getNumProductionControls(), comm);
}

std::size_t packSize(const SpiralICD& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.strength(), comm) +
           packSize(data.length(), comm) +
           packSize(data.densityCalibration(), comm) +
           packSize(data.viscosityCalibration(), comm) +
           packSize(data.criticalValue(), comm) +
           packSize(data.widthTransitionRegion(), comm) +
           packSize(data.maxViscosityRatio(), comm) +
           packSize(data.methodFlowScaling(), comm) +
           packSize(data.maxAbsoluteRate(), comm) +
           packSize(data.status(), comm) +
           packSize(data.scalingFactor(), comm);
}

std::size_t packSize(const Valve& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.conEFlowCoefficient(), comm) +
           packSize(data.conCrossArea(), comm) +
           packSize(data.conMaxCrossArea(), comm) +
           packSize(data.pipeAdditionalLength(), comm) +
           packSize(data.pipeDiameter(), comm) +
           packSize(data.pipeRoughness(), comm) +
           packSize(data.pipeCrossArea(), comm) +
           packSize(data.status(), comm);
}

std::size_t packSize(const Segment& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.segmentNumber(), comm) +
           packSize(data.branchNumber(), comm) +
           packSize(data.outletSegment(), comm) +
           packSize(data.inletSegments(), comm) +
           packSize(data.totalLength(), comm) +
           packSize(data.depth(), comm) +
           packSize(data.internalDiameter(), comm) +
           packSize(data.roughness(), comm) +
           packSize(data.crossArea(), comm) +
           packSize(data.volume(), comm) +
           packSize(data.dataReady(), comm) +
           packSize(data.segmentType(), comm) +
           packSize(data.spiralICD(), comm) +
           packSize(data.getValve(), comm);
}

template<class T>
std::size_t packSize(const std::shared_ptr<T>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(bool(), comm);
    if (data)
         size += packSize(*data, comm);

    return size;
}

template<class T>
std::size_t packSize(const std::unique_ptr<T>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(bool(), comm);
    if (data)
         size += packSize(*data, comm);

    return size;
}

std::size_t packSize(const Dimension& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getSIScalingRaw(), comm) +
           packSize(data.getSIOffset(), comm);
}

std::size_t packSize(const UnitSystem& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getName(), comm) +
           packSize(data.getType(), comm) +
           packSize(data.getDimensions(), comm) +
           packSize(data.use_count(), comm);
}

std::size_t packSize(const WellSegments& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.compPressureDrop(), comm) +
           packSize(data.segments(), comm);
}

std::size_t packSize(const Well& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.name(), comm) +
                       packSize(data.groupName(), comm) +
                       packSize(data.firstTimeStep(), comm) +
                       packSize(data.seqIndex(), comm) +
                       packSize(data.getHeadI(), comm) +
                       packSize(data.getHeadJ(), comm) +
                       packSize(data.getRefDepth(), comm) +
                       packSize(data.wellType(), comm) +
                       packSize(data.getWellConnectionOrdering(), comm) +
                       packSize(data.units(), comm) +
                       packSize(data.udqUndefined(), comm) +
                       packSize(data.getStatus(), comm) +
                       packSize(data.getDrainageRadius(), comm) +
                       packSize(data.getAllowCrossFlow(), comm) +
                       packSize(data.getAutomaticShutIn(), comm) +
                       packSize(data.wellGuideRate(), comm) +
                       packSize(data.getEfficiencyFactor(), comm) +
                       packSize(data.getSolventFraction(), comm) +
                       packSize(data.predictionMode(), comm) +
                       packSize(data.getEconLimits(), comm) +
                       packSize(data.getFoamProperties(), comm) +
                       packSize(data.getPolymerProperties(), comm) +
                       packSize(data.getBrineProperties(), comm) +
                       packSize(data.getTracerProperties(), comm) +
                       packSize(data.getConnections(), comm) +
                       packSize(data.getProductionProperties(), comm) +
                       packSize(data.getInjectionProperties(), comm) +
                       packSize(data.hasSegments(), comm);
    if (data.hasSegments())
        size += packSize(data.getSegments(), comm);

    return size;
}

template<class T>
std::size_t packSize(const IOrderSet<T>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.index(), comm) +
           packSize(data.data(), comm);
}

std::size_t packSize(const Group::GroupInjectionProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.phase, comm) +
           packSize(data.cmode, comm) +
           packSize(data.surface_max_rate, comm) +
           packSize(data.resv_max_rate, comm) +
           packSize(data.target_reinj_fraction, comm) +
           packSize(data.target_void_fraction, comm) +
           packSize(data.reinj_group, comm) +
           packSize(data.voidage_group, comm) +
           packSize(data.injection_controls, comm);
}

std::size_t packSize(const Group::GroupProductionProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.cmode, comm) +
           packSize(data.exceed_action, comm) +
           packSize(data.oil_target, comm) +
           packSize(data.water_target, comm) +
           packSize(data.gas_target, comm) +
           packSize(data.liquid_target, comm) +
           packSize(data.guide_rate, comm) +
           packSize(data.guide_rate_def, comm) +
           packSize(data.resv_target, comm) +
           packSize(data.production_controls, comm);
}

std::size_t packSize(const Group& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name(), comm) +
           packSize(data.insert_index(), comm) +
           packSize(data.initStep(), comm) +
           packSize(data.udqUndefined(), comm) +
           packSize(data.units(), comm) +
           packSize(data.type(), comm) +
           packSize(data.getGroupEfficiencyFactor(), comm) +
           packSize(data.getTransferGroupEfficiencyFactor(), comm) +
           packSize(data.isAvailableForGroupControl(), comm) +
           packSize(data.getGroupNetVFPTable(), comm) +
           packSize(data.parent(), comm) +
           packSize(data.iwells(), comm) +
           packSize(data.igroups(), comm) +
           packSize(data.injectionProperties(), comm) +
           packSize(data.productionProperties(), comm);
}

std::size_t packSize(const WList& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.wellList(), comm);
}

std::size_t packSize(const WListManager& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.lists(), comm);
}

std::size_t packSize(const UDQASTNode& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.var_type, comm) +
           packSize(data.getType(), comm) +
           packSize(data.stringValue(), comm) +
           packSize(data.scalarValue(), comm) +
           packSize(data.getSelectors(), comm) +
           packSize(data.getLeft(), comm) +
           packSize(data.getRight(), comm);
}

std::size_t packSize(const UDQDefine& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.keyword(), comm) +
           packSize(data.getAst(), comm) +
           packSize(data.var_type(), comm) +
           packSize(data.input_string(), comm);
}

std::size_t packSize(const UDQAssign::AssignRecord& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.selector, comm) +
           packSize(data.value, comm);
}

std::size_t packSize(const UDQAssign& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.keyword(), comm) +
           packSize(data.var_type(), comm) +
           packSize(data.getRecords(), comm);
}

std::size_t packSize(const UDQIndex& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.insert_index, comm) +
           packSize(data.typed_insert_index, comm) +
           packSize(data.action, comm) +
           packSize(data.var_type, comm);
}

std::size_t packSize(const UDQConfig& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.params(), comm) +
           packSize(data.definitionMap(), comm) +
           packSize(data.assignmentMap(), comm) +
           packSize(data.unitsMap(), comm) +
           packSize(data.inputIndex(), comm) +
           packSize(data.typeCount(), comm);
}

std::size_t packSize(const UDQActive::InputRecord& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.input_index, comm) +
           packSize(data.udq, comm) +
           packSize(data.wgname, comm) +
           packSize(data.control, comm);
}

std::size_t packSize(const UDQActive::Record& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.udq, comm) +
           packSize(data.input_index, comm) +
           packSize(data.use_index, comm) +
           packSize(data.wgname, comm) +
           packSize(data.control, comm) +
           packSize(data.uad_code, comm) +
           packSize(data.use_count, comm);
}

std::size_t packSize(const UDQActive& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getInputRecords(), comm) +
           packSize(data.getOutputRecords(), comm) +
           packSize(data.getUdqKeys(), comm) +
           packSize(data.getWgKeys(), comm);
}

std::size_t packSize(const GuideRateModel& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timeInterval(), comm) +
           packSize(data.target(), comm) +
           packSize(data.coefs(), comm) +
           packSize(data.allow_increase(), comm) +
           packSize(data.damping_factor(), comm) +
           packSize(data.free_gas(), comm) +
           packSize(data.defaultModel(), comm) +
           packSize(data.udaCoefs(), comm);
}

std::size_t packSize(const GuideRateConfig& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getModel(), comm) +
           packSize(data.getWells(), comm) +
           packSize(data.getGroups(), comm);
}

std::size_t packSize(const GConSale::GCONSALEGroup& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.sales_target, comm) +
           packSize(data.max_sales_rate, comm) +
           packSize(data.min_sales_rate, comm) +
           packSize(data.max_proc, comm) +
           packSize(data.udq_undefined, comm) +
           packSize(data.unit_system, comm);
}

std::size_t packSize(const GConSale& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getGroups(), comm);
}

std::size_t packSize(const GConSump::GCONSUMPGroup& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.consumption_rate, comm) +
           packSize(data.import_rate, comm) +
           packSize(data.network_node, comm) +
           packSize(data.udq_undefined, comm) +
           packSize(data.unit_system, comm);
}

std::size_t packSize(const GConSump& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getGroups(), comm);
}

std::size_t packSize(const RFTConfig& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timeMap(), comm) +
           packSize(data.firstRFTOutput(), comm) +
           packSize(data.wellOpenRftTime(), comm) +
           packSize(data.wellOpenRftName(), comm) +
           packSize(data.wellOpen(), comm) +
           packSize(data.rftConfig(), comm) +
           packSize(data.pltConfig(), comm);
}

std::size_t packSize(const DeckItem& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.dVal(), comm) +
           packSize(data.iVal(), comm) +
           packSize(data.sVal(), comm) +
           packSize(data.uVal(), comm) +
           packSize(data.getType(), comm) +
           packSize(data.name(), comm) +
           packSize(data.valueStatus(), comm) +
           packSize(data.rawData(), comm) +
           packSize(data.activeDimensions(), comm) +
           packSize(data.defaultDimensions(), comm);
}

std::size_t packSize(const DeckRecord& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getItems(), comm);
}

std::size_t packSize(const Location& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.filename, comm) +
           packSize(data.lineno, comm);
}

std::size_t packSize(const DeckKeyword& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name(), comm) +
           packSize(data.location(), comm) +
           packSize(data.records(), comm) +
           packSize(data.isDataKeyword(), comm) +
           packSize(data.isSlashTerminated(), comm);
}

std::size_t packSize(const Deck& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.keywords(), comm) +
           packSize(data.getDefaultUnitSystem(), comm) +
           packSize(data.activeUnitSystem(), comm) +
           packSize(data.getDataFile(), comm) +
           packSize(data.getInputPath(), comm) +
           packSize(data.unitSystemAccessCount(), comm);
}

std::size_t packSize(const Action::ASTNode& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.type, comm) +
           packSize(data.func_type, comm) +
           packSize(data.func, comm) +
           packSize(data.argList(), comm) +
           packSize(data.getNumber(), comm) +
           packSize(data.childrens(), comm);
}

std::size_t packSize(const Action::AST& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getCondition(), comm);
}

std::size_t packSize(const Action::Quantity& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.quantity, comm) +
           packSize(data.args, comm);
}

std::size_t packSize(const Action::Condition& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.lhs, comm) +
           packSize(data.rhs, comm) +
           packSize(data.logic, comm) +
           packSize(data.cmp, comm) +
           packSize(data.cmp_string, comm);
}

std::size_t packSize(const Action::ActionX& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name(), comm) +
           packSize(data.max_run(), comm) +
           packSize(data.min_wait(), comm) +
           packSize(data.start_time(), comm) +
           packSize(data.getKeywords(), comm) +
           packSize(data.getCondition(), comm) +
           packSize(data.conditions(), comm) +
           packSize(data.getRunCount(), comm) +
           packSize(data.getLastRun(), comm);
}

std::size_t packSize(const Action::Actions& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getActions(), comm);
}

template<class Key, class T> using Map2 = std::map<Key,T>;

std::size_t packSize(const Schedule& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getTimeMap(), comm) +
           packSizeDynMap(data.getStaticWells(), comm) +
           packSizeDynMap(data.getGroups(), comm) +
           packSize(data.getOilVapProps(), comm) +
           packSize(data.getEvents(), comm) +
           packSize(data.getModifierDeck(), comm) +
           packSize(data.getTuning(), comm) +
           packSize(data.getMessageLimits(), comm) +
           packSize(data.getRunspec(), comm) +
           packSizeDynMap<Map2>(data.getVFPProdTables(), comm) +
           packSizeDynMap<Map2>(data.getVFPInjTables(), comm) +
           packSize(data.getWellTestConfig(), comm) +
           packSize(data.getWListManager(), comm) +
           packSize(data.getUDQConfig(), comm) +
           packSize(data.getUDQActive(), comm) +
           packSize(data.getGuideRateConfig(), comm) +
           packSize(data.getGConSale(), comm) +
           packSize(data.getGConSump(), comm) +
           packSize(data.getGlobalWhistCtlMode(), comm) +
           packSize(data.getActions(), comm) +
           packSize(data.rftConfig(), comm) +
           packSize(data.getNupCol(), comm) +
           packSize(data.restart(), comm) +
           packSize(data.getWellGroupEvents(), comm);
}

std::size_t packSize(const BrineDensityTable& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getBrineDensityColumn(), comm);
}

std::size_t packSize(const PvtwsaltTable& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.getReferencePressureValue(), comm) +
          packSize(data.getReferenceSaltConcentrationValue(), comm) +
          packSize(data.getTableValues(), comm);
}

std::size_t packSize(const EquilRecord& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.datumDepth(), comm) +
           packSize(data.datumDepthPressure(), comm) +
           packSize(data.waterOilContactDepth(), comm) +
           packSize(data.waterOilContactCapillaryPressure(), comm) +
           packSize(data.gasOilContactDepth(), comm) +
           packSize(data.gasOilContactCapillaryPressure(), comm) +
           packSize(data.liveOilInitConstantRs(), comm) +
           packSize(data.wetGasInitConstantRv(), comm) +
           packSize(data.initializationTargetAccuracy(), comm);
}

std::size_t packSize(const FoamData& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.referenceSurfactantConcentration(), comm) +
           packSize(data.exponent(), comm) +
           packSize(data.minimumSurfactantConcentration(), comm) +
           packSize(data.allowDesorption(), comm) +
           packSize(data.rockDensity(), comm);
}

std::size_t packSize(const RestartSchedule& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timestep, comm) +
           packSize(data.basic, comm) +
           packSize(data.frequency, comm) +
           packSize(data.rptsched_restart_set, comm) +
           packSize(data.rptsched_restart, comm);
}

std::size_t packSize(const TimeStampUTC& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.ymd(), comm) +
           packSize(data.hour(), comm) +
           packSize(data.minutes(), comm) +
           packSize(data.seconds(), comm) +
           packSize(data.microseconds(), comm);
}

std::size_t packSize(const EclHysterConfig& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.active(), comm) +
           packSize(data.pcHysteresisModel(), comm) +
           packSize(data.krHysteresisModel(), comm);
}

std::size_t packSize(const JFunc& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.flag(), comm) +
           packSize(data.owSurfaceTension(), comm) +
           packSize(data.goSurfaceTension(), comm) +
           packSize(data.alphaFactor(), comm) +
           packSize(data.betaFactor(), comm) +
           packSize(data.direction(), comm);
}

std::size_t packSize(const WellPolymerProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.m_polymerConcentration, comm) +
           packSize(data.m_saltConcentration, comm) +
           packSize(data.m_plymwinjtable, comm) +
           packSize(data.m_skprwattable, comm) +
           packSize(data.m_skprpolytable, comm);
}

std::size_t packSize(const Well::WellGuideRate& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.available, comm) +
           packSize(data.guide_rate, comm) +
           packSize(data.guide_phase, comm) +
           packSize(data.scale_factor, comm);
}

std::size_t packSize(const GuideRateConfig::WellTarget& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.guide_rate, comm) +
           packSize(data.target, comm) +
           packSize(data.scaling_factor, comm);
}

std::size_t packSize(const GuideRateConfig::GroupTarget& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.guide_rate, comm) +
           packSize(data.target, comm);
}

std::size_t packSize(const MULTREGTRecord& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.src_value, comm) +
           packSize(data.target_value, comm) +
           packSize(data.trans_mult, comm) +
           packSize(data.directions, comm) +
           packSize(data.nnc_behaviour, comm) +
           packSize(data.region_name, comm);
}

std::size_t packSize(const MULTREGTScanner& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getSize(), comm) +
           packSize(data.getRecords(), comm) +
           packSize(data.getSearchMap(), comm) +
           packSize(data.getRegions(), comm) +
           packSize(data.getDefaultRegion(), comm);
}

std::size_t packSize(const EclipseConfig& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.init(), comm) +
           packSize(data.io(), comm);
}

std::size_t packSize(const TransMult& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getSize(), comm) +
           packSize(data.getTrans(), comm) +
           packSize(data.getNames(), comm) +
           packSize(data.getScanner(), comm);
}

std::size_t packSize(const FaultFace& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getIndices(), comm) +
           packSize(data.getDir(), comm);
}

std::size_t packSize(const Fault& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getName(), comm) +
           packSize(data.getTransMult(), comm) +
           packSize(data.getFaceList(), comm);
}

std::size_t packSize(const FaultCollection& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getFaults(), comm);
}

std::size_t packSize(const SolventDensityTable& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getSolventDensityColumn(), comm);
}

std::size_t packSize(const GridDims& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getNXYZ(), comm);
}

std::size_t packSize(const ShrateTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<ShrateRecord>&>(data), comm);
}

std::size_t packSize(const TlmixparTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<TlmixparRecord>&>(data), comm);
}

std::size_t packSize(const PlmixparTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<PlmixparRecord>&>(data), comm);
}

std::size_t packSize(const PlyvmhTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<PlyvmhRecord>&>(data), comm);
}

std::size_t packSize(const Stone1exTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<Stone1exRecord>&>(data), comm);
}

std::size_t packSize(const PlyshlogTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const SimpleTable&>(data), comm) +
           packSize(data.getRefPolymerConcentration(), comm) +
           packSize(data.getRefSalinity(), comm) +
           packSize(data.getRefTemperature(), comm) +
           packSize(data.hasRefSalinity(), comm) +
           packSize(data.hasRefTemperature(), comm);
}

std::size_t packSize(const RocktabTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const SimpleTable&>(data), comm) +
           packSize(data.isDirectional(), comm);
}

////// pack routines

template<class T>
void pack(const T*, std::size_t, std::vector<char>&, int&,
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

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

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data, l, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.first, buffer, position, comm);
    pack(data.second, buffer, position, comm);
}

template<class T, class A>
void pack(const std::vector<T, A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    if (std::is_pod<T>::value)
    {
        // size written automatically
        pack(data.data(), data.size(), buffer, position, comm);
        return;
    }

    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
        pack(entry, buffer, position, comm);
}

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

template<class T, size_t N>
void pack(const std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    for (const T& entry : data)
        pack(entry, buffer, position, comm);
}

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

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, void>::type
pack_tuple_entry(const Tuple&, std::vector<char>&, int&,
                      Dune::MPIHelper::MPICommunicator)
{
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, void>::type
pack_tuple_entry(const Tuple& tuple, std::vector<char>& buffer,
                 int& position, Dune::MPIHelper::MPICommunicator comm)
{
    pack(std::get<I>(tuple), buffer, position, comm);
    pack_tuple_entry<I+1>(tuple, buffer, position, comm);
}

template<class... Ts>
void pack(const std::tuple<Ts...>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm)
{
    pack_tuple_entry(data, buffer, position, comm);
}

template<class Key, class Value>
void pack(const OrderedMap<Key, Value>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getIndex(), buffer, position, comm);
    pack(data.getStorage(), buffer, position, comm);
}

template<class T>
void pack(const DynamicState<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    auto split = splitDynState(data);
    pack(split.first, buffer, position, comm);
    pack(split.second, buffer, position, comm);
}

template<class T>
void pack(const DynamicVector<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.data(), buffer, position, comm);
}

void pack(const char* str, std::vector<char>& buffer, int& position,
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

void pack(const std::string& str, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(str.c_str(), buffer, position, comm);
}

template<class T1, class T2, class C, class A>
void pack(const std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
    {
        pack(entry, buffer, position, comm);
    }
}

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

template void pack(const std::map<Phase, Group::GroupInjectionProperties>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

void pack(const data::Well& data, std::vector<char>& buffer, int& position,
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

void pack(const RestartKey& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.key, buffer, position, comm);
    pack(data.dim, buffer, position, comm);
    pack(data.required, buffer, position, comm);
}

void pack(const data::CellData& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.dim, buffer, position, comm);
    pack(data.data, buffer, position, comm);
    pack(data.target, buffer, position, comm);
}

void pack(const data::Solution& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    pack(static_cast<const std::map< std::string, data::CellData>&>(data),
         buffer, position, comm);
}

void pack(const data::WellRates& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    pack(static_cast<const std::map< std::string, data::Well>&>(data),
         buffer, position, comm);
}

void pack(const RestartValue& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.solution, buffer, position, comm);
    pack(data.wells, buffer, position, comm);
    pack(data.extra, buffer, position, comm);
}

void pack(const ThresholdPressure& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.active(), buffer, position, comm);
    pack(data.restart(), buffer, position, comm);
    pack(data.thresholdPressureTable(), buffer, position, comm);
    pack(data.pressureTable(), buffer, position, comm);
}

void pack(const AquiferCT::AQUCT_data& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.aquiferID, buffer, position, comm);
    pack(data.inftableID, buffer, position, comm);
    pack(data.pvttableID, buffer, position, comm);
    pack(data.phi_aq, buffer, position, comm);
    pack(data.d0, buffer, position, comm);
    pack(data.C_t, buffer, position, comm);
    pack(data.r_o, buffer, position, comm);
    pack(data.k_a, buffer, position, comm);
    pack(data.c1, buffer, position, comm);
    pack(data.h, buffer, position, comm);
    pack(data.theta, buffer, position, comm);
    pack(data.c2, buffer, position, comm);
    pack(data.p0, buffer, position, comm);
    pack(data.td, buffer, position, comm);
    pack(data.pi, buffer, position, comm);
    pack(data.cell_id, buffer, position, comm);
}

void pack(const AquiferCT& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.data(), buffer, position, comm);
}

void pack(const Aquifetp::AQUFETP_data& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.aquiferID, buffer, position, comm);
    pack(data.pvttableID, buffer, position, comm);
    pack(data.J, buffer, position, comm);
    pack(data.C_t, buffer, position, comm);
    pack(data.V0, buffer, position, comm);
    pack(data.d0, buffer, position, comm);
    pack(data.p0, buffer, position, comm);
}

void pack(const WellType& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.producer(), buffer, position, comm);
    pack(data.preferred_phase(), buffer, position, comm);
}

void pack(const DenT& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.records(), buffer, position, comm);
}

void pack(const Aquifetp& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.data(), buffer, position, comm);
}

void pack(const AquiferConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.fetp(), buffer, position, comm);
    pack(data.ct(), buffer, position, comm);
    pack(data.connections(), buffer, position, comm);
}

void pack(const Aquancon::AquancCell& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.aquiferID, buffer, position, comm);
    pack(data.global_index, buffer, position, comm);
    pack(data.influx_coeff, buffer, position, comm);
    pack(data.influx_mult, buffer, position, comm);
    pack(data.face_dir, buffer, position, comm);
}

void pack(const Aquancon& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.data(), buffer, position, comm);
}

void pack(const BCConfig& bc, std::vector<char>& buffer, int& position,
    Dune::MPIHelper::MPICommunicator comm)
{
    pack(bc.faces(), buffer, position, comm);
}

void pack(const RockConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.active(), buffer, position, comm);
    pack(data.rocknum_property(), buffer, position, comm);
    pack(data.comp(), buffer, position, comm);
    pack(data.num_rock_tables(), buffer, position, comm);
    pack(data.water_compaction(), buffer, position, comm);
    pack(data.hysteresis_mode(), buffer, position, comm);
}

void pack(const NNC& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.data(), buffer, position, comm);
}

void pack(const EDITNNC& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.data(), buffer, position, comm);
}

void pack(const Rock2dTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.pvmultValues(), buffer, position, comm);
    pack(data.pressureValues(), buffer, position, comm);
}

void pack(const Rock2dtrTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.transMultValues(), buffer, position, comm);
    pack(data.pressureValues(), buffer, position, comm);
}

void pack(const ColumnSchema& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.order(), buffer, position, comm);
    pack(data.getDefaultMode(), buffer, position, comm);
    if (data.getDefaultMode() == Table::DEFAULT_CONST)
        pack(data.getDefaultValue(), buffer, position, comm);
}

void pack(const TableSchema& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getColumns(), buffer, position, comm);
}

void pack(const TableColumn& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.schema(), buffer, position, comm);
    pack(data.name(), buffer, position, comm);
    pack(data.values(), buffer, position, comm);
    pack(data.defaults(), buffer, position, comm);
    pack(data.defaultCount(), buffer, position, comm);
}

void pack(const SimpleTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.schema(), buffer, position, comm);
    pack(data.columns(), buffer, position, comm);
    pack(data.jfunc(), buffer, position, comm);
}

void pack(const TableContainer& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.max(), buffer, position, comm);
    size_t entries = 0;
    for (const auto& it : data.tables()) {
        if (it.second) {
          ++entries;
        }
    }
    pack(entries, buffer, position, comm);
    for (const auto& it : data.tables()) {
        if (it.second) {
          pack(it.first, buffer, position, comm);
          pack(*it.second, buffer, position, comm);
        }
    }
}

void pack(const Equil& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.records(), buffer, position, comm);
}

void pack(const FoamConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.records(), buffer, position, comm);
    pack(data.getTransportPhase(), buffer, position, comm);
    pack(data.getMobilityModel(), buffer, position, comm);
}

void pack(const InitConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getEquil(), buffer, position, comm);
    pack(data.getFoamConfig(), buffer, position, comm);
    pack(data.filleps(), buffer, position, comm);
    pack(data.hasGravity(), buffer, position, comm);
    pack(data.restartRequested(), buffer, position, comm);
    pack(data.getRestartStep(), buffer, position, comm);
    pack(data.getRestartRootName(), buffer, position, comm);
}

void pack(const SimulationConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getThresholdPressure(), buffer, position, comm);
    pack(data.bcconfig(), buffer, position, comm);
    pack(data.rock_config(), buffer, position, comm);
    pack(data.useCPR(), buffer, position, comm);
    pack(data.hasDISGAS(), buffer, position, comm);
    pack(data.hasVAPOIL(), buffer, position, comm);
    pack(data.isThermal(), buffer, position, comm);
}

void pack(const TimeMap& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timeList(), buffer, position, comm);
}

void pack(const RestartConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timeMap(), buffer, position, comm);
    pack(data.getFirstRestartStep(), buffer, position, comm);
    pack(data.writeInitialRst(), buffer, position, comm);
    pack(data.restartSchedule(), buffer, position, comm);
    pack(data.restartKeywords(), buffer, position, comm);
    pack(data.saveKeywords(), buffer, position, comm);
}

void pack(const IOConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getWriteINITFile(), buffer, position, comm);
    pack(data.getWriteEGRIDFile(), buffer, position, comm);
    pack(data.getUNIFIN(), buffer, position, comm);
    pack(data.getUNIFOUT(), buffer, position, comm);
    pack(data.getFMTIN(), buffer, position, comm);
    pack(data.getFMTOUT(), buffer, position, comm);
    pack(data.getDeckFileName(), buffer, position, comm);
    pack(data.getOutputEnabled(), buffer, position, comm);
    pack(data.getOutputDir(), buffer, position, comm);
    pack(data.getNoSim(), buffer, position, comm);
    pack(data.getBaseName(), buffer, position, comm);
    pack(data.getEclCompatibleRST(), buffer, position, comm);
}

void pack(const Phases& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getBits(), buffer, position, comm);
}

void pack(const EndpointScaling& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getBits(), buffer, position, comm);
}

void pack(const UDQParams& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.reseed(), buffer, position, comm);
    pack(data.rand_seed(), buffer, position, comm);
    pack(data.range(), buffer, position, comm);
    pack(data.undefinedValue(), buffer, position, comm);
    pack(data.cmpEpsilon(), buffer, position, comm);
}

void pack(const Runspec& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.phases(), buffer, position, comm);
    pack(data.tabdims(), buffer, position, comm);
    pack(data.endpointScaling(), buffer, position, comm);
    pack(data.wellDimensions(), buffer, position, comm);
    pack(data.wellSegmentDimensions(), buffer, position, comm);
    pack(data.udqParams(), buffer, position, comm);
    pack(data.hysterPar(), buffer, position, comm);
    pack(data.actdims(), buffer, position, comm);
    pack(data.saturationFunctionControls(), buffer, position, comm);
}

void pack(const PvtxTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getOuterColumnSchema(), buffer, position, comm);
    pack(data.getOuterColumn(), buffer, position, comm);
    pack(data.getUnderSaturatedSchema(), buffer, position, comm);
    pack(data.getSaturatedSchema(), buffer, position, comm);
    pack(data.getUnderSaturatedTables(), buffer, position, comm);
    pack(data.getSaturatedTable(), buffer, position, comm);
}

void pack(const PvtgTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PvtxTable&>(data), buffer, position, comm);
}

void pack(const PvtoTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PvtxTable&>(data), buffer, position, comm);
}

void pack(const PvtwTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<PVTWRecord>&>(data), buffer, position, comm);
}

void pack(const PvcdoTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<PVCDORecord>&>(data), buffer, position, comm);
}

void pack(const DensityTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<DENSITYRecord>&>(data), buffer, position, comm);
}

void pack(const ViscrefTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<VISCREFRecord>&>(data), buffer, position, comm);
}

void pack(const WatdentTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<WATDENTRecord>&>(data), buffer, position, comm);
}

void pack(const PolyInjTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getThroughputs(), buffer, position, comm);
    pack(data.getVelocities(), buffer, position, comm);
    pack(data.getTableNumber(), buffer, position, comm);
    pack(data.getTableData(), buffer, position, comm);
}

void pack(const PlymwinjTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PolyInjTable&>(data), buffer, position, comm);
}

void pack(const SkprpolyTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PolyInjTable&>(data), buffer, position, comm);
    pack(data.referenceConcentration(), buffer, position, comm);
}

void pack(const SkprwatTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PolyInjTable&>(data), buffer, position, comm);
}

void pack(const RockTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<ROCKRecord>&>(data), buffer, position, comm);
}

void pack(const TableManager& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    auto simpleTables = data.getSimpleTables();
    auto splitTab = splitSimpleTable(simpleTables);

    pack(simpleTables, buffer, position, comm);
    pack(splitTab.plyshMax, buffer, position, comm);
    pack(splitTab.plyshMap, buffer, position, comm);
    pack(splitTab.rockMax, buffer, position, comm);
    pack(splitTab.rockMap, buffer, position, comm);
    pack(data.getPvtgTables(), buffer, position, comm);
    pack(data.getPvtoTables(), buffer, position, comm);
    pack(data.getRock2dTables(), buffer, position, comm);
    pack(data.getRock2dtrTables(), buffer, position, comm);
    pack(data.getPvtwTable(), buffer, position, comm);
    pack(data.getPvcdoTable(), buffer, position, comm);
    pack(data.getDensityTable(), buffer, position, comm);
    pack(data.getPlyvmhTable(), buffer, position, comm);
    pack(data.getRockTable(), buffer, position, comm);
    pack(data.getPlmixparTable(), buffer, position, comm);
    pack(data.getShrateTable(), buffer, position, comm);
    pack(data.getStone1exTable(), buffer, position, comm);
    pack(data.getTlmixparTable(), buffer, position, comm);
    pack(data.getViscrefTable(), buffer, position, comm);
    pack(data.getWatdentTable(), buffer, position, comm);
    pack(data.getPvtwSaltTables(), buffer, position, comm);
    pack(data.getBrineDensityTables(), buffer, position, comm);
    pack(data.getSolventDensityTables(), buffer, position, comm);
    pack(data.getPlymwinjTables(), buffer, position, comm);
    pack(data.getSkprwatTables(), buffer, position, comm);
    pack(data.getSkprpolyTables(), buffer, position, comm);
    pack(data.getTabdims(), buffer, position, comm);
    pack(data.getRegdims(), buffer, position, comm);
    pack(data.getEqldims(), buffer, position, comm);
    pack(data.getAqudims(), buffer, position, comm);
    pack(data.useImptvd(), buffer, position, comm);
    pack(data.useEnptvd(), buffer, position, comm);
    pack(data.useEqlnum(), buffer, position, comm);
    pack(data.useShrate(), buffer, position, comm);
    pack(data.useJFunc(), buffer, position, comm);
    if (data.useJFunc())
        pack(data.getJFunc(), buffer, position, comm);
    pack(data.OilDenT(), buffer, position, comm);
    pack(data.GasDenT(), buffer, position, comm);
    pack(data.WatDenT(), buffer, position, comm);
    pack(data.stCond(), buffer, position, comm);
    pack(data.gas_comp_index(), buffer, position, comm);
    pack(data.rtemp(), buffer, position, comm);
}

void pack(const OilVaporizationProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getType(), buffer, position, comm);
    pack(data.vap1(), buffer, position, comm);
    pack(data.vap2(), buffer, position, comm);
    pack(data.maxDRSDT(), buffer, position, comm);
    pack(data.maxDRSDT_allCells(), buffer, position, comm);
    pack(data.maxDRVDT(), buffer, position, comm);

}

void pack(const Events& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.events(), buffer, position, comm);
}

void pack(const MessageLimits& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getLimits(), buffer, position, comm);
}
void pack(const VFPInjTable& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getTableNum(), buffer, position, comm);
    pack(data.getDatumDepth(), buffer, position, comm);
    pack(data.getFloType(), buffer, position, comm);
    pack(data.getFloAxis(), buffer, position, comm);
    pack(data.getTHPAxis(), buffer, position, comm);
    pack(data.getTable(), buffer, position, comm);
}

void pack(const VFPProdTable& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getTableNum(), buffer, position, comm);
    pack(data.getDatumDepth(), buffer, position, comm);
    pack(data.getFloType(), buffer, position, comm);
    pack(data.getWFRType(), buffer, position, comm);
    pack(data.getGFRType(), buffer, position, comm);
    pack(data.getALQType(), buffer, position, comm);
    pack(data.getFloAxis(), buffer, position, comm);
    pack(data.getTHPAxis(), buffer, position, comm);
    pack(data.getWFRAxis(), buffer, position, comm);
    pack(data.getGFRAxis(), buffer, position, comm);
    pack(data.getALQAxis(), buffer, position, comm);
    pack(data.getTable(), buffer, position, comm);
}

void pack(const WellTestConfig::WTESTWell& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name, buffer, position, comm);
    pack(data.shut_reason, buffer, position, comm);
    pack(data.test_interval, buffer, position, comm);
    pack(data.num_test, buffer, position, comm);
    pack(data.startup_time, buffer, position, comm);
    pack(data.begin_report_step, buffer, position, comm);
}

void pack(const WellTestConfig& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getWells(), buffer, position, comm);
}

void pack(const WellTracerProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getConcentrations(), buffer, position, comm);
}

void pack(const UDAValue& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.get_dim(), buffer, position, comm);
    pack(data.is<double>(), buffer, position, comm);
    if (data.is<double>())
        pack(data.get<double>(), buffer, position, comm);
    else
        pack(data.get<std::string>(), buffer, position, comm);
}

void pack(const Connection& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.dir(), buffer, position, comm);
    pack(data.depth(), buffer, position, comm);
    pack(data.state(), buffer, position, comm);
    pack(data.satTableId(), buffer, position, comm);
    pack(data.complnum(), buffer, position, comm);
    pack(data.CF(), buffer, position, comm);
    pack(data.Kh(), buffer, position, comm);
    pack(data.rw(), buffer, position, comm);
    pack(data.r0(), buffer, position, comm);
    pack(data.skinFactor(), buffer, position, comm);
    pack(data.getI(), buffer, position, comm);
    pack(data.getJ(), buffer, position, comm);
    pack(data.getK(), buffer, position, comm);
    pack(data.kind(), buffer, position, comm);
    pack(data.getSeqIndex(), buffer, position, comm);
    pack(data.getSegDistStart(), buffer, position, comm);
    pack(data.getSegDistEnd(), buffer, position, comm);
    pack(data.getDefaultSatTabId(), buffer, position, comm);
    pack(data.getCompSegSeqIndex(), buffer, position, comm);
    pack(data.segment(), buffer, position, comm);
    pack(data.wellPi(), buffer, position, comm);
}

void pack(const Well::WellInjectionProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name, buffer, position, comm);
    pack(data.surfaceInjectionRate, buffer, position, comm);
    pack(data.reservoirInjectionRate, buffer, position, comm);
    pack(data.BHPTarget, buffer, position, comm);
    pack(data.THPTarget, buffer, position, comm);
    pack(data.bhp_hist_limit, buffer, position, comm);
    pack(data.thp_hist_limit, buffer, position, comm);
    pack(data.temperature, buffer, position, comm);
    pack(data.BHPH, buffer, position, comm);
    pack(data.THPH, buffer, position, comm);
    pack(data.VFPTableNumber, buffer, position, comm);
    pack(data.predictionMode, buffer, position, comm);
    pack(data.injectionControls, buffer, position, comm);
    pack(data.injectorType, buffer, position, comm);
    pack(data.controlMode, buffer, position, comm);
}

void pack(const WellEconProductionLimits& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.minOilRate(), buffer, position, comm);
    pack(data.minGasRate(), buffer, position, comm);
    pack(data.maxWaterCut(), buffer, position, comm);
    pack(data.maxGasOilRatio(), buffer, position, comm);
    pack(data.maxWaterGasRatio(), buffer, position, comm);
    pack(data.workover(), buffer, position, comm);
    pack(data.endRun(), buffer, position, comm);
    pack(data.followonWell(), buffer, position, comm);
    pack(data.quantityLimit(), buffer, position, comm);
    pack(data.maxSecondaryMaxWaterCut(), buffer, position, comm);
    pack(data.workoverSecondary(), buffer, position, comm);
    pack(data.maxGasLiquidRatio(), buffer, position, comm);
    pack(data.minLiquidRate(), buffer, position, comm);
    pack(data.maxTemperature(), buffer, position, comm);
    pack(data.minReservoirFluidRate(), buffer, position, comm);
}

void pack(const WellConnections& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getHeadI(), buffer, position, comm);
    pack(data.getHeadJ(), buffer, position, comm);
    pack(data.getNumRemoved(), buffer, position, comm);
    pack(data.getConnections(), buffer, position, comm);
}

void pack(const Well::WellProductionProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name, buffer, position, comm);
    pack(data.OilRate, buffer, position, comm);
    pack(data.WaterRate, buffer, position, comm);
    pack(data.GasRate, buffer, position, comm);
    pack(data.LiquidRate, buffer, position, comm);
    pack(data.ResVRate, buffer, position, comm);
    pack(data.BHPTarget, buffer, position, comm);
    pack(data.THPTarget, buffer, position, comm);
    pack(data.bhp_hist_limit, buffer, position, comm);
    pack(data.thp_hist_limit, buffer, position, comm);
    pack(data.BHPH, buffer, position, comm);
    pack(data.THPH, buffer, position, comm);
    pack(data.VFPTableNumber, buffer, position, comm);
    pack(data.ALQValue, buffer, position, comm);
    pack(data.predictionMode, buffer, position, comm);
    pack(data.controlMode, buffer, position, comm);
    pack(data.whistctl_cmode, buffer, position, comm);
    pack(data.getNumProductionControls(), buffer, position, comm);
}

void pack(const SpiralICD& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.strength(), buffer, position, comm);
    pack(data.length(), buffer, position, comm);
    pack(data.densityCalibration(), buffer, position, comm);
    pack(data.viscosityCalibration(), buffer, position, comm);
    pack(data.criticalValue(), buffer, position, comm);
    pack(data.widthTransitionRegion(), buffer, position, comm);
    pack(data.maxViscosityRatio(), buffer, position, comm);
    pack(data.methodFlowScaling(), buffer, position, comm);
    pack(data.maxAbsoluteRate(), buffer, position, comm);
    pack(data.status(), buffer, position, comm);
    pack(data.scalingFactor(), buffer, position, comm);
}

void pack(const Valve& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.conEFlowCoefficient(), buffer, position, comm);
    pack(data.conCrossArea(), buffer, position, comm);
    pack(data.conMaxCrossArea(), buffer, position, comm);
    pack(data.pipeAdditionalLength(), buffer, position, comm);
    pack(data.pipeDiameter(), buffer, position, comm);
    pack(data.pipeRoughness(), buffer, position, comm);
    pack(data.pipeCrossArea(), buffer, position, comm);
    pack(data.status(), buffer, position, comm);
}

void pack(const Segment& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.segmentNumber(), buffer, position, comm);
    pack(data.branchNumber(), buffer, position, comm);
    pack(data.outletSegment(), buffer, position, comm);
    pack(data.inletSegments(), buffer, position, comm);
    pack(data.totalLength(), buffer, position, comm);
    pack(data.depth(), buffer, position, comm);
    pack(data.internalDiameter(), buffer, position, comm);
    pack(data.roughness(), buffer, position, comm);
    pack(data.crossArea(), buffer, position, comm);
    pack(data.volume(), buffer, position, comm);
    pack(data.dataReady(), buffer, position, comm);
    pack(data.segmentType(), buffer, position, comm);
    pack(data.spiralICD(), buffer, position, comm);
    pack(data.getValve(), buffer, position, comm);
}

template<class T>
void pack(const std::shared_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data != nullptr, buffer, position, comm);
    if (data)
        pack(*data, buffer, position, comm);
}

template<class T>
void pack(const std::unique_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data != nullptr, buffer, position, comm);
    if (data)
        pack(*data, buffer, position, comm);
}

void pack(const Dimension& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getSIScalingRaw(), buffer, position, comm);
    pack(data.getSIOffset(), buffer, position, comm);
}

void pack(const UnitSystem& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getName(), buffer, position, comm);
    pack(data.getType(), buffer, position, comm);
    pack(data.getDimensions(), buffer, position, comm);
    pack(data.use_count(), buffer, position, comm);
}

void pack(const WellSegments& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.compPressureDrop(), buffer, position, comm);
    pack(data.segments(), buffer, position, comm);
}

void pack(const Well& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.groupName(), buffer, position, comm);
    pack(data.firstTimeStep(), buffer, position, comm);
    pack(data.seqIndex(), buffer, position, comm);
    pack(data.getHeadI(), buffer, position, comm);
    pack(data.getHeadJ(), buffer, position, comm);
    pack(data.getRefDepth(), buffer, position, comm);
    pack(data.wellType(), buffer, position, comm);
    pack(data.getWellConnectionOrdering(), buffer, position, comm);
    pack(data.units(), buffer, position, comm);
    pack(data.udqUndefined(), buffer, position, comm);
    pack(data.getStatus(), buffer, position, comm);
    pack(data.getDrainageRadius(), buffer, position, comm);
    pack(data.getAllowCrossFlow(), buffer, position, comm);
    pack(data.getAutomaticShutIn(), buffer, position, comm);
    pack(data.wellGuideRate(), buffer, position, comm);
    pack(data.getEfficiencyFactor(), buffer, position, comm);
    pack(data.getSolventFraction(), buffer, position, comm);
    pack(data.predictionMode(), buffer, position, comm);
    pack(data.getEconLimits(), buffer, position, comm);
    pack(data.getFoamProperties(), buffer, position, comm);
    pack(data.getPolymerProperties(), buffer, position, comm);
    pack(data.getBrineProperties(), buffer, position, comm);
    pack(data.getTracerProperties(), buffer, position, comm);
    pack(data.getConnections(), buffer, position, comm);
    pack(data.getProductionProperties(), buffer, position, comm);
    pack(data.getInjectionProperties(), buffer, position, comm);
    pack(data.hasSegments(), buffer, position, comm);
    if (data.hasSegments())
        pack(data.getSegments(), buffer, position, comm);
}

template<class T>
void pack(const IOrderSet<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.index(), buffer, position, comm);
    pack(data.data(), buffer, position, comm);
}

void pack(const Group::GroupInjectionProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.phase, buffer, position, comm);
    pack(data.cmode, buffer, position, comm);
    pack(data.surface_max_rate, buffer, position, comm);
    pack(data.resv_max_rate, buffer, position, comm);
    pack(data.target_reinj_fraction, buffer, position, comm);
    pack(data.target_void_fraction, buffer, position, comm);
    pack(data.reinj_group, buffer, position, comm);
    pack(data.voidage_group, buffer, position, comm);
    pack(data.injection_controls, buffer, position, comm);
}

void pack(const Group::GroupProductionProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.cmode, buffer, position, comm);
    pack(data.exceed_action, buffer, position, comm);
    pack(data.oil_target, buffer, position, comm);
    pack(data.water_target, buffer, position, comm);
    pack(data.gas_target, buffer, position, comm);
    pack(data.liquid_target, buffer, position, comm);
    pack(data.guide_rate, buffer, position, comm);
    pack(data.guide_rate_def, buffer, position, comm);
    pack(data.resv_target, buffer, position, comm);
    pack(data.production_controls, buffer, position, comm);
}

void pack(const Group& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.insert_index(), buffer, position, comm);
    pack(data.initStep(), buffer, position, comm);
    pack(data.udqUndefined(), buffer, position, comm);
    pack(data.units(), buffer, position, comm);
    pack(data.type(), buffer, position, comm);
    pack(data.getGroupEfficiencyFactor(), buffer, position, comm);
    pack(data.getTransferGroupEfficiencyFactor(), buffer, position, comm);
    pack(data.isAvailableForGroupControl(), buffer, position, comm);
    pack(data.getGroupNetVFPTable(), buffer, position, comm);
    pack(data.parent(), buffer, position, comm);
    pack(data.iwells(), buffer, position, comm);
    pack(data.igroups(), buffer, position, comm);
    pack(data.injectionProperties(), buffer, position, comm);
    pack(data.productionProperties(), buffer, position, comm);
}

void pack(const WList& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.wellList(), buffer, position, comm);
}

void pack(const WListManager& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.lists(), buffer, position, comm);
}

void pack(const UDQASTNode& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.var_type, buffer, position, comm);
    pack(data.getType(), buffer, position, comm);
    pack(data.stringValue(), buffer, position, comm);
    pack(data.scalarValue(), buffer, position, comm);
    pack(data.getSelectors(), buffer, position, comm);
    pack(data.getLeft(), buffer, position, comm);
    pack(data.getRight(), buffer, position, comm);
}

void pack(const UDQDefine& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.keyword(), buffer, position, comm);
    pack(data.getAst(), buffer, position, comm);
    pack(data.var_type(), buffer, position, comm);
    pack(data.input_string(), buffer, position, comm);
}

void pack(const UDQAssign::AssignRecord& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.selector, buffer, position, comm);
    pack(data.value, buffer, position, comm);
}

void pack(const UDQAssign& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.keyword(), buffer, position, comm);
    pack(data.var_type(), buffer, position, comm);
    pack(data.getRecords(), buffer, position, comm);
}

void pack(const UDQIndex& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.insert_index, buffer, position, comm);
    pack(data.typed_insert_index, buffer, position, comm);
    pack(data.action, buffer, position, comm);
    pack(data.var_type, buffer, position, comm);
}

void pack(const UDQConfig& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.params(), buffer, position, comm);
    pack(data.definitionMap(), buffer, position, comm);
    pack(data.assignmentMap(), buffer, position, comm);
    pack(data.unitsMap(), buffer, position, comm);
    pack(data.inputIndex(), buffer, position, comm);
    pack(data.typeCount(), buffer, position, comm);
}

void pack(const UDQActive::InputRecord& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.input_index, buffer, position, comm);
    pack(data.udq, buffer, position, comm);
    pack(data.wgname, buffer, position, comm);
    pack(data.control, buffer, position, comm);
}

void pack(const UDQActive::Record& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.udq, buffer, position, comm);
    pack(data.input_index, buffer, position, comm);
    pack(data.use_index, buffer, position, comm);
    pack(data.wgname, buffer, position, comm);
    pack(data.control, buffer, position, comm);
    pack(data.uad_code, buffer, position, comm);
    pack(data.use_count, buffer, position, comm);
}

void pack(const UDQActive& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getInputRecords(), buffer, position, comm);
    pack(data.getOutputRecords(), buffer, position, comm);
    pack(data.getUdqKeys(), buffer, position, comm);
    pack(data.getWgKeys(), buffer, position, comm);
}

void pack(const GuideRateModel& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timeInterval(), buffer, position, comm);
    pack(data.target(), buffer, position, comm);
    pack(data.coefs(), buffer, position, comm);
    pack(data.allow_increase(), buffer, position, comm);
    pack(data.damping_factor(), buffer, position, comm);
    pack(data.free_gas(), buffer, position, comm);
    pack(data.defaultModel(), buffer, position, comm);
    pack(data.udaCoefs(), buffer, position, comm);
}

void pack(const GuideRateConfig& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getModel(), buffer, position, comm);
    pack(data.getWells(), buffer, position, comm);
    pack(data.getGroups(), buffer, position, comm);
}

void pack(const GConSale::GCONSALEGroup& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.sales_target, buffer, position, comm);
    pack(data.max_sales_rate, buffer, position, comm);
    pack(data.min_sales_rate, buffer, position, comm);
    pack(data.max_proc, buffer, position, comm);
    pack(data.udq_undefined, buffer, position, comm);
    pack(data.unit_system, buffer, position, comm);
}

void pack(const GConSale& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getGroups(), buffer, position, comm);
}

void pack(const GConSump::GCONSUMPGroup& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.consumption_rate, buffer, position, comm);
    pack(data.import_rate, buffer, position, comm);
    pack(data.network_node, buffer, position, comm);
    pack(data.udq_undefined, buffer, position, comm);
    pack(data.unit_system, buffer, position, comm);
}

void pack(const GConSump& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getGroups(), buffer, position, comm);
}

void pack(const RFTConfig& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timeMap(), buffer, position, comm);
    pack(data.firstRFTOutput(), buffer, position, comm);
    pack(data.wellOpenRftTime(), buffer, position, comm);
    pack(data.wellOpenRftName(), buffer, position, comm);
    pack(data.wellOpen(), buffer, position, comm);
    pack(data.rftConfig(), buffer, position, comm);
    pack(data.pltConfig(), buffer, position, comm);
}

void pack(const DeckItem& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.dVal(), buffer, position, comm);
    pack(data.iVal(), buffer, position, comm);
    pack(data.sVal(), buffer, position, comm);
    pack(data.uVal(), buffer, position, comm);
    pack(data.getType(), buffer, position, comm);
    pack(data.name(), buffer, position, comm);
    pack(data.valueStatus(), buffer, position, comm);
    pack(data.rawData(), buffer, position, comm);
    pack(data.activeDimensions(), buffer, position, comm);
    pack(data.defaultDimensions(), buffer, position, comm);
}

void pack(const DeckRecord& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getItems(), buffer, position, comm);
}

void pack(const Location& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.filename, buffer, position, comm);
    pack(data.lineno, buffer, position, comm);
}

void pack(const DeckKeyword& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.location(), buffer, position, comm);
    pack(data.records(), buffer, position, comm);
    pack(data.isDataKeyword(), buffer, position, comm);
    pack(data.isSlashTerminated(), buffer, position, comm);
}

void pack(const Deck& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.keywords(), buffer, position, comm);
    pack(data.getDefaultUnitSystem(), buffer, position, comm);
    pack(data.activeUnitSystem(), buffer, position, comm);
    pack(data.getDataFile(), buffer, position, comm);
    pack(data.getInputPath(), buffer, position, comm);
    pack(data.unitSystemAccessCount(), buffer, position, comm);
}

void pack(const Action::ASTNode& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.type, buffer, position, comm);
    pack(data.func_type, buffer, position, comm);
    pack(data.func, buffer, position, comm);
    pack(data.argList(), buffer, position, comm);
    pack(data.getNumber(), buffer, position, comm);
    pack(data.childrens(), buffer, position, comm);
}

void pack(const Action::AST& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getCondition(), buffer, position, comm);
}

void pack(const Action::Quantity& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.quantity, buffer, position, comm);
    pack(data.args, buffer, position, comm);
}

void pack(const Action::Condition& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.lhs, buffer, position, comm);
    pack(data.rhs, buffer, position, comm);
    pack(data.logic, buffer, position, comm);
    pack(data.cmp, buffer, position, comm);
    pack(data.cmp_string, buffer, position, comm);
}

void pack(const Action::ActionX& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.max_run(), buffer, position, comm);
    pack(data.min_wait(), buffer, position, comm);
    pack(data.start_time(), buffer, position, comm);
    pack(data.getKeywords(), buffer, position, comm);
    pack(data.getCondition(), buffer, position, comm);
    pack(data.conditions(), buffer, position, comm);
    pack(data.getRunCount(), buffer, position, comm);
    pack(data.getLastRun(), buffer, position, comm);
}

void pack(const Action::Actions& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getActions(), buffer, position, comm);
}

void pack(const Schedule& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getTimeMap(), buffer, position, comm);
    packDynMap(data.getStaticWells(), buffer, position, comm);
    packDynMap(data.getGroups(), buffer, position, comm);
    pack(data.getOilVapProps(), buffer, position, comm);
    pack(data.getEvents(), buffer, position, comm);
    pack(data.getModifierDeck(), buffer, position, comm);
    pack(data.getTuning(), buffer, position, comm);
    pack(data.getMessageLimits(), buffer, position, comm);
    pack(data.getRunspec(), buffer, position, comm);
    packDynMap<Map2>(data.getVFPProdTables(), buffer, position, comm);
    packDynMap<Map2>(data.getVFPInjTables(), buffer, position, comm);
    pack(data.getWellTestConfig(), buffer, position, comm);
    pack(data.getWListManager(), buffer, position, comm);
    pack(data.getUDQConfig(), buffer, position, comm);
    pack(data.getUDQActive(), buffer, position, comm);
    pack(data.getGuideRateConfig(), buffer, position, comm);
    pack(data.getGConSale(), buffer, position, comm);
    pack(data.getGConSump(), buffer, position, comm);
    pack(data.getGlobalWhistCtlMode(), buffer, position, comm);
    pack(data.getActions(), buffer, position, comm);
    pack(data.rftConfig(), buffer, position, comm);
    pack(data.getNupCol(), buffer, position, comm);
    pack(data.restart(), buffer, position, comm);
    pack(data.getWellGroupEvents(), buffer, position, comm);
}

void pack(const BrineDensityTable& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getBrineDensityColumn(), buffer, position, comm);
}

void pack(const PvtwsaltTable& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getReferencePressureValue(), buffer, position, comm);
    pack(data.getReferenceSaltConcentrationValue(), buffer, position, comm);
    pack(data.getTableValues(), buffer, position, comm);
}

void pack(const EquilRecord& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.datumDepth(), buffer, position, comm);
    pack(data.datumDepthPressure(), buffer, position, comm);
    pack(data.waterOilContactDepth(), buffer, position, comm);
    pack(data.waterOilContactCapillaryPressure(), buffer, position, comm);
    pack(data.gasOilContactDepth(), buffer, position, comm);
    pack(data.gasOilContactCapillaryPressure(), buffer, position, comm);
    pack(data.liveOilInitConstantRs(), buffer, position, comm);
    pack(data.wetGasInitConstantRv(), buffer, position, comm);
    pack(data.initializationTargetAccuracy(), buffer, position, comm);
}

void pack(const FoamData& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.referenceSurfactantConcentration(), buffer, position, comm);
    pack(data.exponent(), buffer, position, comm);
    pack(data.minimumSurfactantConcentration(), buffer, position, comm);
    pack(data.allowDesorption(), buffer, position, comm);
    pack(data.rockDensity(), buffer, position, comm);
}

void pack(const RestartSchedule& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timestep, buffer, position, comm);
    pack(data.basic, buffer, position, comm);
    pack(data.frequency, buffer, position, comm);
    pack(data.rptsched_restart_set, buffer, position, comm);
    pack(data.rptsched_restart, buffer, position, comm);
}

void pack(const TimeStampUTC& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.ymd(), buffer, position, comm);
    pack(data.hour(), buffer, position, comm);
    pack(data.minutes(), buffer, position, comm);
    pack(data.seconds(), buffer, position, comm);
    pack(data.microseconds(), buffer, position, comm);
}

void pack(const EclHysterConfig& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.active(), buffer, position, comm);
    pack(data.pcHysteresisModel(), buffer, position, comm);
    pack(data.krHysteresisModel(), buffer, position, comm);
}

void pack(const JFunc& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.flag(), buffer, position, comm);
    pack(data.owSurfaceTension(), buffer, position, comm);
    pack(data.goSurfaceTension(), buffer, position, comm);
    pack(data.alphaFactor(), buffer, position, comm);
    pack(data.betaFactor(), buffer, position, comm);
    pack(data.direction(), buffer, position, comm);
}

void pack(const WellPolymerProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.m_polymerConcentration, buffer, position, comm);
    pack(data.m_saltConcentration, buffer, position, comm);
    pack(data.m_plymwinjtable, buffer, position, comm);
    pack(data.m_skprwattable, buffer, position, comm);
    pack(data.m_skprpolytable, buffer, position, comm);
}

void pack(const Well::WellGuideRate& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.available, buffer, position, comm);
    pack(data.guide_rate, buffer, position, comm);
    pack(data.guide_phase, buffer, position, comm);
    pack(data.scale_factor, buffer, position, comm);
}

void pack(const GuideRateConfig::WellTarget& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.guide_rate, buffer, position, comm);
    pack(data.target, buffer, position, comm);
    pack(data.scaling_factor, buffer, position, comm);
}

void pack(const GuideRateConfig::GroupTarget& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.guide_rate, buffer, position, comm);
    pack(data.target, buffer, position, comm);
}

void pack(const MULTREGTRecord& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.src_value, buffer, position, comm);
    pack(data.target_value, buffer, position, comm);
    pack(data.trans_mult, buffer, position, comm);
    pack(data.directions, buffer, position, comm);
    pack(data.nnc_behaviour, buffer, position, comm);
    pack(data.region_name, buffer, position, comm);
}

void pack(const MULTREGTScanner& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getSize(), buffer, position, comm);
    pack(data.getRecords(), buffer, position, comm);
    pack(data.getSearchMap(), buffer, position, comm);
    pack(data.getRegions(), buffer, position, comm);
    pack(data.getDefaultRegion(), buffer, position, comm);
}

void pack(const EclipseConfig& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.init(), buffer, position, comm);
    pack(data.io(), buffer, position, comm);
}

void pack(const TransMult& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getSize(), buffer, position, comm);
    pack(data.getTrans(), buffer, position, comm);
    pack(data.getNames(), buffer, position, comm);
    pack(data.getScanner(), buffer, position, comm);
}

void pack(const FaultFace& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getIndices(), buffer, position, comm);
    pack(data.getDir(), buffer, position, comm);
}

void pack(const Fault& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getName(), buffer, position, comm);
    pack(data.getTransMult(), buffer, position, comm);
    pack(data.getFaceList(), buffer, position, comm);
}

void pack(const FaultCollection& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getFaults(), buffer, position, comm);
}

void pack(const SolventDensityTable& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getSolventDensityColumn(), buffer, position, comm);
}

void pack(const GridDims& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getNXYZ(), buffer, position, comm);
}

void pack(const ShrateTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<ShrateRecord>&>(data), buffer, position, comm);
}

void pack(const TlmixparTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<TlmixparRecord>&>(data), buffer, position, comm);
}

void pack(const PlmixparTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<PlmixparRecord>&>(data), buffer, position, comm);
}

void pack(const PlyvmhTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<PlyvmhRecord>&>(data), buffer, position, comm);
}

void pack(const Stone1exTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<Stone1exRecord>&>(data), buffer, position, comm);
}

void pack(const PlyshlogTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const SimpleTable&>(data), buffer, position, comm);
    pack(data.getRefPolymerConcentration(), buffer, position, comm);
    pack(data.getRefSalinity(), buffer, position, comm);
    pack(data.getRefTemperature(), buffer, position, comm);
    pack(data.hasRefSalinity(), buffer, position, comm);
    pack(data.hasRefTemperature(), buffer, position, comm);
}

void pack(const RocktabTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const SimpleTable&>(data), buffer, position, comm);
    pack(data.isDirectional(), buffer, position, comm);
}

/// unpack routines

template<class T>
void unpack(T*, const std::size_t&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    EWOMS_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
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

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data, l, buffer, position, comm, typename std::is_pod<T>::type());
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
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t length = 0;
    unpack(length, buffer, position, comm);
    data.resize(length);

    if (std::is_pod<T>::value)
    {
        unpack(data.data(), data.size(), buffer, position, comm);
        return;
    }

    for (auto& entry: data)
        unpack(entry, buffer, position, comm);
}

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

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, void>::type
unpack_tuple_entry(Tuple&, std::vector<char>&, int&,
                   Dune::MPIHelper::MPICommunicator)
{
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, void>::type
unpack_tuple_entry(Tuple& tuple, std::vector<char>& buffer,
                   int& position, Dune::MPIHelper::MPICommunicator comm)
{
    unpack(std::get<I>(tuple), buffer, position, comm);
    unpack_tuple_entry<I+1>(tuple, buffer, position, comm);
}

template<class... Ts>
void unpack(std::tuple<Ts...>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm)
{
    unpack_tuple_entry(data, buffer, position, comm);
}

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

template<class T, size_t N>
void unpack(std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    for (T& entry : data)
        unpack(entry, buffer, position, comm);
}

template<class Key, class Value>
void unpack(OrderedMap<Key,Value>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
  typename OrderedMap<Key,Value>::index_type index;
  typename OrderedMap<Key,Value>::storage_type storage;
  unpack(index, buffer, position, comm);
  unpack(storage, buffer, position, comm);
  data = OrderedMap<Key,Value>(index, storage);
}

template<class T>
void unpack(DynamicState<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<T> unique;
    std::vector<int> indices;
    Ewoms::Mpi::unpack(unique, buffer, position, comm);
    Ewoms::Mpi::unpack(indices, buffer, position, comm);
    reconstructDynState(unique, indices, data);
}

template<class T>
void unpack(DynamicVector<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<T> ddata;
    unpack(ddata, buffer, position, comm);
    data = DynamicVector<T>(ddata);
}

void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
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

void unpack(std::string& str, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t length=0;
    unpack(length, buffer, position, comm);
    std::vector<char> cStr(length, '\0');
    unpack(cStr.data(), length, buffer, position, comm);
    assert(str.empty());
    str.append(cStr.data());
}

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

template<class T1, class T2, class H, class P, class A>
void unpack(std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
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

void unpack(data::Well& data, std::vector<char>& buffer, int& position,
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

void unpack(RestartKey& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.key, buffer, position, comm);
    unpack(data.dim, buffer, position, comm);
    unpack(data.required, buffer, position, comm);
}

void unpack(data::CellData& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.dim, buffer, position, comm);
    unpack(data.data, buffer, position, comm);
    unpack(data.target, buffer, position, comm);
}

void unpack(data::Solution& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    unpack(static_cast<std::map< std::string, data::CellData>&>(data),
           buffer, position, comm);
}

void unpack(data::WellRates& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    unpack(static_cast<std::map< std::string, data::Well>&>(data),
           buffer, position, comm);
}

void unpack(RestartValue& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.solution, buffer, position, comm);
    unpack(data.wells, buffer, position, comm);
    unpack(data.extra, buffer, position, comm);
}

void unpack(RockConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    RockConfig rock_config;
    bool active;
    std::vector<RockConfig::RockComp> rock_comp;
    std::string rocknum_property;
    std::size_t num_rock_tables;
    bool water_compaction;
    RockConfig::Hysteresis hyst_mode;

    unpack(active, buffer, position, comm);
    unpack(rocknum_property, buffer, position, comm);
    unpack(rock_comp, buffer, position, comm);
    unpack(num_rock_tables, buffer, position, comm);
    unpack(water_compaction, buffer, position, comm);
    unpack(hyst_mode, buffer, position, comm);
    data = RockConfig(active, rock_comp, rocknum_property, num_rock_tables, water_compaction, hyst_mode);
}

void unpack(ThresholdPressure& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    ThresholdPressure::ThresholdPressureTable thpTable;
    ThresholdPressure::PressureTable pTable;
    bool active, restart;
    unpack(active, buffer, position, comm);
    unpack(restart, buffer, position, comm);
    unpack(thpTable, buffer, position, comm);
    unpack(pTable, buffer, position, comm);

    data = ThresholdPressure(active, restart, thpTable, pTable);
}

void unpack(AquiferCT::AQUCT_data& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm)
{
    int aquiferID;
    int inftableID, pvttableID;
    double phi_aq, d0, C_t, r_o, k_a, c1, h, theta, c2;
    std::pair<bool, double> p0;
    std::vector<double> td, pi;
    std::vector<int> cell_id;

    unpack(aquiferID, buffer, position, comm);
    unpack(inftableID, buffer, position, comm);
    unpack(pvttableID, buffer, position, comm);
    unpack(phi_aq, buffer, position, comm);
    unpack(d0, buffer, position, comm);
    unpack(C_t, buffer, position, comm);
    unpack(r_o, buffer, position, comm);
    unpack(k_a, buffer, position, comm);
    unpack(c1, buffer, position, comm);
    unpack(h, buffer, position, comm);
    unpack(theta, buffer, position, comm);
    unpack(c2, buffer, position, comm);
    unpack(p0, buffer, position, comm);
    unpack(td, buffer, position, comm);
    unpack(pi, buffer, position, comm);
    unpack(cell_id, buffer, position, comm);

    data = AquiferCT::AQUCT_data(aquiferID,
                                 inftableID,
                                 pvttableID,
                                 phi_aq,
                                 d0,
                                 C_t,
                                 r_o,
                                 k_a,
                                 c1,
                                 h,
                                 theta,
                                 c2,
                                 p0,
                                 td,
                                 pi,
                                 cell_id);
}

void unpack(WellType& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm)
{
    Phase preferred_phase;
    bool producer;
    unpack(producer, buffer, position, comm);
    unpack(preferred_phase, buffer, position, comm);
    data = WellType( producer, preferred_phase );
}

void unpack(DenT& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<DenT::entry> records;
    unpack(records, buffer, position, comm);
    data = DenT( records );
}

void unpack(AquiferCT& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<AquiferCT::AQUCT_data> aquiferList;
    unpack(aquiferList, buffer, position, comm);
    data = AquiferCT(aquiferList);
}

void unpack(Aquifetp::AQUFETP_data& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm)
{
    int aquiferID;
    int pvttableID;
    double  J, C_t, V0, d0;
    std::pair<bool, double> p0;

    unpack(aquiferID, buffer, position, comm);
    unpack(pvttableID, buffer, position, comm);
    unpack(J, buffer, position, comm);
    unpack(C_t, buffer, position, comm);
    unpack(V0, buffer, position, comm);
    unpack(d0, buffer, position, comm);
    unpack(p0, buffer, position, comm);
    data = Aquifetp::AQUFETP_data(aquiferID, pvttableID, J, C_t, V0, d0, p0);
}

void unpack(Aquifetp& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Aquifetp::AQUFETP_data> aquiferList;
    unpack(aquiferList, buffer, position, comm);
    data = Aquifetp(aquiferList);
}

void unpack(Aquancon::AquancCell& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm)
{
    int aquiferID;
    std::size_t globalIndex;
    std::pair<bool, double> influxCoeff;
    double influxMult;
    FaceDir::DirEnum faceDir;

    unpack(aquiferID, buffer, position, comm);
    unpack(globalIndex, buffer, position, comm);
    unpack(influxCoeff, buffer, position, comm);
    unpack(influxMult, buffer, position, comm);
    unpack(faceDir, buffer, position, comm);

    data = Aquancon::AquancCell(aquiferID, globalIndex, influxCoeff, influxMult, faceDir);
}

void unpack(AquiferConfig& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm) {
    Aquifetp fetp;
    AquiferCT ct;
    Aquancon conn;

    unpack(fetp, buffer, position, comm);
    unpack(ct, buffer, position, comm);
    unpack(conn, buffer, position, comm);
    data = AquiferConfig(fetp, ct, conn);
}

void unpack(Aquancon& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm)
{
    std::unordered_map<int, std::vector<Aquancon::AquancCell>> aquiferCells;
    unpack(aquiferCells, buffer, position, comm);
    data = Aquancon(aquiferCells);
}

void unpack(BCConfig& bc, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<BCConfig::BCFace> faces;

    unpack(faces, buffer, position, comm);
    bc = BCConfig(faces);
}

void unpack(NNC& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<NNCdata> res;
    unpack(res, buffer, position, comm);
    data = NNC(res);
}

void unpack(EDITNNC& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<NNCdata> res;
    unpack(res, buffer, position, comm);
    data = EDITNNC(res);
}

void unpack(Rock2dTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<std::vector<double>> pvmultValues;
    std::vector<double> pressureValues;
    unpack(pvmultValues, buffer, position, comm);
    unpack(pressureValues, buffer, position, comm);
    data = Rock2dTable(pvmultValues, pressureValues);
}

void unpack(Rock2dtrTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<std::vector<double>> transMultValues;
    std::vector<double> pressureValues;
    unpack(transMultValues, buffer, position, comm);
    unpack(pressureValues, buffer, position, comm);
    data = Rock2dtrTable(transMultValues, pressureValues);
}

void unpack(ColumnSchema& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    Table::ColumnOrderEnum order;
    Table::DefaultAction action;
    unpack(name, buffer, position, comm);
    unpack(order, buffer, position, comm);
    unpack(action, buffer, position, comm);
    if (action == Table::DEFAULT_CONST) {
        double value;
        unpack(value, buffer, position, comm);
        data = ColumnSchema(name, order, value);
    } else
        data = ColumnSchema(name, order, action);
}

void unpack(TableSchema& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    OrderedMap<std::string, ColumnSchema> columns;
    unpack(columns, buffer, position, comm);
    data = TableSchema(columns);
}

void unpack(TableColumn& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    ColumnSchema schema;
    std::string name;
    std::vector<double> values;
    std::vector<bool> defaults;
    size_t defaultCount;
    unpack(schema, buffer, position, comm);
    unpack(name, buffer, position, comm);
    unpack(values, buffer, position, comm);
    unpack(defaults, buffer, position, comm);
    unpack(defaultCount, buffer, position, comm);
    data = TableColumn(schema, name, values, defaults, defaultCount);
}

void unpack(SimpleTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TableSchema schema;
    OrderedMap<std::string, TableColumn> columns;
    bool jf;
    unpack(schema, buffer, position, comm);
    unpack(columns, buffer, position, comm);
    unpack(jf, buffer, position, comm);
    data = SimpleTable(schema, columns, jf);
}

void unpack(TableContainer& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    size_t max;
    unpack(max, buffer, position, comm);
    data = TableContainer(max);
    size_t entries;
    unpack(entries, buffer, position, comm);
    for (size_t i = 0; i < entries; ++i) {
        size_t id;
        unpack(id, buffer, position, comm);
        SimpleTable table;
        unpack(table, buffer, position, comm);
        data.addTable(id, std::make_shared<SimpleTable>(table));
    }
}

void unpack(Equil& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<EquilRecord> records;
    unpack(records, buffer, position, comm);
    data = Equil(records);
}

void unpack(FoamConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<FoamData> records;
    Phase transport_phase;
    FoamConfig::MobilityModel mobility_model;
    unpack(records, buffer, position, comm);
    unpack(transport_phase, buffer, position, comm);
    unpack(mobility_model, buffer, position, comm);
    data = FoamConfig(records, transport_phase, mobility_model);
}

void unpack(InitConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    Equil equil;
    FoamConfig foam;
    bool filleps, hasGravity, restartRequested;
    int restartStep;
    std::string restartRootName;
    unpack(equil, buffer, position, comm);
    unpack(foam, buffer, position, comm);
    unpack(filleps, buffer, position, comm);
    unpack(hasGravity, buffer, position, comm);
    unpack(restartRequested, buffer, position, comm);
    unpack(restartStep, buffer, position, comm);
    unpack(restartRootName, buffer, position, comm);
    data = InitConfig(equil, foam, filleps, hasGravity,
                      restartRequested, restartStep, restartRootName);
}

void unpack(SimulationConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    ThresholdPressure thresholdPressure;
    BCConfig bc;
    RockConfig rock_config;
    bool useCPR, DISGAS, VAPOIL, isThermal;
    unpack(thresholdPressure, buffer, position, comm);
    unpack(bc, buffer, position, comm);
    unpack(rock_config, buffer, position, comm);
    unpack(useCPR, buffer, position, comm);
    unpack(DISGAS, buffer, position, comm);
    unpack(VAPOIL, buffer, position, comm);
    unpack(isThermal, buffer, position, comm);
    data = SimulationConfig(thresholdPressure, bc, rock_config, useCPR, DISGAS, VAPOIL, isThermal);
}

void unpack(TimeMap& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<std::time_t> timeList;
    unpack(timeList, buffer, position, comm);

    data = TimeMap(timeList);
}

void unpack(RestartConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TimeMap timemap;
    int firstRstStep;
    bool writeInitialRst;
    DynamicState<RestartSchedule> restart_sched;
    DynamicState<std::map<std::string,int>> restart_keyw;
    std::vector<bool> save_keyw;
    unpack(timemap, buffer, position, comm);
    unpack(firstRstStep, buffer, position, comm);
    unpack(writeInitialRst, buffer, position, comm);
    unpack(restart_sched, buffer, position, comm);
    unpack(restart_keyw, buffer, position, comm);
    unpack(save_keyw, buffer, position, comm);
    data = RestartConfig(timemap, firstRstStep, writeInitialRst, restart_sched,
                         restart_keyw, save_keyw);
}

void unpack(IOConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    bool write_init, write_egrid, unifin, unifout, fmtin, fmtout;
    std::string deck_name, output_dir, base_name;
    bool output_enabled, no_sim, ecl_compatible_rst;

    unpack(write_init, buffer, position, comm);
    unpack(write_egrid, buffer, position, comm);
    unpack(unifin, buffer, position, comm);
    unpack(unifout, buffer, position, comm);
    unpack(fmtin, buffer, position, comm);
    unpack(fmtout, buffer, position, comm);
    unpack(deck_name, buffer, position, comm);
    unpack(output_enabled, buffer, position, comm);
    unpack(output_dir, buffer, position, comm);
    unpack(no_sim, buffer, position, comm);
    unpack(base_name, buffer, position, comm);
    unpack(ecl_compatible_rst, buffer, position, comm);
    data = IOConfig(write_init, write_egrid, unifin, unifout, fmtin, fmtout,
                    deck_name, output_enabled, output_dir,
                    no_sim, base_name, ecl_compatible_rst);
}

void unpack(Phases& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unsigned long bits;
    unpack(bits, buffer, position, comm);
    data = Phases(std::bitset<NUM_PHASES_IN_ENUM>(bits));
}

void unpack(EndpointScaling& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unsigned long bits;
    unpack(bits, buffer, position, comm);
    data = EndpointScaling(std::bitset<4>(bits));
}

void unpack(UDQParams& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    bool reseed;
    int rand_seed;
    double range, undefVal, cmp_eps;

    unpack(reseed, buffer, position, comm);
    unpack(rand_seed, buffer, position, comm);
    unpack(range, buffer, position, comm);
    unpack(undefVal, buffer, position, comm);
    unpack(cmp_eps, buffer, position, comm);
    data = UDQParams(reseed, rand_seed, range, undefVal, cmp_eps);
}

void unpack(Runspec& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    Phases phases;
    Tabdims tabdims;
    EndpointScaling endScale;
    Welldims wellDims;
    WellSegmentDims wsegDims;
    UDQParams udqparams;
    EclHysterConfig hystPar;
    Actdims actdims;
    SatFuncControls sfuncctrl;

    unpack(phases, buffer, position, comm);
    unpack(tabdims, buffer, position, comm);
    unpack(endScale, buffer, position, comm);
    unpack(wellDims, buffer, position, comm);
    unpack(wsegDims, buffer, position, comm);
    unpack(udqparams, buffer, position, comm);
    unpack(hystPar, buffer, position, comm);
    unpack(actdims, buffer, position, comm);
    unpack(sfuncctrl, buffer, position, comm);
    data = Runspec(phases, tabdims, endScale, wellDims, wsegDims,
                   udqparams, hystPar, actdims, sfuncctrl);
}

template<class PVTType>
void unpack_pvt(PVTType& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    ColumnSchema outer_schema;
    TableColumn outer_column;
    TableSchema undersat_schema, sat_schema;
    std::vector<SimpleTable> undersat_tables;
    SimpleTable sat_table;
    unpack(outer_schema, buffer, position, comm);
    unpack(outer_column, buffer, position, comm);
    unpack(undersat_schema, buffer, position, comm);
    unpack(sat_schema, buffer, position, comm);
    unpack(undersat_tables, buffer, position, comm);
    unpack(sat_table, buffer, position, comm);
    data = PVTType(outer_schema, outer_column, undersat_schema, sat_schema,
                   undersat_tables, sat_table);
}

void unpack(PvtgTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack_pvt(data, buffer, position, comm);
}

void unpack(PvtoTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack_pvt(data, buffer, position, comm);
}

void unpack(PvtwTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<PVTWRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = PvtwTable(pdata);
}

void unpack(PvcdoTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<PVCDORecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = PvcdoTable(pdata);
}

void unpack(DensityTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<DENSITYRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = DensityTable(pdata);
}

void unpack(ViscrefTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<VISCREFRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = ViscrefTable(pdata);
}

void unpack(WatdentTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<WATDENTRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = WatdentTable(pdata);
}

void unpack(PolyInjTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<double> throughputs, velocities;
    int tableNumber;
    std::vector<std::vector<double>> tableData;
    unpack(throughputs, buffer, position, comm);
    unpack(velocities, buffer, position, comm);
    unpack(tableNumber, buffer, position, comm);
    unpack(tableData, buffer, position, comm);
    data = PolyInjTable(throughputs, velocities, tableNumber, tableData);
}

void unpack(PlymwinjTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(static_cast<PolyInjTable&>(data), buffer, position, comm);
}

void unpack(SkprpolyTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(static_cast<PolyInjTable&>(data), buffer, position, comm);
    double refConcentration;
    unpack(refConcentration, buffer, position, comm);
    data.setReferenceConcentration(refConcentration);
}

void unpack(SkprwatTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(static_cast<PolyInjTable&>(data), buffer, position, comm);
}

void unpack(RockTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<ROCKRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = RockTable(pdata);
}

void unpack(TableManager& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    std::map<std::string, TableContainer> simpleTables;
    SplitSimpleTables split;
    std::vector<PvtgTable> pvtgTables;
    std::vector<PvtoTable> pvtoTables;
    std::vector<Rock2dTable> rock2dTables;
    std::vector<Rock2dtrTable> rock2dtrTables;
    PvtwTable pvtwTable;
    PvcdoTable pvcdoTable;
    DensityTable densityTable;
    PlyvmhTable plyvmhTable;
    RockTable rockTable;
    ViscrefTable viscrefTable;
    PlmixparTable plmixparTable;
    ShrateTable shrateTable;
    Stone1exTable stone1exTable;
    TlmixparTable tlmixparTable;
    WatdentTable watdentTable;
    std::vector<PvtwsaltTable> pvtwsaltTables;
    std::vector<BrineDensityTable> bdensityTables;
    std::vector<SolventDensityTable> sdensityTables;
    std::map<int, PlymwinjTable> plymwinjTables;
    std::map<int, SkprwatTable> skprwatTables;
    std::map<int, SkprpolyTable> skprpolyTables;
    Tabdims tabdims;
    Regdims regdims;
    Eqldims eqldims;
    Aqudims aqudims;
    bool hasImptvd;
    bool hasEntpvd;
    bool hasEqlnum;
    bool hasShrate;
    DenT oilDenT, gasDenT, watDenT;
    StandardCond stcond;
    std::size_t gas_comp_index;
    std::shared_ptr<JFunc> jfunc;
    double rtemp;

    unpack(simpleTables, buffer, position, comm);
    unpack(split.plyshMax, buffer, position, comm);
    unpack(split.plyshMap, buffer, position, comm);
    unpack(split.rockMax, buffer, position, comm);
    unpack(split.rockMap, buffer, position, comm);
    unpack(pvtgTables, buffer, position, comm);
    unpack(pvtoTables, buffer, position, comm);
    unpack(rock2dTables, buffer, position, comm);
    unpack(rock2dtrTables, buffer, position, comm);
    unpack(pvtwTable, buffer, position, comm);
    unpack(pvcdoTable, buffer, position, comm);
    unpack(densityTable, buffer, position, comm);
    unpack(plyvmhTable, buffer, position, comm);
    unpack(rockTable, buffer, position, comm);
    unpack(plmixparTable, buffer, position, comm);
    unpack(shrateTable, buffer, position, comm);
    unpack(stone1exTable, buffer, position, comm);
    unpack(tlmixparTable, buffer, position, comm);
    unpack(viscrefTable, buffer, position, comm);
    unpack(watdentTable, buffer, position, comm);
    unpack(pvtwsaltTables, buffer, position, comm);
    unpack(bdensityTables, buffer, position, comm);
    unpack(sdensityTables, buffer, position, comm);
    unpack(plymwinjTables, buffer, position, comm);
    unpack(skprwatTables, buffer, position, comm);
    unpack(skprpolyTables, buffer, position, comm);
    unpack(tabdims, buffer, position, comm);
    unpack(regdims, buffer, position, comm);
    unpack(eqldims, buffer, position, comm);
    unpack(aqudims, buffer, position, comm);
    unpack(hasImptvd, buffer, position, comm);
    unpack(hasEntpvd, buffer, position, comm);
    unpack(hasEqlnum, buffer, position, comm);
    unpack(hasShrate, buffer, position, comm);
    bool hasJf;
    unpack(hasJf, buffer, position, comm);
    if (hasJf) {
        jfunc = std::make_shared<JFunc>();
        unpack(*jfunc, buffer, position, comm);
    }
    unpack(oilDenT, buffer, position, comm);
    unpack(gasDenT, buffer, position, comm);
    unpack(watDenT, buffer, position, comm);
    unpack(stcond, buffer, position, comm);
    unpack(gas_comp_index, buffer, position, comm);
    unpack(rtemp, buffer, position, comm);

    if (split.plyshMax > 0) {
        TableContainer container(split.plyshMax);
        for (const auto& it : split.plyshMap) {
            container.addTable(it.first, it.second);
        }
        simpleTables.insert(std::make_pair("PLYSHLOG", container));
    }
    if (split.rockMax > 0) {
        TableContainer container(split.rockMax);
        for (const auto& it : split.rockMap) {
            container.addTable(it.first, it.second);
        }
        simpleTables.insert(std::make_pair("ROCKTAB", container));
    }

    data = TableManager(simpleTables, pvtgTables, pvtoTables, rock2dTables,
                        rock2dtrTables, pvtwTable, pvcdoTable, densityTable,
                        plyvmhTable, rockTable, plmixparTable, shrateTable, stone1exTable,
                        tlmixparTable, viscrefTable, watdentTable, pvtwsaltTables,
                        bdensityTables, sdensityTables, plymwinjTables, skprwatTables,
                        skprpolyTables, tabdims, regdims, eqldims, aqudims, hasImptvd,
                        hasEntpvd, hasEqlnum, hasShrate, jfunc, oilDenT, gasDenT,
                        watDenT, stcond, gas_comp_index, rtemp);
}

void unpack(OilVaporizationProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    OilVaporizationProperties::OilVaporization type;
    double vap1, vap2;
    std::vector<double> maxDRSDT, maxDRVDT;
    std::vector<bool> maxDRSDT_allCells;
    unpack(type, buffer, position, comm);
    unpack(vap1, buffer, position, comm);
    unpack(vap2, buffer, position, comm);
    unpack(maxDRSDT, buffer, position, comm);
    unpack(maxDRSDT_allCells, buffer, position, comm);
    unpack(maxDRVDT, buffer, position, comm);
    data = OilVaporizationProperties(type, vap1, vap2, maxDRSDT,
                                     maxDRSDT_allCells, maxDRVDT);
}

void unpack(Events& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    DynamicVector<uint64_t> events;
    unpack(events, buffer, position, comm);
    data = Events(events);
}

void unpack(MessageLimits& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    DynamicState<MLimits> limits;
    unpack(limits, buffer, position, comm);
    data = MessageLimits(limits);
}

void unpack(VFPInjTable& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int tableNum;
    double datumDepth;
    VFPInjTable::FLO_TYPE floType;
    std::vector<double> floAxis, thpAxis;
    VFPInjTable::array_type table;

    unpack(tableNum, buffer, position, comm);
    unpack(datumDepth, buffer, position, comm);
    unpack(floType, buffer, position, comm);
    unpack(floAxis, buffer, position, comm);
    unpack(thpAxis, buffer, position, comm);
    unpack(table, buffer, position, comm);

    data = VFPInjTable(tableNum, datumDepth, floType,
                       floAxis, thpAxis, table);
}

void unpack(VFPProdTable& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int tableNum;
    double datumDepth;
    VFPProdTable::FLO_TYPE floType;
    VFPProdTable::WFR_TYPE wfrType;
    VFPProdTable::GFR_TYPE gfrType;
    VFPProdTable::ALQ_TYPE alqType;
    std::vector<double> floAxis, thpAxis, wfrAxis, gfrAxis, alqAxis;
    VFPProdTable::array_type table;

    unpack(tableNum, buffer, position, comm);
    unpack(datumDepth, buffer, position, comm);
    unpack(floType, buffer, position, comm);
    unpack(wfrType, buffer, position, comm);
    unpack(gfrType, buffer, position, comm);
    unpack(alqType, buffer, position, comm);
    unpack(floAxis, buffer, position, comm);
    unpack(thpAxis, buffer, position, comm);
    unpack(wfrAxis, buffer, position, comm);
    unpack(gfrAxis, buffer, position, comm);
    unpack(alqAxis, buffer, position, comm);
    unpack(table, buffer, position, comm);

    data = VFPProdTable(tableNum, datumDepth, floType, wfrType,
                        gfrType, alqType, floAxis, thpAxis,
                        wfrAxis, gfrAxis, alqAxis, table);
}

void unpack(WellTestConfig::WTESTWell& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.name, buffer, position, comm);
    unpack(data.shut_reason, buffer, position, comm);
    unpack(data.test_interval, buffer, position, comm);
    unpack(data.num_test, buffer, position, comm);
    unpack(data.startup_time, buffer, position, comm);
    unpack(data.begin_report_step, buffer, position, comm);
}

void unpack(WellTestConfig& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<WellTestConfig::WTESTWell> ddata;
    unpack(ddata, buffer, position, comm);
    data = WellTestConfig(ddata);
}

void unpack(WellTracerProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    WellTracerProperties::ConcentrationMap ddata;
    unpack(ddata, buffer, position, comm);
    data = WellTracerProperties(ddata);
}

void unpack(UDAValue& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    bool isDouble;
    Dimension dim;
    unpack(dim, buffer, position, comm);
    unpack(isDouble, buffer, position, comm);
    if (isDouble) {
        double val;
        unpack(val, buffer, position, comm);
        data = UDAValue(val, dim);
    } else {
        std::string val;
        unpack(val, buffer, position, comm);
        data = UDAValue(val, dim);
    }
}

void unpack(Connection& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    Connection::Direction dir;
    double depth;
    Connection::State state;
    int satTableId, complnum;
    double CF, Kh, rw, r0, skinFactor;
    int I, J, K;
    size_t seqIndex;
    double segDistStart, segDistEnd;
    bool defaultSatTabId;
    size_t compSegSeqIndex;
    int segment;
    double wellPi;
    Connection::CTFKind kind;

    unpack(dir, buffer, position, comm);
    unpack(depth, buffer, position, comm);
    unpack(state, buffer, position, comm);
    unpack(satTableId, buffer, position, comm);
    unpack(complnum, buffer, position, comm);
    unpack(CF, buffer, position, comm);
    unpack(Kh, buffer, position, comm);
    unpack(rw, buffer, position, comm);
    unpack(r0, buffer, position, comm);
    unpack(skinFactor, buffer, position, comm);
    unpack(I, buffer, position, comm);
    unpack(J, buffer, position, comm);
    unpack(K, buffer, position, comm);
    unpack(kind, buffer, position, comm);
    unpack(seqIndex, buffer, position, comm);
    unpack(segDistStart, buffer, position, comm);
    unpack(segDistEnd, buffer, position, comm);
    unpack(defaultSatTabId, buffer, position, comm);
    unpack(compSegSeqIndex, buffer, position, comm);
    unpack(segment, buffer, position, comm);
    unpack(wellPi, buffer, position, comm);

    data = Connection(dir, depth, state, satTableId,
                      complnum, CF, Kh, rw, r0,
                      skinFactor, {I,J,K}, kind, seqIndex,
                      segDistStart, segDistEnd,
                      defaultSatTabId, compSegSeqIndex,
                      segment, wellPi);
}

void unpack(Well::WellInjectionProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.name, buffer, position, comm);
    unpack(data.surfaceInjectionRate, buffer, position, comm);
    unpack(data.reservoirInjectionRate, buffer, position, comm);
    unpack(data.BHPTarget, buffer, position, comm);
    unpack(data.THPTarget, buffer, position, comm);
    unpack(data.bhp_hist_limit, buffer, position, comm);
    unpack(data.thp_hist_limit, buffer, position, comm);
    unpack(data.temperature, buffer, position, comm);
    unpack(data.BHPH, buffer, position, comm);
    unpack(data.THPH, buffer, position, comm);
    unpack(data.VFPTableNumber, buffer, position, comm);
    unpack(data.predictionMode, buffer, position, comm);
    unpack(data.injectionControls, buffer, position, comm);
    unpack(data.injectorType, buffer, position, comm);
    unpack(data.controlMode, buffer, position, comm);
}

void unpack(WellEconProductionLimits& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double minOilRate, minGasRate, maxWaterCut, maxGasOilRatio, maxWaterGasRatio;
    WellEconProductionLimits::EconWorkover workover, workoverSecondary;
    bool endRun;
    std::string followonWell;
    WellEconProductionLimits::QuantityLimit quantityLimit;
    double secondaryMaxWaterCut, maxGasLiquidRatio, minLiquidRate,
           maxTemperature, minReservoirFluidRate;
    unpack(minOilRate, buffer, position, comm);
    unpack(minGasRate, buffer, position, comm);
    unpack(maxWaterCut, buffer, position, comm);
    unpack(maxGasOilRatio, buffer, position, comm);
    unpack(maxWaterGasRatio, buffer, position, comm);
    unpack(workover, buffer, position, comm);
    unpack(endRun, buffer, position, comm);
    unpack(followonWell, buffer, position, comm);
    unpack(quantityLimit, buffer, position, comm);
    unpack(secondaryMaxWaterCut, buffer, position, comm);
    unpack(workoverSecondary, buffer, position, comm);
    unpack(maxGasLiquidRatio, buffer, position, comm);
    unpack(minLiquidRate, buffer, position, comm);
    unpack(maxTemperature, buffer, position, comm);
    unpack(minReservoirFluidRate, buffer, position, comm);
    data = WellEconProductionLimits(minOilRate, minGasRate, maxWaterCut,
                                    maxGasOilRatio, maxWaterGasRatio,
                                    workover, endRun, followonWell,
                                    quantityLimit, secondaryMaxWaterCut,
                                    workoverSecondary, maxGasLiquidRatio,
                                    minLiquidRate, maxTemperature,
                                    minReservoirFluidRate);
}

void unpack(WellConnections& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int headI, headJ;
    size_t numRemoved;
    std::vector<Connection> connections;

    unpack(headI, buffer, position, comm),
    unpack(headJ, buffer, position, comm),
    unpack(numRemoved, buffer, position, comm),
    unpack(connections, buffer, position, comm),

    data = WellConnections(headI, headJ, numRemoved, connections);
}

void unpack(Well::WellProductionProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    UDAValue OilRate, WaterRate, GasRate, LiquidRate, ResVRate;
    UDAValue BHPTarget, THPTarget;
    double bhp_hist_limit, thp_hist_limit;
    double BHPH, THPH;
    int VFPTableNumber;
    double ALQValue;
    bool predictionMode;
    Well::ProducerCMode controlMode, whistctl_cmode;
    int prodCtrls;

    unpack(name, buffer, position, comm);
    unpack(OilRate, buffer, position, comm);
    unpack(WaterRate, buffer, position, comm);
    unpack(GasRate, buffer, position, comm);
    unpack(LiquidRate, buffer, position, comm);
    unpack(ResVRate, buffer, position, comm);
    unpack(BHPTarget, buffer, position, comm);
    unpack(THPTarget, buffer, position, comm);
    unpack(bhp_hist_limit, buffer, position, comm);
    unpack(thp_hist_limit, buffer, position, comm);
    unpack(BHPH, buffer, position, comm);
    unpack(THPH, buffer, position, comm);
    unpack(VFPTableNumber, buffer, position, comm);
    unpack(ALQValue, buffer, position, comm);
    unpack(predictionMode, buffer, position, comm);
    unpack(controlMode, buffer, position, comm);
    unpack(whistctl_cmode, buffer, position, comm);
    unpack(prodCtrls, buffer, position, comm);
    data = Well::WellProductionProperties(name, OilRate, WaterRate, GasRate,
                                          LiquidRate, ResVRate, BHPTarget,
                                          THPTarget, bhp_hist_limit, thp_hist_limit,
                                          BHPH, THPH, VFPTableNumber,
                                          ALQValue, predictionMode, controlMode,
                                          whistctl_cmode, prodCtrls);
}

void unpack(SpiralICD& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double strength, length, densityCalibration,
           viscosityCalibration, criticalValue,
           widthTransitionRegion, maxViscosityRatio;
    int methodFlowScaling;
    double maxAbsoluteRate;
    ICDStatus status;
    double scalingFactor;

    unpack(strength, buffer, position, comm);
    unpack(length, buffer, position, comm);
    unpack(densityCalibration, buffer, position, comm);
    unpack(viscosityCalibration, buffer, position, comm);
    unpack(criticalValue, buffer, position, comm);
    unpack(widthTransitionRegion, buffer, position, comm);
    unpack(maxViscosityRatio, buffer, position, comm);
    unpack(methodFlowScaling, buffer, position, comm);
    unpack(maxAbsoluteRate, buffer, position, comm);
    unpack(status, buffer, position, comm);
    unpack(scalingFactor, buffer, position, comm);

    data = SpiralICD(strength, length, densityCalibration,
                     viscosityCalibration, criticalValue,
                     widthTransitionRegion, maxViscosityRatio,
                     methodFlowScaling, maxAbsoluteRate,
                     status, scalingFactor);
}

void unpack(Valve& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double conEFlowCoefficient;
    double conCrossArea;
    double conMaxCrossArea;
    double pipeAdditionalLength;
    double pipeDiameter;
    double pipeRoughness;
    double pipeCrossArea;
    ICDStatus status;

    unpack(conEFlowCoefficient, buffer, position, comm);
    unpack(conCrossArea, buffer, position, comm);
    unpack(conMaxCrossArea, buffer, position, comm);
    unpack(pipeAdditionalLength, buffer, position, comm);
    unpack(pipeDiameter, buffer, position, comm);
    unpack(pipeRoughness, buffer, position, comm);
    unpack(pipeCrossArea, buffer, position, comm);
    unpack(status, buffer, position, comm);
    data = Valve(conEFlowCoefficient, conCrossArea, conMaxCrossArea,
                 pipeAdditionalLength, pipeDiameter, pipeRoughness,
                 pipeCrossArea, status);
}

void unpack(Segment& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int segmentNumber, branchNumber, outletSegment;
    std::vector<int> inletSegments;
    double totalLength, depth, internalDiameter, roughness, crossArea, volume;
    bool dataReady;
    Segment::SegmentType segmentType;
    std::shared_ptr<SpiralICD> spiralICD;
    std::shared_ptr<Valve> valve;

    unpack(segmentNumber, buffer, position, comm);
    unpack(branchNumber, buffer, position, comm);
    unpack(outletSegment, buffer, position, comm);
    unpack(inletSegments, buffer, position, comm);
    unpack(totalLength, buffer, position, comm);
    unpack(depth, buffer, position, comm);
    unpack(internalDiameter, buffer, position, comm);
    unpack(roughness, buffer, position, comm);
    unpack(crossArea, buffer, position, comm);
    unpack(volume, buffer, position, comm);
    unpack(dataReady, buffer, position, comm);
    unpack(segmentType, buffer, position, comm);
    unpack(spiralICD, buffer, position, comm);
    unpack(valve, buffer, position, comm);
    data = Segment(segmentNumber, branchNumber, outletSegment,
                   inletSegments, totalLength, depth,
                   internalDiameter, roughness, crossArea,
                   volume, dataReady, segmentType, spiralICD, valve);
}

template<class T>
void unpack(std::shared_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    bool hasVal;
    unpack(hasVal, buffer, position, comm);
    if (hasVal) {
        data = std::make_shared<T>();
        unpack(*data, buffer, position, comm);
    }
}

template<class T>
void unpack(std::unique_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    bool hasVal;
    unpack(hasVal, buffer, position, comm);
    if (hasVal) {
        data.reset(new T);
        unpack(*data, buffer, position, comm);
    }
}

void unpack(Dimension& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double siScaling, siOffset;

    unpack(siScaling, buffer, position, comm);
    unpack(siOffset, buffer, position, comm);
    data = Dimension(siScaling, siOffset);
}

void unpack(UnitSystem& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    UnitSystem::UnitType type;
    std::map<std::string, Dimension> dimensions;
    size_t use_count;
    unpack(name, buffer, position, comm);
    unpack(type, buffer, position, comm);
    unpack(dimensions, buffer, position, comm);
    unpack(use_count, buffer, position, comm);

    data = UnitSystem(name, type, dimensions, use_count);
}

void unpack(WellSegments& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    WellSegments::CompPressureDrop compPressureDrop;
    std::vector<Segment> segments;

    unpack(compPressureDrop, buffer, position, comm);
    unpack(segments, buffer, position, comm);

    data = WellSegments(compPressureDrop, segments);
}

void unpack(Well& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name, groupName;
    std::size_t firstTimeStep, seqIndex;
    int headI, headJ;
    double ref_depth;
    WellType wtype;
    Connection::Order ordering;
    UnitSystem units;
    double udq_undefined;
    Well::Status status;
    double drainageRadius;
    bool allowCrossFlow, automaticShutIn;
    Well::WellGuideRate guideRate;
    double efficiencyFactor;
    double solventFraction;
    bool prediction_mode;
    auto econLimits = std::make_shared<WellEconProductionLimits>();
    auto foamProperties = std::make_shared<WellFoamProperties>();
    auto polymerProperties = std::make_shared<WellPolymerProperties>();
    auto brineProperties = std::make_shared<WellBrineProperties>();
    auto tracerProperties = std::make_shared<WellTracerProperties>();
    auto connection = std::make_shared<WellConnections>();
    auto production = std::make_shared<Well::WellProductionProperties>();
    auto injection = std::make_shared<Well::WellInjectionProperties>();
    std::shared_ptr<WellSegments> segments;

    unpack(name, buffer, position, comm);
    unpack(groupName, buffer, position, comm);
    unpack(firstTimeStep, buffer, position, comm);
    unpack(seqIndex, buffer, position, comm);
    unpack(headI, buffer, position, comm);
    unpack(headJ, buffer, position, comm);
    unpack(ref_depth, buffer, position, comm);
    unpack(wtype, buffer, position, comm);
    unpack(ordering, buffer, position, comm);
    unpack(units, buffer, position, comm);
    unpack(udq_undefined, buffer, position, comm);
    unpack(status, buffer, position, comm);
    unpack(drainageRadius, buffer, position, comm);
    unpack(allowCrossFlow, buffer, position, comm);
    unpack(automaticShutIn, buffer, position, comm);
    unpack(guideRate, buffer, position, comm);
    unpack(efficiencyFactor, buffer, position, comm);
    unpack(solventFraction, buffer, position, comm);
    unpack(prediction_mode, buffer, position, comm);
    unpack(*econLimits, buffer, position, comm);
    unpack(*foamProperties, buffer, position, comm);
    unpack(*polymerProperties, buffer, position, comm);
    unpack(*brineProperties, buffer, position, comm);
    unpack(*tracerProperties, buffer, position, comm);
    unpack(*connection, buffer, position, comm);
    unpack(*production, buffer, position, comm);
    unpack(*injection, buffer, position, comm);
    bool hasSegments;
    unpack(hasSegments, buffer, position, comm);
    if (hasSegments) {
        segments = std::make_shared<WellSegments>();
        unpack(*segments, buffer, position, comm);
    }
    data = Well(name, groupName, firstTimeStep, seqIndex, headI, headJ,
                ref_depth, wtype, ordering, units, udq_undefined, status,
                drainageRadius, allowCrossFlow, automaticShutIn,
                guideRate, efficiencyFactor, solventFraction, prediction_mode,
                econLimits, foamProperties, polymerProperties, brineProperties,
                tracerProperties, connection, production, injection, segments);
}

template<class T>
void unpack(IOrderSet<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    typename IOrderSet<T>::index_type index;
    typename IOrderSet<T>::storage_type storage;
    unpack(index, buffer, position, comm);
    unpack(storage, buffer, position, comm);
    data = IOrderSet<T>(index, storage);
}

template void unpack(std::map<Phase,Group::GroupInjectionProperties>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

void unpack(Group::GroupInjectionProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.phase, buffer, position, comm);
    unpack(data.cmode, buffer, position, comm);
    unpack(data.surface_max_rate, buffer, position, comm);
    unpack(data.resv_max_rate, buffer, position, comm);
    unpack(data.target_reinj_fraction, buffer, position, comm);
    unpack(data.target_void_fraction, buffer, position, comm);
    unpack(data.reinj_group, buffer, position, comm);
    unpack(data.voidage_group, buffer, position, comm);
    unpack(data.injection_controls, buffer, position, comm);
}

void unpack(Group::GroupProductionProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.cmode, buffer, position, comm);
    unpack(data.exceed_action, buffer, position, comm);
    unpack(data.oil_target, buffer, position, comm);
    unpack(data.water_target, buffer, position, comm);
    unpack(data.gas_target, buffer, position, comm);
    unpack(data.liquid_target, buffer, position, comm);
    unpack(data.guide_rate, buffer, position, comm);
    unpack(data.guide_rate_def, buffer, position, comm);
    unpack(data.resv_target, buffer, position, comm);
    unpack(data.production_controls, buffer, position, comm);
}

void unpack(Group& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    std::size_t insert_index, initStep;
    double udqUndefined;
    UnitSystem units;
    Group::GroupType type;
    double groupEfficiencyFactor;
    bool transferGroupEfficiencyFactor;
    bool availableForGroupControl;
    int groupNetVFPTable;
    std::string parent;
    IOrderSet<std::string> wells, groups;
    std::map<Phase, Group::GroupInjectionProperties> injection;
    Group::GroupProductionProperties production;

    unpack(name, buffer, position, comm);
    unpack(insert_index, buffer, position, comm);
    unpack(initStep, buffer, position, comm);
    unpack(udqUndefined, buffer, position, comm);
    unpack(units, buffer, position, comm);
    unpack(type, buffer, position, comm);
    unpack(groupEfficiencyFactor, buffer, position, comm);
    unpack(transferGroupEfficiencyFactor, buffer, position, comm);
    unpack(availableForGroupControl, buffer, position, comm);
    unpack(groupNetVFPTable, buffer, position, comm);
    unpack(parent, buffer, position, comm);
    unpack(wells, buffer, position, comm);
    unpack(groups, buffer, position, comm);
    unpack(injection, buffer, position, comm);
    unpack(production, buffer, position, comm);
    data = Group(name, insert_index, initStep, udqUndefined,
                 units, type, groupEfficiencyFactor,
                 transferGroupEfficiencyFactor,
                 availableForGroupControl,
                 groupNetVFPTable, parent, wells, groups,
                 injection, production);
}

void unpack(WList& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    WList::storage ddata;
    unpack(ddata, buffer, position, comm);
    data = WList(ddata);
}

void unpack(WListManager& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::map<std::string,WList> lists;
    unpack(lists, buffer, position, comm);
    data = WListManager(lists);
}

void unpack(UDQASTNode& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    UDQVarType var_type;
    UDQTokenType type;
    std::string stringValue;
    double scalarValue;
    std::vector<std::string> selectors;
    std::shared_ptr<UDQASTNode> left;
    std::shared_ptr<UDQASTNode> right;

    unpack(var_type, buffer, position, comm);
    unpack(type, buffer, position, comm);
    unpack(stringValue, buffer, position, comm);
    unpack(scalarValue, buffer, position, comm);
    unpack(selectors, buffer, position, comm);
    unpack(left, buffer, position, comm);
    unpack(right, buffer, position, comm);
    data = UDQASTNode(var_type, type, stringValue, scalarValue, selectors, left, right);
}

void unpack(UDQDefine& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string keyword;
    std::shared_ptr<UDQASTNode> ast;
    UDQVarType varType;
    std::string string_data;

    unpack(keyword, buffer, position, comm);
    unpack(ast, buffer, position, comm);
    unpack(varType, buffer, position, comm);
    unpack(string_data, buffer, position, comm);
    data = UDQDefine(keyword, ast, varType, string_data);
}

void unpack(UDQAssign::AssignRecord& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.selector, buffer, position, comm);
    unpack(data.value, buffer, position, comm);
}

void unpack(UDQAssign& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string keyword;
    UDQVarType varType;
    std::vector<UDQAssign::AssignRecord> records;
    unpack(keyword, buffer, position, comm);
    unpack(varType, buffer, position, comm);
    unpack(records, buffer, position, comm);
    data = UDQAssign(keyword, varType, records);
}

void unpack(UDQIndex& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.insert_index, buffer, position, comm);
    unpack(data.typed_insert_index, buffer, position, comm);
    unpack(data.action, buffer, position, comm);
    unpack(data.var_type, buffer, position, comm);
}

void unpack(UDQConfig& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    UDQParams params;
    UDQFunctionTable function_table;
    std::unordered_map<std::string,UDQDefine> definitionsMap;
    std::unordered_map<std::string,UDQAssign> assignmentsMap;
    std::unordered_map<std::string,std::string> units;
    OrderedMap<std::string,UDQIndex> inputIndex;
    std::map<UDQVarType,std::size_t> typeCount;

    unpack(params, buffer, position, comm);
    function_table = UDQFunctionTable(params);
    unpack(definitionsMap, buffer, position, comm);
    unpack(assignmentsMap, buffer, position, comm);
    unpack(units, buffer, position, comm);
    unpack(inputIndex, buffer, position, comm);
    unpack(typeCount, buffer, position, comm);
    data = UDQConfig(params, function_table, definitionsMap,
                     assignmentsMap, units, inputIndex, typeCount);
}

void unpack(UDQActive::InputRecord& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.input_index, buffer, position, comm);
    unpack(data.udq, buffer, position, comm);
    unpack(data.wgname, buffer, position, comm);
    unpack(data.control, buffer, position, comm);
}

void unpack(UDQActive::Record& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.udq, buffer, position, comm);
    unpack(data.input_index, buffer, position, comm);
    unpack(data.use_index, buffer, position, comm);
    unpack(data.wgname, buffer, position, comm);
    unpack(data.control, buffer, position, comm);
    unpack(data.uad_code, buffer, position, comm);
    unpack(data.use_count, buffer, position, comm);
}

void unpack(UDQActive& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<UDQActive::InputRecord> inputRecords;
    std::vector<UDQActive::Record> outputRecords;
    std::unordered_map<std::string,std::size_t> udqKeys, wgKeys;

    unpack(inputRecords, buffer, position, comm);
    unpack(outputRecords, buffer, position, comm);
    unpack(udqKeys, buffer, position, comm);
    unpack(wgKeys, buffer, position, comm);
    data = UDQActive(inputRecords, outputRecords, udqKeys, wgKeys);
}

void unpack(GuideRateModel& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double timeInterval;
    GuideRateModel::Target target;
    std::array<double,6> coefs;
    bool allow_increase, free_gas, defaultModel;
    double damping_factor;
    std::array<UDAValue,3> udaCoefs;

    unpack(timeInterval, buffer, position, comm);
    unpack(target, buffer, position, comm);
    unpack(coefs, buffer, position, comm);
    unpack(allow_increase, buffer, position, comm);
    unpack(damping_factor, buffer, position, comm);
    unpack(free_gas, buffer, position, comm);
    unpack(defaultModel, buffer, position, comm);
    unpack(udaCoefs, buffer, position, comm);
    data = GuideRateModel(timeInterval, target, coefs, allow_increase,
                          damping_factor, free_gas, defaultModel, udaCoefs);
}

void unpack(GuideRateConfig& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::shared_ptr<GuideRateModel> model;
    std::unordered_map<std::string, GuideRateConfig::WellTarget> wells;
    std::unordered_map<std::string, GuideRateConfig::GroupTarget> groups;

    unpack(model, buffer, position, comm);
    unpack(wells, buffer, position, comm);
    unpack(groups, buffer, position, comm);
    data = GuideRateConfig(model, wells, groups);
}

void unpack(GConSale::GCONSALEGroup& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.sales_target, buffer, position, comm);
    unpack(data.max_sales_rate, buffer, position, comm);
    unpack(data.min_sales_rate, buffer, position, comm);
    unpack(data.max_proc, buffer, position, comm);
    unpack(data.udq_undefined, buffer, position, comm);
    unpack(data.unit_system, buffer, position, comm);
}

void unpack(GConSale& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::map<std::string,GConSale::GCONSALEGroup> groups;
    unpack(groups, buffer, position, comm);
    data = GConSale(groups);
}

void unpack(GConSump::GCONSUMPGroup& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.consumption_rate, buffer, position, comm);
    unpack(data.import_rate, buffer, position, comm);
    unpack(data.network_node, buffer, position, comm);
    unpack(data.udq_undefined, buffer, position, comm);
    unpack(data.unit_system, buffer, position, comm);
}

void unpack(GConSump& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::map<std::string,GConSump::GCONSUMPGroup> groups;
    unpack(groups, buffer, position, comm);
    data = GConSump(groups);
}

void unpack(RFTConfig& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TimeMap timeMap;
    std::size_t first_rft;
    std::pair<bool, std::size_t> wellOpenRftTime;
    RFTConfig::WellOpenTimeMap wellOpenRftName;
    RFTConfig::WellOpenTimeMap wellOpen;
    RFTConfig::RFTMap rftConfig;
    RFTConfig::PLTMap pltConfig;

    unpack(timeMap, buffer, position, comm);
    unpack(first_rft, buffer, position, comm);
    unpack(wellOpenRftTime, buffer, position, comm);
    unpack(wellOpenRftName, buffer, position, comm);
    unpack(wellOpen, buffer, position, comm);
    unpack(rftConfig, buffer, position, comm);
    unpack(pltConfig, buffer, position, comm);
    data = RFTConfig(timeMap, first_rft, wellOpenRftTime, wellOpenRftName,
                     wellOpen, rftConfig, pltConfig);
}

void unpack(DeckItem& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<double> dVal;
    std::vector<int> iVal;
    std::vector<std::string> sVal;
    std::vector<UDAValue> uVal;
    type_tag type;
    std::string name;
    std::vector<value::status> valueStatus;
    bool rawData;
    std::vector<Dimension> activeDimensions, defaultDimensions;

    unpack(dVal, buffer, position, comm);
    unpack(iVal, buffer, position, comm);
    unpack(sVal, buffer, position, comm);
    unpack(uVal, buffer, position, comm);
    unpack(type, buffer, position, comm);
    unpack(name, buffer, position, comm);
    unpack(valueStatus, buffer, position, comm);
    unpack(rawData, buffer, position, comm);
    unpack(activeDimensions, buffer, position, comm);
    unpack(defaultDimensions, buffer, position, comm);
    data = DeckItem(dVal, iVal, sVal, uVal, type, name,
                    valueStatus, rawData, activeDimensions, defaultDimensions);
}

void unpack(DeckRecord& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<DeckItem> items;
    unpack(items, buffer, position, comm);
    data = DeckRecord(std::move(items));
}

void unpack(Location& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    data.filename.clear();
    unpack(data.filename, buffer, position, comm);
    unpack(data.lineno, buffer, position, comm);
}

void unpack(DeckKeyword& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    Location location;
    std::vector<DeckRecord> records;
    bool isDataKeyword, isSlashTerminated;

    unpack(name, buffer, position, comm);
    unpack(location, buffer, position, comm);
    unpack(records, buffer, position, comm);
    unpack(isDataKeyword, buffer, position, comm);
    unpack(isSlashTerminated, buffer, position, comm);
    data = DeckKeyword(name, location, records,
                       isDataKeyword, isSlashTerminated);
}

void unpack(Deck& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<DeckKeyword> keywords;
    UnitSystem defaultUnitSystem;
    std::unique_ptr<UnitSystem> activeUnitSystem;
    std::string dataFile, inputPath;
    size_t accessCount;

    unpack(keywords, buffer, position, comm);
    unpack(defaultUnitSystem, buffer, position, comm);
    unpack(activeUnitSystem, buffer, position, comm);
    unpack(dataFile, buffer, position, comm);
    unpack(inputPath, buffer, position, comm);
    unpack(accessCount, buffer, position, comm);
    data = Deck(keywords, defaultUnitSystem,
                activeUnitSystem.get(), dataFile, inputPath, accessCount);
}

void unpack(Action::ASTNode& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TokenType token;
    FuncType func_type;
    std::string func;
    std::vector<std::string> argList;
    double number;
    std::vector<Action::ASTNode> children;

    unpack(token, buffer, position, comm);
    unpack(func_type, buffer, position, comm);
    unpack(func, buffer, position, comm);
    unpack(argList, buffer, position, comm);
    unpack(number, buffer, position, comm);
    unpack(children, buffer, position, comm);
    data = Action::ASTNode(token, func_type, func, argList, number, children);
}

void unpack(Action::AST& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::shared_ptr<Action::ASTNode> condition;
    unpack(condition, buffer, position, comm);
    data = Action::AST(condition);
}

void unpack(Action::Quantity& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.quantity, buffer, position, comm);
    unpack(data.args, buffer, position, comm);
}

void unpack(Action::Condition& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.lhs, buffer, position, comm);
    unpack(data.rhs, buffer, position, comm);
    unpack(data.logic, buffer, position, comm);
    unpack(data.cmp, buffer, position, comm);
    unpack(data.cmp_string, buffer, position, comm);
}

void unpack(Action::ActionX& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    size_t max_run;
    double min_wait;
    std::time_t start_time;
    std::vector<DeckKeyword> keywords;
    Action::AST condition;
    std::vector<Action::Condition> conditions;
    size_t run_count;
    std::time_t last_run;

    unpack(name, buffer, position, comm);
    unpack(max_run, buffer, position, comm);
    unpack(min_wait, buffer, position, comm);
    unpack(start_time, buffer, position, comm);
    unpack(keywords, buffer, position, comm);
    unpack(condition, buffer, position, comm);
    unpack(conditions, buffer, position, comm);
    unpack(run_count, buffer, position, comm);
    unpack(last_run, buffer, position, comm);
    data = Action::ActionX(name, max_run, min_wait, start_time, keywords,
                           condition, conditions, run_count, last_run);
}

void unpack(Action::Actions& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Action::ActionX> actions;
    unpack(actions, buffer, position, comm);
    data = Action::Actions(actions);
}

void unpack(Schedule& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TimeMap timeMap;
    Schedule::WellMap staticWells;
    Schedule::GroupMap groups;
    DynamicState<OilVaporizationProperties> oilVapProps;
    Events events;
    DynamicVector<Deck> modifierDeck;
    DynamicState<Tuning> tuning;
    MessageLimits messageLimits;
    Runspec runspec;
    Schedule::VFPProdMap vfpProdTables;
    Schedule::VFPInjMap vfpInjTables;
    DynamicState<std::shared_ptr<WellTestConfig>> wellTestConfig;
    DynamicState<std::shared_ptr<WListManager>> wListManager;
    DynamicState<std::shared_ptr<UDQConfig>> udqConfig;
    DynamicState<std::shared_ptr<UDQActive>> udqActive;
    DynamicState<std::shared_ptr<GuideRateConfig>> guideRateConfig;
    DynamicState<std::shared_ptr<GConSale>> gconSale;
    DynamicState<std::shared_ptr<GConSump>> gconSump;
    DynamicState<Well::ProducerCMode> globalWhistCtlMode;
    DynamicState<std::shared_ptr<Action::Actions>> actions;
    RFTConfig rftConfig;
    DynamicState<int> nupCol;
    RestartConfig restartConfig;
    std::map<std::string,Events> wellGroupEvents;

    unpack(timeMap, buffer, position, comm);
    unpackDynMap(staticWells, buffer, position, comm);
    unpackDynMap(groups, buffer, position, comm);
    unpack(oilVapProps, buffer, position, comm);
    unpack(events, buffer, position, comm);
    unpack(modifierDeck, buffer, position, comm);
    unpack(tuning, buffer, position, comm);
    unpack(messageLimits, buffer, position, comm);
    unpack(runspec, buffer, position, comm);
    unpackDynMap<Map2>(vfpProdTables, buffer, position, comm);
    unpackDynMap<Map2>(vfpInjTables, buffer, position, comm);
    unpack(wellTestConfig, buffer, position, comm);
    unpack(wListManager, buffer, position, comm);
    unpack(udqConfig, buffer, position, comm);
    unpack(udqActive, buffer, position, comm);
    unpack(guideRateConfig, buffer, position, comm);
    unpack(gconSale, buffer, position, comm);
    unpack(gconSump, buffer, position, comm);
    unpack(globalWhistCtlMode, buffer, position, comm);
    unpack(actions, buffer, position, comm);

    unpack(rftConfig, buffer, position, comm);
    unpack(nupCol, buffer, position, comm);
    unpack(restartConfig, buffer, position, comm);
    unpack(wellGroupEvents, buffer, position, comm);
    data = Schedule(timeMap, staticWells, groups, oilVapProps, events,
                    modifierDeck, tuning, messageLimits, runspec,
                    vfpProdTables, vfpInjTables, wellTestConfig,
                    wListManager, udqConfig, udqActive, guideRateConfig,
                    gconSale, gconSump, globalWhistCtlMode, actions,
                    rftConfig, nupCol, restartConfig, wellGroupEvents);
}

void unpack(BrineDensityTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<double> tableValues;

    unpack(tableValues, buffer, position, comm);
    data = BrineDensityTable(tableValues);
}

void unpack(PvtwsaltTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double refPressValue, refSaltConValue;
    std::vector<double> tableValues;

    unpack(refPressValue, buffer, position, comm);
    unpack(refSaltConValue, buffer, position, comm);
    unpack(tableValues, buffer, position, comm);
    data = PvtwsaltTable(refPressValue, refSaltConValue, tableValues);
}

void unpack(EquilRecord& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double datumDepth, datumDepthPressure, waterOilContactDepth;
    double waterOilContactCapillaryPressure, gasOilContactDepth;
    double gasOilContactCapillaryPressure;
    bool liveOilInitConstantRs, wetGasInitConstantRv;
    int initializationTargetAccuracy;

    unpack(datumDepth, buffer, position, comm);
    unpack(datumDepthPressure, buffer, position, comm);
    unpack(waterOilContactDepth, buffer, position, comm);
    unpack(waterOilContactCapillaryPressure, buffer, position, comm);
    unpack(gasOilContactDepth, buffer, position, comm);
    unpack(gasOilContactCapillaryPressure, buffer, position, comm);
    unpack(liveOilInitConstantRs, buffer, position, comm);
    unpack(wetGasInitConstantRv, buffer, position, comm);
    unpack(initializationTargetAccuracy, buffer, position, comm);
    data = EquilRecord(datumDepth, datumDepthPressure, waterOilContactDepth,
                       waterOilContactCapillaryPressure, gasOilContactDepth,
                       gasOilContactCapillaryPressure, liveOilInitConstantRs,
                       wetGasInitConstantRv, initializationTargetAccuracy);
}

void unpack(FoamData& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double referenceSurfactantConcentration, exponent;
    double minimumSurfactantConcentration;
    bool allowDesorption;
    double rockDensity;

    unpack(referenceSurfactantConcentration, buffer, position, comm);
    unpack(exponent, buffer, position, comm);
    unpack(minimumSurfactantConcentration, buffer, position, comm);
    unpack(allowDesorption, buffer, position, comm);
    unpack(rockDensity, buffer, position, comm);
    data = FoamData(referenceSurfactantConcentration, exponent,
                    minimumSurfactantConcentration, allowDesorption, rockDensity);
}

void unpack(RestartSchedule& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.timestep, buffer, position, comm);
    unpack(data.basic, buffer, position, comm);
    unpack(data.frequency, buffer, position, comm);
    unpack(data.rptsched_restart_set, buffer, position, comm);
    unpack(data.rptsched_restart, buffer, position, comm);
}

void unpack(TimeStampUTC& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TimeStampUTC::YMD ymd;
    int hour, minutes, seconds, usec;

    unpack(ymd, buffer, position, comm);
    unpack(hour, buffer, position, comm);
    unpack(minutes, buffer, position, comm);
    unpack(seconds, buffer, position, comm);
    unpack(usec, buffer, position, comm);
    data = TimeStampUTC(ymd, hour, minutes, seconds, usec);
}

void unpack(EclHysterConfig& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    bool active;
    int pcHysteresisModel, krHysteresisModel;

    unpack(active, buffer, position, comm);
    unpack(pcHysteresisModel, buffer, position, comm);
    unpack(krHysteresisModel, buffer, position, comm);
    data = EclHysterConfig(active, pcHysteresisModel, krHysteresisModel);
}

void unpack(JFunc& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    JFunc::Flag flag;
    double owSurfaceTension, goSurfaceTension;
    double alphaFactor, betaFactor;
    JFunc::Direction dir;

    unpack(flag, buffer, position, comm);
    unpack(owSurfaceTension, buffer, position, comm);
    unpack(goSurfaceTension, buffer, position, comm);
    unpack(alphaFactor, buffer, position, comm);
    unpack(betaFactor, buffer, position, comm);
    unpack(dir, buffer, position, comm);
    data = JFunc(flag, owSurfaceTension, goSurfaceTension,
                  alphaFactor, betaFactor, dir);
}

void unpack(WellPolymerProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.m_polymerConcentration, buffer, position, comm);
    unpack(data.m_saltConcentration, buffer, position, comm);
    unpack(data.m_plymwinjtable, buffer, position, comm);
    unpack(data.m_skprwattable, buffer, position, comm);
    unpack(data.m_skprpolytable, buffer, position, comm);
}

void unpack(Well::WellGuideRate& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.available, buffer, position, comm);
    unpack(data.guide_rate, buffer, position, comm);
    unpack(data.guide_phase, buffer, position, comm);
    unpack(data.scale_factor, buffer, position, comm);
}

void unpack(GuideRateConfig::WellTarget& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.guide_rate, buffer, position, comm);
    unpack(data.target, buffer, position, comm);
    unpack(data.scaling_factor, buffer, position, comm);
}

void unpack(GuideRateConfig::GroupTarget& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.guide_rate, buffer, position, comm);
    unpack(data.target, buffer, position, comm);
}

void unpack(MULTREGTRecord& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.src_value, buffer, position, comm);
    unpack(data.target_value, buffer, position, comm);
    unpack(data.trans_mult, buffer, position, comm);
    unpack(data.directions, buffer, position, comm);
    unpack(data.nnc_behaviour, buffer, position, comm);
    unpack(data.region_name, buffer, position, comm);
}

void unpack(MULTREGTScanner& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::array<size_t, 3> size;
    std::vector<MULTREGTRecord> records;
    MULTREGTScanner::ExternalSearchMap searchMap;
    std::map<std::string, std::vector<int>> regions;
    std::string defaultRegion;

    unpack(size, buffer, position, comm);
    unpack(records, buffer, position, comm);
    unpack(searchMap, buffer, position, comm);
    unpack(regions, buffer, position, comm);
    unpack(defaultRegion, buffer, position, comm);

    data = MULTREGTScanner(size, records, searchMap, regions, defaultRegion);
}

void unpack(EclipseConfig& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    InitConfig init;
    IOConfig io;

    unpack(init, buffer, position, comm);
    unpack(io, buffer, position, comm);
    data = EclipseConfig(init, io);
}

void unpack(TransMult& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::array<size_t, 3> size;
    std::map<FaceDir::DirEnum, std::vector<double>> trans;
    std::map<FaceDir::DirEnum, std::string> names;
    MULTREGTScanner scanner;

    unpack(size, buffer, position, comm);
    unpack(trans, buffer, position, comm);
    unpack(names, buffer, position, comm);
    unpack(scanner, buffer, position, comm);
    data = TransMult(size, trans, names, scanner);
}

void unpack(FaultFace& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<size_t> indices;
    FaceDir::DirEnum dir;

    unpack(indices, buffer, position, comm);
    unpack(dir, buffer, position, comm);
    data = FaultFace(indices, dir);
}

void unpack(Fault& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    double transMult;
    std::vector<FaultFace> faceList;

    unpack(name, buffer, position, comm);
    unpack(transMult, buffer, position, comm);
    unpack(faceList, buffer, position, comm);
    data = Fault(name, transMult, faceList);
}

void unpack(FaultCollection& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    OrderedMap<std::string, Fault> faults;

    unpack(faults, buffer, position, comm);
    data = FaultCollection(faults);
}

void unpack(SolventDensityTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<double> tableValues;

    unpack(tableValues, buffer, position, comm);
    data = SolventDensityTable(tableValues);
}

void unpack(GridDims& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::array<int,3> NXYZ;

    unpack(NXYZ, buffer, position, comm);
    data = GridDims(NXYZ);
}

void unpack(ShrateTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<ShrateRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = ShrateTable(pdata);
}

void unpack(TlmixparTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<TlmixparRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = TlmixparTable(pdata);
}

void unpack(PlmixparTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<PlmixparRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = PlmixparTable(pdata);
}

void unpack(PlyvmhTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<PlyvmhRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = PlyvmhTable(pdata);
}

void unpack(Stone1exTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Stone1exRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = Stone1exTable(pdata);
}

void unpack(PlyshlogTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TableSchema schema;
    OrderedMap<std::string, TableColumn> columns;
    bool jfunc;
    double refPolymerConcentration;
    double refSalinity;
    double refTemperature;
    bool hasRefSalinity;
    bool hasRefTemperature;

    unpack(schema, buffer, position, comm);
    unpack(columns, buffer, position, comm);
    unpack(jfunc, buffer, position, comm);
    unpack(refPolymerConcentration, buffer, position, comm);
    unpack(refSalinity, buffer, position, comm);
    unpack(refTemperature, buffer, position, comm);
    unpack(hasRefSalinity, buffer, position, comm);
    unpack(hasRefTemperature, buffer, position, comm);

    data = PlyshlogTable(schema, columns, jfunc,
                         refPolymerConcentration, refSalinity,
                         refTemperature, hasRefSalinity, hasRefTemperature);
}

void unpack(RocktabTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TableSchema schema;
    OrderedMap<std::string, TableColumn> columns;
    bool jfunc;
    bool isDirectional;

    unpack(schema, buffer, position, comm);
    unpack(columns, buffer, position, comm);
    unpack(jfunc, buffer, position, comm);
    unpack(isDirectional, buffer, position, comm);

    data = RocktabTable(schema, columns, jfunc, isDirectional);
}

#define INSTANTIATE_PACK_VECTOR(...) \
template std::size_t packSize(const std::vector<__VA_ARGS__>& data, \
                              Dune::MPIHelper::MPICommunicator comm); \
template void pack(const std::vector<__VA_ARGS__>& data, \
                   std::vector<char>& buffer, int& position, \
                   Dune::MPIHelper::MPICommunicator comm); \
template void unpack(std::vector<__VA_ARGS__>& data, \
                     std::vector<char>& buffer, int& position, \
                     Dune::MPIHelper::MPICommunicator comm);

INSTANTIATE_PACK_VECTOR(double)
INSTANTIATE_PACK_VECTOR(std::vector<double>)
INSTANTIATE_PACK_VECTOR(bool)
INSTANTIATE_PACK_VECTOR(std::string)
INSTANTIATE_PACK_VECTOR(char)
INSTANTIATE_PACK_VECTOR(int)
INSTANTIATE_PACK_VECTOR(std::array<double, 3>)

#undef INSTANTIATE_PACK_VECTOR

#define INSTANTIATE_PACK_SET(...) \
template std::size_t packSize(const std::set<__VA_ARGS__>& data, \
                              Dune::MPIHelper::MPICommunicator comm); \
template void pack(const std::set<__VA_ARGS__>& data, \
                   std::vector<char>& buffer, int& position, \
                   Dune::MPIHelper::MPICommunicator comm); \
template void unpack(std::set<__VA_ARGS__>& data, \
                     std::vector<char>& buffer, int& position, \
                     Dune::MPIHelper::MPICommunicator comm);

INSTANTIATE_PACK_SET(std::string)

#undef INSTANTIATE_PACK_SET

#define INSTANTIATE_PACK_SHARED_PTR(...) \
template std::size_t packSize(const std::shared_ptr<__VA_ARGS__>& data, \
                              Dune::MPIHelper::MPICommunicator comm); \
template void pack(const std::shared_ptr<__VA_ARGS__>& data, \
                   std::vector<char>& buffer, int& position, \
                   Dune::MPIHelper::MPICommunicator comm); \
template void unpack(std::shared_ptr<__VA_ARGS__>& data, \
                     std::vector<char>& buffer, int& position, \
                     Dune::MPIHelper::MPICommunicator comm);

INSTANTIATE_PACK_SHARED_PTR(SpiralICD)
#undef INSTANTIATE_PACK_SHARED_PTR

#define INSTANTIATE_PACK(...) \
template std::size_t packSize(const __VA_ARGS__& data, \
                              Dune::MPIHelper::MPICommunicator comm); \
template void pack(const __VA_ARGS__& data, \
                   std::vector<char>& buffer, int& position, \
                   Dune::MPIHelper::MPICommunicator comm); \
template void unpack(__VA_ARGS__& data, \
                     std::vector<char>& buffer, int& position, \
                     Dune::MPIHelper::MPICommunicator comm);

INSTANTIATE_PACK(double)
INSTANTIATE_PACK(std::size_t)
INSTANTIATE_PACK(bool)
INSTANTIATE_PACK(int)
INSTANTIATE_PACK(std::array<short,3>)
INSTANTIATE_PACK(std::array<bool,3>)
INSTANTIATE_PACK(unsigned char)
#undef INSTANTIATE_PACK

} // end namespace Mpi

RestartValue loadParallelRestart(const EclipseIO* eclIO, SummaryState& summaryState,
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
        restartValues = eclIO->loadRestart(summaryState, solutionKeys, extraKeys);
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
