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

#include <ewoms/common/tabulated1dfunction.hh>
#include <ewoms/common/intervaltabulated2dfunction.hh>
#include <ewoms/common/uniformxtabulated2dfunction.hh>
#include <ewoms/eclio/output/restartvalue.hh>
#include <ewoms/eclio/output/eclipseio.hh>
#include <ewoms/eclio/output/summary.hh>
#include <ewoms/eclio/parser/eclipsestate/aquiferconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/aquancon.hh>
#include <ewoms/eclio/parser/eclipsestate/aquiferct.hh>
#include <ewoms/eclio/parser/eclipsestate/aquifetp.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/dynamicstate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/dynamicvector.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/gconsale.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/gconsump.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/group.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/guiderateconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/timemap.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqassign.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqactive.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/well.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/welltestconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/bcconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/rockconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/dent.hh>
#include <ewoms/eclio/parser/eclipsestate/util/orderedmap.hh>

#include <dune/common/parallel/mpihelper.hh>

#include <set>
#include <tuple>
#include <vector>
#include <map>
#include <unordered_map>

namespace Ewoms
{

class Actdims;

namespace Action {
    class Actions;
    class ActionX;
    class AST;
    class ASTNode;
    class Condition;
    class Quantity;
}

class Aqudims;
class BrineDensityTable;
class ColumnSchema;
class Connection;
class DeckItem;
class DeckRecord;
struct DENSITYRecord;
struct DensityTable;
class Dimension;
class EclHysterConfig;
class EclipseConfig;
class Eqldims;
template<class Scalar> struct EclEpsScalingPointsInfo;
class EDITNNC;
class EndpointScaling;
class Equil;
class EquilRecord;
class Events;
class Fault;
class FaultCollection;
class FaultFace;
class FoamConfig;
class FoamData;
class GridDims;
class InitConfig;
class IOConfig;
template<class T> class IOrderSet;
class JFunc;
class Location;
class MessageLimits;
struct MLimits;
struct MULTREGTRecord;
class MULTREGTScanner;
class NNC;
struct NNCdata;
class OilVaporizationProperties;
class Phases;
class PlymwinjTable;
class PlyshlogTable;
class PlyvmhRecord;
class PlyvmhTable;
class PolyInjTable;
class PVCDORecord;
class PvcdoTable;
class PlmixparRecord;
class PlmixparTable;
class PvtgTable;
class PvtoTable;
class PVTWRecord;
class PvtwsaltTable;
class PvtwTable;
class Regdims;
class RestartConfig;
class RestartSchedule;
class RFTConfig;
class ROCKRecord;
class RockTable;
class RocktabTable;
class Rock2dTable;
class Rock2dtrTable;
class Runspec;
class Schedule;
class Segment;
class ShrateRecord;
class ShrateTable;
class SimulationConfig;
class SimpleTable;
class SkprpolyTable;
class SkprwatTable;
class SolventDensityTable;
class SpiralICD;
class StandardCond;
class Stone1exRecord;
class Stone1exTable;
class SummaryConfig;
class SummaryNode;
class Tabdims;
class TableColumn;
class TableContainer;
class TableManager;
class TableSchema;
class ThresholdPressure;
class TimeStampUTC;
class TlmixparRecord;
class TlmixparTable;
class TransMult;
struct Tuning;
class UDAValue;
class UDQASTNode;
class UDQConfig;
class UDQDefine;
class UDQIndex;
class UDQParams;
class UnitSystem;
class Valve;
class VFPInjTable;
class VFPProdTable;
struct VISCREFRecord;
struct ViscrefTable;
struct WATDENTRecord;
struct WatdentTable;
class WellConnections;
class Welldims;
class WellEconProductionLimits;
struct WellFoamProperties;
struct WellPolymerProperties;
class WellSegmentDims;
class WellSegments;
class WellTracerProperties;
class WList;
class WListManager;

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
                     std::integral_constant<bool, false>);

template<class T>
std::size_t packSize(const T&, Dune::MPIHelper::MPICommunicator comm,
                     std::integral_constant<bool, true>);

template<class T>
std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, Dune::MPIHelper::MPICommunicator comm);

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

template<class T>
std::size_t packSize(const std::shared_ptr<T>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class T, std::size_t N>
std::size_t packSize(const std::array<T,N>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const std::unique_ptr<T>& data,
                     Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const char* str, Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const std::string& str, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
std::size_t packSize(const std::map<T1,T2,C,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
std::size_t packSize(const std::unordered_map<T1,T2,H,P,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Key, class Value>
std::size_t packSize(const OrderedMap<Key,Value>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const DynamicVector<T>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const DynamicState<T>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const Tabulated1DFunction<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const IntervalTabulated2DFunction<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const UniformXTabulated2DFunction<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const IOrderSet<T>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const EclEpsScalingPointsInfo<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm);

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
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>);

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2>
void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

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

template<class T>
void pack(const std::shared_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T, size_t N>
void pack(const std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const std::unique_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
void pack(const std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
void pack(const std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class Key, class Value>
void pack(const OrderedMap<Key,Value>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const DynamicState<T>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const DynamicVector<T>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const Tabulated1DFunction<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const IntervalTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const UniformXTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const IOrderSet<T>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const EclEpsScalingPointsInfo<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

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
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>);

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2>
void unpack(std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

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

template<class T>
void unpack(std::shared_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T, size_t N>
void unpack(std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(std::unique_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
void unpack(std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
void unpack(std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class Key, class Value>
void unpack(OrderedMap<Key,Value>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(DynamicState<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(DynamicVector<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(Tabulated1DFunction<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(IntervalTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(UniformXTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(IOrderSet<T>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(EclEpsScalingPointsInfo<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

/// prototypes for complex types

#define ADD_PACK_PROTOTYPES(T) \
  std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm); \
  void pack(const T& data, std::vector<char>& buffer, int& position, \
          Dune::MPIHelper::MPICommunicator comm); \
  void unpack(T& data, std::vector<char>& buffer, int& position, \
              Dune::MPIHelper::MPICommunicator comm);

ADD_PACK_PROTOTYPES(Actdims)
ADD_PACK_PROTOTYPES(Action::Actions)
ADD_PACK_PROTOTYPES(Action::ActionX)
ADD_PACK_PROTOTYPES(Action::AST)
ADD_PACK_PROTOTYPES(Action::ASTNode)
ADD_PACK_PROTOTYPES(Action::Condition)
ADD_PACK_PROTOTYPES(Action::Quantity)
ADD_PACK_PROTOTYPES(Aqudims)
ADD_PACK_PROTOTYPES(AquiferConfig)
ADD_PACK_PROTOTYPES(Aquancon)
ADD_PACK_PROTOTYPES(Aquancon::AquancCell)
ADD_PACK_PROTOTYPES(AquiferCT)
ADD_PACK_PROTOTYPES(AquiferCT::AQUCT_data)
ADD_PACK_PROTOTYPES(Aquifetp)
ADD_PACK_PROTOTYPES(Aquifetp::AQUFETP_data)
ADD_PACK_PROTOTYPES(BCConfig)
ADD_PACK_PROTOTYPES(BCConfig::BCFace)
ADD_PACK_PROTOTYPES(BrineDensityTable)
ADD_PACK_PROTOTYPES(ColumnSchema)
ADD_PACK_PROTOTYPES(Connection)
ADD_PACK_PROTOTYPES(data::CellData)
ADD_PACK_PROTOTYPES(data::Connection)
ADD_PACK_PROTOTYPES(data::CurrentControl)
ADD_PACK_PROTOTYPES(data::Rates)
ADD_PACK_PROTOTYPES(data::Segment)
ADD_PACK_PROTOTYPES(data::Solution)
ADD_PACK_PROTOTYPES(data::Well)
ADD_PACK_PROTOTYPES(data::WellRates)
ADD_PACK_PROTOTYPES(Deck)
ADD_PACK_PROTOTYPES(DeckItem)
ADD_PACK_PROTOTYPES(DeckKeyword)
ADD_PACK_PROTOTYPES(DeckRecord)
ADD_PACK_PROTOTYPES(DENSITYRecord)
ADD_PACK_PROTOTYPES(DensityTable)
ADD_PACK_PROTOTYPES(DenT)
ADD_PACK_PROTOTYPES(DenT::entry)
ADD_PACK_PROTOTYPES(Dimension)
ADD_PACK_PROTOTYPES(EclHysterConfig)
ADD_PACK_PROTOTYPES(EclipseConfig)
ADD_PACK_PROTOTYPES(EDITNNC)
ADD_PACK_PROTOTYPES(EndpointScaling)
ADD_PACK_PROTOTYPES(Equil)
ADD_PACK_PROTOTYPES(Eqldims)
ADD_PACK_PROTOTYPES(EquilRecord)
ADD_PACK_PROTOTYPES(Events)
ADD_PACK_PROTOTYPES(Fault)
ADD_PACK_PROTOTYPES(FaultCollection)
ADD_PACK_PROTOTYPES(FaultFace)
ADD_PACK_PROTOTYPES(FoamConfig)
ADD_PACK_PROTOTYPES(FoamData)
ADD_PACK_PROTOTYPES(GConSale)
ADD_PACK_PROTOTYPES(GConSale::GCONSALEGroup)
ADD_PACK_PROTOTYPES(GConSump)
ADD_PACK_PROTOTYPES(GConSump::GCONSUMPGroup)
ADD_PACK_PROTOTYPES(GridDims)
ADD_PACK_PROTOTYPES(GuideRateConfig)
ADD_PACK_PROTOTYPES(GuideRateConfig::GroupTarget)
ADD_PACK_PROTOTYPES(GuideRateConfig::WellTarget)
ADD_PACK_PROTOTYPES(GuideRateModel)
ADD_PACK_PROTOTYPES(Group)
ADD_PACK_PROTOTYPES(Group::GroupInjectionProperties)
ADD_PACK_PROTOTYPES(Group::GroupProductionProperties)
ADD_PACK_PROTOTYPES(InitConfig)
ADD_PACK_PROTOTYPES(IOConfig)
ADD_PACK_PROTOTYPES(JFunc)
ADD_PACK_PROTOTYPES(Location)
ADD_PACK_PROTOTYPES(MessageLimits)
ADD_PACK_PROTOTYPES(MLimits)
ADD_PACK_PROTOTYPES(MULTREGTRecord)
ADD_PACK_PROTOTYPES(MULTREGTScanner)
ADD_PACK_PROTOTYPES(NNC)
ADD_PACK_PROTOTYPES(NNCdata)
ADD_PACK_PROTOTYPES(OilVaporizationProperties)
ADD_PACK_PROTOTYPES(Phases)
ADD_PACK_PROTOTYPES(PlmixparRecord)
ADD_PACK_PROTOTYPES(PlmixparTable)
ADD_PACK_PROTOTYPES(PlymwinjTable)
ADD_PACK_PROTOTYPES(PlyshlogTable)
ADD_PACK_PROTOTYPES(PlyvmhRecord)
ADD_PACK_PROTOTYPES(PlyvmhTable)
ADD_PACK_PROTOTYPES(PolyInjTable)
ADD_PACK_PROTOTYPES(PVCDORecord)
ADD_PACK_PROTOTYPES(PvcdoTable)
ADD_PACK_PROTOTYPES(PvtgTable)
ADD_PACK_PROTOTYPES(PvtoTable)
ADD_PACK_PROTOTYPES(PVTWRecord)
ADD_PACK_PROTOTYPES(PvtwsaltTable)
ADD_PACK_PROTOTYPES(PvtwTable)
ADD_PACK_PROTOTYPES(Regdims)
ADD_PACK_PROTOTYPES(RestartConfig)
ADD_PACK_PROTOTYPES(RestartKey)
ADD_PACK_PROTOTYPES(RestartSchedule)
ADD_PACK_PROTOTYPES(RestartValue)
ADD_PACK_PROTOTYPES(RFTConfig)
ADD_PACK_PROTOTYPES(RockConfig)
ADD_PACK_PROTOTYPES(RockConfig::RockComp)
ADD_PACK_PROTOTYPES(ROCKRecord)
ADD_PACK_PROTOTYPES(RockTable)
ADD_PACK_PROTOTYPES(Rock2dTable)
ADD_PACK_PROTOTYPES(Rock2dtrTable)
ADD_PACK_PROTOTYPES(RocktabTable)
ADD_PACK_PROTOTYPES(Runspec)
ADD_PACK_PROTOTYPES(Schedule)
ADD_PACK_PROTOTYPES(Segment)
ADD_PACK_PROTOTYPES(ShrateRecord)
ADD_PACK_PROTOTYPES(ShrateTable)
ADD_PACK_PROTOTYPES(SimulationConfig)
ADD_PACK_PROTOTYPES(SimpleTable)
ADD_PACK_PROTOTYPES(SkprpolyTable)
ADD_PACK_PROTOTYPES(SkprwatTable)
ADD_PACK_PROTOTYPES(SolventDensityTable)
ADD_PACK_PROTOTYPES(SpiralICD)
ADD_PACK_PROTOTYPES(std::string)
ADD_PACK_PROTOTYPES(Stone1exRecord)
ADD_PACK_PROTOTYPES(Stone1exTable)
ADD_PACK_PROTOTYPES(SummaryConfig)
ADD_PACK_PROTOTYPES(SummaryNode)
ADD_PACK_PROTOTYPES(Tabdims)
ADD_PACK_PROTOTYPES(TableColumn)
ADD_PACK_PROTOTYPES(TableContainer)
ADD_PACK_PROTOTYPES(TableManager)
ADD_PACK_PROTOTYPES(TableSchema)
ADD_PACK_PROTOTYPES(ThresholdPressure)
ADD_PACK_PROTOTYPES(TimeMap)
ADD_PACK_PROTOTYPES(TimeStampUTC)
ADD_PACK_PROTOTYPES(TlmixparRecord)
ADD_PACK_PROTOTYPES(TlmixparTable)
ADD_PACK_PROTOTYPES(TransMult)
ADD_PACK_PROTOTYPES(Tuning)
ADD_PACK_PROTOTYPES(UDAValue)
ADD_PACK_PROTOTYPES(UDQActive)
ADD_PACK_PROTOTYPES(UDQActive::InputRecord)
ADD_PACK_PROTOTYPES(UDQActive::Record)
ADD_PACK_PROTOTYPES(UDQAssign)
ADD_PACK_PROTOTYPES(UDQAssign::AssignRecord)
ADD_PACK_PROTOTYPES(UDQASTNode)
ADD_PACK_PROTOTYPES(UDQConfig)
ADD_PACK_PROTOTYPES(UDQDefine)
ADD_PACK_PROTOTYPES(UDQIndex)
ADD_PACK_PROTOTYPES(UDQParams)
ADD_PACK_PROTOTYPES(UnitSystem)
ADD_PACK_PROTOTYPES(Valve)
ADD_PACK_PROTOTYPES(VFPInjTable)
ADD_PACK_PROTOTYPES(VFPProdTable)
ADD_PACK_PROTOTYPES(VISCREFRecord)
ADD_PACK_PROTOTYPES(ViscrefTable)
ADD_PACK_PROTOTYPES(WATDENTRecord)
ADD_PACK_PROTOTYPES(WatdentTable)
ADD_PACK_PROTOTYPES(Well)
ADD_PACK_PROTOTYPES(Well::WellGuideRate)
ADD_PACK_PROTOTYPES(Well::WellInjectionProperties)
ADD_PACK_PROTOTYPES(Well::WellProductionProperties)
ADD_PACK_PROTOTYPES(WellBrineProperties)
ADD_PACK_PROTOTYPES(WellConnections)
ADD_PACK_PROTOTYPES(Welldims)
ADD_PACK_PROTOTYPES(WellEconProductionLimits)
ADD_PACK_PROTOTYPES(WellFoamProperties)
ADD_PACK_PROTOTYPES(WellPolymerProperties)
ADD_PACK_PROTOTYPES(WellSegmentDims)
ADD_PACK_PROTOTYPES(WellSegments)
ADD_PACK_PROTOTYPES(WellTestConfig)
ADD_PACK_PROTOTYPES(WellTestConfig::WTESTWell)
ADD_PACK_PROTOTYPES(WellTracerProperties)
ADD_PACK_PROTOTYPES(WList)
ADD_PACK_PROTOTYPES(WListManager)

template<class T, class C>
const T& packAndSend(const T& in, const C& comm)
{
    if (comm.size() == 1)
        return in;

    std::size_t size = packSize(in, comm);
    std::vector<char> buffer(size);
    int pos = 0;
    Mpi::pack(in, buffer, pos, comm);
    comm.broadcast(&pos, 1, 0);
    comm.broadcast(buffer.data(), pos, 0);
    return in;
}

template<class T, class C>
void receiveAndUnpack(T& result, const C& comm)
{
    int size;
    comm.broadcast(&size, 1, 0);
    std::vector<char> buffer(size);
    comm.broadcast(buffer.data(), size, 0);
    int pos = 0;
    unpack(result, buffer, pos, comm);
}
} // end namespace Mpi

RestartValue loadParallelRestart(const EclipseIO* eclIO, SummaryState& summaryState,
                                 const std::vector<Ewoms::RestartKey>& solutionKeys,
                                 const std::vector<Ewoms::RestartKey>& extraKeys,
                                 Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm);

} // end namespace Ewoms
#endif // PARALLEL_RESTART_HH
