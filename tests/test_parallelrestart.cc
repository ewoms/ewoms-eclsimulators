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

#define BOOST_TEST_MODULE TestParallelRestart
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <tuple>
#include <utility>

#include <ewoms/eclio/opmlog/keywordlocation.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/deck/deckitem.hh>
#include <ewoms/eclio/parser/eclipsestate/aquancon.hh>
#include <ewoms/eclio/parser/eclipsestate/aquiferct.hh>
#include <ewoms/eclio/parser/eclipsestate/aquiferconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/aquifetp.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipseconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/runspec.hh>
#include <ewoms/eclio/parser/eclipsestate/tracerconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/edit/editnnc.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/facedir.hh>
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
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actionast.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/pyaction.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actions.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actionx.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/astnode.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/condition.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/events.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/gconsale.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/group.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/guideratemodel.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/messagelimits.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/icd.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/sicd.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/valve.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/network/node.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/oilvaporizationproperties.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/rftconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/scheduletypes.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/timemap.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/tuning.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqactive.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqassign.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqastnode.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqdefine.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqfunction.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqfunctiontable.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqinput.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/vfpinjtable.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/vfpprodtable.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/connection.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wellfoamproperties.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wellpolymerproperties.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/welltracerproperties.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/well.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wlist.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wlistmanager.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/bcconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/rockconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/simulationconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/thresholdpressure.hh>
#include <ewoms/eclio/parser/eclipsestate/summaryconfig/summaryconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/aqudims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/columnschema.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/eqldims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/flattable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/dent.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/jfunc.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/plymwinjtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/plyshlogtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/pvtgtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/pvtotable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/regdims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rock2dtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rock2dtrtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rocktabtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/simpletable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/skprpolytable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/skprwattable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tabdims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablecolumn.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablecontainer.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablemanager.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tableschema.hh>
#include <ewoms/eclio/output/restartvalue.hh>
#include <ewoms/eclsimulators/utils/parallelrestart.hh>
#include <eebos/eclmpiserializer.hh>

namespace {

Ewoms::data::Solution getSolution()
{
    Ewoms::data::Solution sol1;
    sol1.insert("testdata", Ewoms::UnitSystem::measure::length,
                {1.0, 2.0, 3.0}, Ewoms::data::TargetType::RESTART_SOLUTION);
    return sol1;
}

Ewoms::data::Rates getRates()
{
    Ewoms::data::Rates rat1;
    rat1.set(Ewoms::data::Rates::opt::wat, 1.0);
    rat1.set(Ewoms::data::Rates::opt::oil, 2.0);
    rat1.set(Ewoms::data::Rates::opt::gas, 3.0);
    rat1.set(Ewoms::data::Rates::opt::polymer, 4.0);
    rat1.set(Ewoms::data::Rates::opt::solvent, 5.0);
    rat1.set(Ewoms::data::Rates::opt::energy, 6.0);
    rat1.set(Ewoms::data::Rates::opt::dissolved_gas, 7.0);
    rat1.set(Ewoms::data::Rates::opt::vaporized_oil, 8.0);
    rat1.set(Ewoms::data::Rates::opt::reservoir_water, 9.0);
    rat1.set(Ewoms::data::Rates::opt::reservoir_oil, 10.0);
    rat1.set(Ewoms::data::Rates::opt::reservoir_gas, 11.0);
    rat1.set(Ewoms::data::Rates::opt::productivity_index_water, 12.0);
    rat1.set(Ewoms::data::Rates::opt::productivity_index_oil, 13.0);
    rat1.set(Ewoms::data::Rates::opt::productivity_index_gas, 14.0);
    rat1.set(Ewoms::data::Rates::opt::well_potential_water, 15.0);
    rat1.set(Ewoms::data::Rates::opt::well_potential_oil, 16.0);
    rat1.set(Ewoms::data::Rates::opt::well_potential_gas, 17.0);

    return rat1;
}

Ewoms::data::Connection getConnection()
{
    Ewoms::data::Connection con1;
    con1.rates = getRates();
    con1.index = 1;
    con1.pressure = 2.0;
    con1.reservoir_rate = 3.0;
    con1.cell_pressure = 4.0;
    con1.cell_saturation_water = 5.0;
    con1.cell_saturation_gas = 6.0;
    con1.effective_Kh = 7.0;
    return con1;
}

Ewoms::data::Segment getSegment()
{
    Ewoms::data::Segment seg1;
    seg1.rates = getRates();
    seg1.segNumber = 1;
    const auto pres_idx = Ewoms::data::SegmentPressures::Value::Pressure;
    seg1.pressures[pres_idx] = 2.0;
    return seg1;
}

Ewoms::data::CurrentControl getCurrentControl()
{
    Ewoms::data::CurrentControl curr;
    curr.isProducer = true;
    curr.prod = ::Ewoms::Well::ProducerCMode::CRAT;
    return curr;
}

Ewoms::data::GuideRateValue getWellGuideRate()
{
    using Item = Ewoms::data::GuideRateValue::Item;

    return Ewoms::data::GuideRateValue{}.set(Item::Oil  , 1.23)
                                      .set(Item::Gas  , 2.34)
                                      .set(Item::Water, 3.45)
                                      .set(Item::ResV , 4.56);
}

Ewoms::data::Well getWell()
{
    Ewoms::data::Well well1;
    well1.rates = getRates();
    well1.bhp = 1.0;
    well1.thp = 2.0;
    well1.temperature = 3.0;
    well1.control = 4;
    well1.connections.push_back(getConnection());
    well1.segments.insert({0, getSegment()});
    well1.current_control = getCurrentControl();
    well1.guide_rates = getWellGuideRate();
    return well1;
}

Ewoms::data::GroupGuideRates getGroupGuideRates()
{
    using Item = Ewoms::data::GuideRateValue::Item;

    auto gr = Ewoms::data::GroupGuideRates{};

    gr.production.set(Item::Oil ,   999.888)
                 .set(Item::Gas ,  8888.777)
                 .set(Item::ResV, 12345.678);

    gr.injection.set(Item::Gas  , 9876.543)
                .set(Item::Water, 2345.987);

    return gr;
}

Ewoms::data::GroupConstraints getGroupConstraints()
{
    using PMode = ::Ewoms::Group::ProductionCMode;
    using IMode = ::Ewoms::Group::InjectionCMode;

    return Ewoms::data::GroupConstraints{}
    .set(PMode::ORAT,           // Production
         IMode::VREP,           // Gas Injection
         IMode::NONE);          // Water Injection
}

Ewoms::data::GroupData getGroupData()
{
    return Ewoms::data::GroupData {
        getGroupConstraints(),
        getGroupGuideRates()
    };
}
}

template<class T>
std::tuple<T,int,int> PackUnpack(const T& in)
{
    auto comm = Dune::MPIHelper::getCollectiveCommunication();
    std::size_t packSize = Ewoms::Mpi::packSize(in, comm);
    std::vector<char> buffer(packSize);
    int pos1 = 0;
    Ewoms::Mpi::pack(in, buffer, pos1, comm);
    int pos2 = 0;
    T out;
    Ewoms::Mpi::unpack(out, buffer, pos2, comm);

    return std::make_tuple(out, pos1, pos2);
}

template<class T>
std::tuple<T,int,int> PackUnpack2(T& in)
{
    auto comm = Dune::MPIHelper::getCollectiveCommunication();
    Ewoms::EclMpiSerializer ser(comm);
    ser.pack(in);
    size_t pos1 = ser.position();
    T out;
    ser.unpack(out);
    size_t pos2 = ser.position();

    return std::make_tuple(out, pos1, pos2);
}

#define DO_CHECKS(TYPE_NAME) \
    BOOST_CHECK_MESSAGE(std::get<1>(val2) == std::get<2>(val2), "Packed size differ from unpack size for " #TYPE_NAME);  \
    BOOST_CHECK_MESSAGE(val1 == std::get<0>(val2), "Deserialized " #TYPE_NAME " differ");

BOOST_AUTO_TEST_CASE(Solution)
{
    Ewoms::data::Solution val1 = getSolution();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Solution)
}

BOOST_AUTO_TEST_CASE(Rates)
{
    Ewoms::data::Rates val1 = getRates();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Rates)
}

BOOST_AUTO_TEST_CASE(dataGuideRateValue)
{
    using Item = Ewoms::data::GuideRateValue::Item;

    const auto val1 = Ewoms::data::GuideRateValue{}
    .set(Item::Oil ,   999.888)
    .set(Item::Gas ,  8888.777)
    .set(Item::ResV, 12345.678);

    const auto val2 = PackUnpack(val1);

    BOOST_CHECK_MESSAGE(! std::get<0>(val2).has(Item::Water),
                        "Water Must Not Appear From "
                        "Serializing GuideRateValues");

    DO_CHECKS(data::GuideRateValue)
}

BOOST_AUTO_TEST_CASE(dataConnection)
{
    Ewoms::data::Connection val1 = getConnection();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Connection)
}

BOOST_AUTO_TEST_CASE(dataCurrentControl)
{
    Ewoms::data::CurrentControl val1 = getCurrentControl();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::CurrentControl)
}

BOOST_AUTO_TEST_CASE(dataSegment)
{
    Ewoms::data::Segment val1 = getSegment();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Segment)
}

BOOST_AUTO_TEST_CASE(dataWell)
{
    Ewoms::data::Well val1 = getWell();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Well)
}

BOOST_AUTO_TEST_CASE(WellRates)
{
    Ewoms::data::WellRates val1;
    val1.insert({"test_well", getWell()});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::WellRates)
}

BOOST_AUTO_TEST_CASE(dataGroupConstraints)
{
    const auto val1 = getGroupConstraints();
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::GroupConstraints)
}

BOOST_AUTO_TEST_CASE(dataGroupGuideRates)
{
    const auto val1 = getGroupData().guideRates;
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::GroupGuideRates)
}

BOOST_AUTO_TEST_CASE(dataGroupData)
{
    const auto val1 = getGroupData();
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::GroupData)
}

BOOST_AUTO_TEST_CASE(CellData)
{
    Ewoms::data::CellData val1;
    val1.dim = Ewoms::UnitSystem::measure::length;
    val1.data = {1.0, 2.0, 3.0};
    val1.target = Ewoms::data::TargetType::RESTART_SOLUTION;
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::cellData)
}

BOOST_AUTO_TEST_CASE(RestartKey)
{
    Ewoms::RestartKey val1("key", Ewoms::UnitSystem::measure::length, true);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RestartKey)
}

BOOST_AUTO_TEST_CASE(RestartValue)
{
    auto wells1 = Ewoms::data::WellRates {{
        { "test_well", getWell() },
    }};
    auto groups1 = Ewoms::data::GroupValues {{
        { "test_group1", getGroupData() },
    }};

    const auto val1 = Ewoms::RestartValue {
        getSolution(), std::move(wells1), std::move(groups1)
    };
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(RestartValue)
}

#define TEST_FOR_TYPE(TYPE) \
BOOST_AUTO_TEST_CASE(TYPE) \
{ \
    auto val1 = Ewoms::TYPE::serializeObject(); \
    auto val2 = PackUnpack2(val1); \
    DO_CHECKS(TYPE) \
}

#define TEST_FOR_TYPE2(TYPE1, TYPE2) \
BOOST_AUTO_TEST_CASE(TYPE2) \
{ \
    auto val1 = Ewoms::TYPE1::TYPE2::serializeObject(); \
    auto val2 = PackUnpack2(val1); \
    DO_CHECKS(TYPE1::TYPE2) \
}

TEST_FOR_TYPE(Actdims)
TEST_FOR_TYPE(Aqudims)
TEST_FOR_TYPE(Aquancon)
TEST_FOR_TYPE(AquiferConfig)
TEST_FOR_TYPE(AquiferCT)
TEST_FOR_TYPE(Aquifetp)
TEST_FOR_TYPE2(Action, Actions)
TEST_FOR_TYPE2(Action, ActionX)
TEST_FOR_TYPE2(Action, AST)
TEST_FOR_TYPE2(Action, ASTNode)
TEST_FOR_TYPE(BCConfig)
TEST_FOR_TYPE(BrineDensityTable)
TEST_FOR_TYPE(ColumnSchema)
TEST_FOR_TYPE(Connection)
TEST_FOR_TYPE(Deck)
TEST_FOR_TYPE(DeckItem)
TEST_FOR_TYPE(DeckKeyword)
TEST_FOR_TYPE(DeckRecord)
TEST_FOR_TYPE(DensityTable)
TEST_FOR_TYPE(DenT)
TEST_FOR_TYPE(Dimension)
TEST_FOR_TYPE(EclHysterConfig)
TEST_FOR_TYPE(EclipseConfig)
TEST_FOR_TYPE(EDITNNC)
TEST_FOR_TYPE(EndpointScaling)
TEST_FOR_TYPE(Eqldims)
TEST_FOR_TYPE(Equil)
TEST_FOR_TYPE(TLMixpar)
TEST_FOR_TYPE(Events)
TEST_FOR_TYPE(Fault)
TEST_FOR_TYPE(FaultCollection)
TEST_FOR_TYPE(FaultFace)
TEST_FOR_TYPE(FoamConfig)
TEST_FOR_TYPE(FoamData)
TEST_FOR_TYPE(GConSale)
TEST_FOR_TYPE(GConSump)
TEST_FOR_TYPE(GridDims)
TEST_FOR_TYPE(Group)
TEST_FOR_TYPE2(Group, GroupInjectionProperties)
TEST_FOR_TYPE2(Group, GroupProductionProperties)
TEST_FOR_TYPE(GuideRateConfig)
TEST_FOR_TYPE(GuideRateModel)
TEST_FOR_TYPE(InitConfig)
TEST_FOR_TYPE(IOConfig)
TEST_FOR_TYPE(JFunc)
TEST_FOR_TYPE(KeywordLocation)
TEST_FOR_TYPE(MessageLimits)
TEST_FOR_TYPE(MLimits)
TEST_FOR_TYPE(MULTREGTScanner)
TEST_FOR_TYPE(NNC)
TEST_FOR_TYPE2(Network, Node)
TEST_FOR_TYPE(OilVaporizationProperties)
TEST_FOR_TYPE(Phases)
TEST_FOR_TYPE(PlymwinjTable)
TEST_FOR_TYPE(PlyshlogTable)
TEST_FOR_TYPE(PvcdoTable)
TEST_FOR_TYPE(PvtgTable)
TEST_FOR_TYPE(PvtoTable)
TEST_FOR_TYPE(PvtwsaltTable)
TEST_FOR_TYPE(PvtwTable)
TEST_FOR_TYPE(Regdims)
TEST_FOR_TYPE(RestartConfig)
TEST_FOR_TYPE(RestartSchedule)
TEST_FOR_TYPE(RFTConfig)
TEST_FOR_TYPE(RockConfig)
TEST_FOR_TYPE(RockTable)
TEST_FOR_TYPE(RocktabTable)
TEST_FOR_TYPE(Rock2dtrTable)
TEST_FOR_TYPE(Rock2dTable)
TEST_FOR_TYPE(Runspec)
TEST_FOR_TYPE(Schedule)
TEST_FOR_TYPE(Segment)
TEST_FOR_TYPE(SimpleTable)
TEST_FOR_TYPE(SimulationConfig)
TEST_FOR_TYPE(SkprpolyTable)
TEST_FOR_TYPE(SkprwatTable)
TEST_FOR_TYPE(SICD)
TEST_FOR_TYPE(SolventDensityTable)
TEST_FOR_TYPE(SummaryConfig)
TEST_FOR_TYPE(SummaryConfigNode)
TEST_FOR_TYPE(Tabdims)
TEST_FOR_TYPE(TableColumn)
TEST_FOR_TYPE(TableContainer)
TEST_FOR_TYPE(TableManager)
TEST_FOR_TYPE(TableSchema)
TEST_FOR_TYPE(ThresholdPressure)
TEST_FOR_TYPE(TimeMap)
TEST_FOR_TYPE(TracerConfig)
TEST_FOR_TYPE(TransMult)
TEST_FOR_TYPE(Tuning)
TEST_FOR_TYPE(UDAValue)
TEST_FOR_TYPE(UDQAssign)
TEST_FOR_TYPE(UDQActive)
TEST_FOR_TYPE(UDQASTNode)
TEST_FOR_TYPE(UDQConfig)
TEST_FOR_TYPE(UDQDefine)
TEST_FOR_TYPE(UDQIndex)
TEST_FOR_TYPE(UDQParams)
TEST_FOR_TYPE(UnitSystem)
TEST_FOR_TYPE(Valve)
TEST_FOR_TYPE(VFPInjTable)
TEST_FOR_TYPE(VFPProdTable)
TEST_FOR_TYPE(ViscrefTable)
TEST_FOR_TYPE(WatdentTable)
TEST_FOR_TYPE(Well)
TEST_FOR_TYPE(Welldims)
TEST_FOR_TYPE(WellBrineProperties)
TEST_FOR_TYPE(WellConnections)
TEST_FOR_TYPE(WellEconProductionLimits)
TEST_FOR_TYPE(WellFoamProperties)
TEST_FOR_TYPE2(Well, WellGuideRate)
TEST_FOR_TYPE2(Well, WellInjectionProperties)
TEST_FOR_TYPE(WellPolymerProperties)
TEST_FOR_TYPE2(Well, WellProductionProperties)
TEST_FOR_TYPE(WellTracerProperties)
TEST_FOR_TYPE(WellSegmentDims)
TEST_FOR_TYPE(WellSegments)
TEST_FOR_TYPE(WellTestConfig)
TEST_FOR_TYPE(WellType)
TEST_FOR_TYPE(WListManager)

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
