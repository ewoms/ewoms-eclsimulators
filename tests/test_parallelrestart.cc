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

#include <ewoms/eclio/opmlog/location.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/deck/deckitem.hh>
#include <ewoms/eclio/parser/eclipsestate/aquancon.hh>
#include <ewoms/eclio/parser/eclipsestate/aquiferct.hh>
#include <ewoms/eclio/parser/eclipsestate/aquiferconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/aquifetp.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipseconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/runspec.hh>
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
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/spiralicd.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/valve.hh>
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

#if HAVE_MPI
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
    seg1.pressure = 2.0;
    return seg1;
}

Ewoms::data::CurrentControl getCurrentControl()
{
    Ewoms::data::CurrentControl curr;
    curr.isProducer = true;
    curr.prod = ::Ewoms::Well::ProducerCMode::CRAT;
    return curr;
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
    return well1;
}

Ewoms::ThresholdPressure getThresholdPressure()
{
    return Ewoms::ThresholdPressure(false, true, {{true, 1.0}, {false, 2.0}},
                                  {{{1,2},{false,3.0}},{{2,3},{true,4.0}}});
}

Ewoms::RockConfig getRockConfig()
{
    return Ewoms::RockConfig(true, {{100, 0.25}, {200, 0.30}}, "ROCKNUM", 10, false, Ewoms::RockConfig::Hysteresis::HYSTER);
}

Ewoms::BCConfig getBCConfig()
{
    return Ewoms::BCConfig({{10,11,12,13,14,15,Ewoms::BCType::RATE,Ewoms::FaceDir::XPlus, Ewoms::BCComponent::GAS, 100.0}});
}

Ewoms::TableSchema getTableSchema()
{
    Ewoms::OrderedMap<std::string, Ewoms::ColumnSchema> data;
    data.insert({"test1", Ewoms::ColumnSchema("test1", Ewoms::Table::INCREASING,
                                                     Ewoms::Table::DEFAULT_LINEAR)});
    data.insert({"test2", Ewoms::ColumnSchema("test2", Ewoms::Table::INCREASING, 1.0)});
    return Ewoms::TableSchema(data);
}

Ewoms::TableColumn getTableColumn()
{
    return Ewoms::TableColumn(Ewoms::ColumnSchema("test1", Ewoms::Table::INCREASING,
                                              Ewoms::Table::DEFAULT_LINEAR),
                            "test2", {1.0, 2.0}, {false, true}, 2);
}

Ewoms::DenT getDenT()
{
    std::vector<Ewoms::DenT::entry> records;
    records.emplace_back(1,2,3);
    records.emplace_back(4,5,6);
    return Ewoms::DenT(records);
}

Ewoms::SimpleTable getSimpleTable()
{
    Ewoms::OrderedMap<std::string, Ewoms::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    return Ewoms::SimpleTable(getTableSchema(), data, true);
}

Ewoms::EquilRecord getEquilRecord()
{
    return Ewoms::EquilRecord(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, true, false, 1);
}

Ewoms::FoamData getFoamData()
{
    return Ewoms::FoamData(1.0, 2.0, 3.0, true, 4.0);
}

Ewoms::FoamConfig getFoamConfig()
{
    return Ewoms::FoamConfig({getFoamData(), getFoamData()},
                           Ewoms::Phase::GAS,
                           Ewoms::FoamConfig::MobilityModel::TAB);
}

Ewoms::TimeMap getTimeMap()
{
    return Ewoms::TimeMap({123});
}

Ewoms::RestartConfig getRestartConfig()
{
    Ewoms::DynamicState<Ewoms::RestartSchedule> rsched({Ewoms::RestartSchedule(1, 2, 3)}, 2);
    Ewoms::DynamicState<std::map<std::string,int>> rkw({{{"test",3}}}, 3);
    Ewoms::IOConfig io(true, false, true, false, false, true, "test1", true,
                     "test2", true, "test3", false);
    return Ewoms::RestartConfig(getTimeMap(), 1, true, rsched, rkw, {false, true});
}

Ewoms::PvtgTable getPvtgTable()
{
    return Ewoms::PvtgTable(Ewoms::ColumnSchema("test1", Ewoms::Table::INCREASING,
                                            Ewoms::Table::DEFAULT_LINEAR),
                          getTableColumn(),
                          getTableSchema(),
                          getTableSchema(),
                          {getSimpleTable()},
                          getSimpleTable());
}

Ewoms::PvtoTable getPvtoTable()
{
    return Ewoms::PvtoTable(Ewoms::ColumnSchema("test1", Ewoms::Table::INCREASING,
                                            Ewoms::Table::DEFAULT_LINEAR),
                          getTableColumn(),
                          getTableSchema(),
                          getTableSchema(),
                          {getSimpleTable()},
                          getSimpleTable());
}

Ewoms::TableContainer getTableContainer()
{
    Ewoms::OrderedMap<std::string, Ewoms::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    Ewoms::SimpleTable tab1(getTableSchema(), data, true);
    Ewoms::TableContainer result(2);
    result.addTable(0, std::make_shared<Ewoms::SimpleTable>(tab1));
    result.addTable(1, std::make_shared<Ewoms::SimpleTable>(tab1));
    return result;
}

Ewoms::Well getFullWell()
{
    Ewoms::UnitSystem unitSystem;
    return Ewoms::Well("test1", "test2", 1, 2, 3, 4, 5.0,
                     Ewoms::WellType(Ewoms::Phase::WATER),
                     unitSystem, 6.0, Ewoms::Well::Status::SHUT,
                     7.0, true, false,
                     Ewoms::Well::WellGuideRate{true, 1.0, Ewoms::Well::GuideRateTarget::COMB, 2.0},
                     8.0, 9.0, false,
                     std::make_shared<Ewoms::WellEconProductionLimits>(),
                     std::make_shared<Ewoms::WellFoamProperties>(),
                     std::make_shared<Ewoms::WellPolymerProperties>(),
                     std::make_shared<Ewoms::WellBrineProperties>(),
                     std::make_shared<Ewoms::WellTracerProperties>(),
                     std::make_shared<Ewoms::WellConnections>(),
                     std::make_shared<Ewoms::Well::WellProductionProperties>(),
                     std::make_shared<Ewoms::Well::WellInjectionProperties>(),
                     std::make_shared<Ewoms::WellSegments>());
}

Ewoms::VFPInjTable getVFPInjTable()
{
    Ewoms::VFPInjTable::array_type table;
    table.resize(3*2);
    std::iota(table.begin(), table.end(), 1.0);

    return Ewoms::VFPInjTable(1, 2.0, Ewoms::VFPInjTable::FLO_WAT, {1.0, 2.0},
                            {3.0, 4.0, 5.0}, table);
}

Ewoms::VFPProdTable getVFPProdTable()
{
    Ewoms::VFPProdTable::array_type table;
    table.resize(1*2*3*4*5);
    std::iota(table.begin(), table.end(), 1.0);

    return Ewoms::VFPProdTable(1, 2.0, Ewoms::VFPProdTable::FLO_OIL,
                             Ewoms::VFPProdTable::WFR_WOR,
                             Ewoms::VFPProdTable::GFR_GLR,
                             Ewoms::VFPProdTable::ALQ_TGLR,
                             {1.0, 2.0, 3.0, 4.0, 5.0},
                             {1.0},
                             {1.0, 2.0},
                             {1.0, 2.0, 3.0},
                             {1.0, 2.0, 3.0, 4.0}, table);
}

Ewoms::UDQConfig getUDQConfig()
{
    Ewoms::UDQParams params(true, 1, 2.0, 3.0, 4.0);
    std::shared_ptr<Ewoms::UDQASTNode> n0;
    Ewoms::UDQASTNode n1(Ewoms::UDQVarType::NONE,
                       Ewoms::UDQTokenType::error,
                       "test", 1.0, {"test1", "test2"}, n0, n0);
    Ewoms::UDQDefine def("test", std::make_shared<Ewoms::UDQASTNode>(n1),
                       Ewoms::UDQVarType::NONE, "test2");
    Ewoms::UDQAssign ass("test", Ewoms::UDQVarType::NONE,
                       {Ewoms::UDQAssign::AssignRecord{{"test1"}, 1.0},
                        Ewoms::UDQAssign::AssignRecord{{"test2"}, 2.0}});
    Ewoms::OrderedMap<std::string, Ewoms::UDQIndex> omap;
    omap.insert({"test8", Ewoms::UDQIndex(1, 2, Ewoms::UDQAction::ASSIGN,
                                        Ewoms::UDQVarType::WELL_VAR)});
    omap.insert({"test9", Ewoms::UDQIndex(3, 4, Ewoms::UDQAction::ASSIGN,
                                        Ewoms::UDQVarType::WELL_VAR)});
    return Ewoms::UDQConfig(params,
                          Ewoms::UDQFunctionTable(params),
                          {{"test1", def}, {"test2", def}},
                          {{"test3", ass}, {"test4", ass}},
                          {{"test5", "test6"}, {"test7", "test8"}},
                          omap,
                          {{Ewoms::UDQVarType::SCALAR, 5}, {Ewoms::UDQVarType::WELL_VAR, 6}});
}

Ewoms::GuideRateModel getGuideRateModel()
{
    return Ewoms::GuideRateModel(1.0, Ewoms::GuideRateModel::Target::WAT,
                               {2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
                               true, 8.0, false, false,
                               {Ewoms::UDAValue(9.0),
                               Ewoms::UDAValue(10.0),
                               Ewoms::UDAValue(11.0)});
}

Ewoms::GuideRateConfig::GroupTarget getGuideRateConfigGroup()
{
    return Ewoms::GuideRateConfig::GroupTarget{1.0, Ewoms::Group::GuideRateTarget::COMB};
}

Ewoms::GuideRateConfig::WellTarget getGuideRateConfigWell()
{
    return Ewoms::GuideRateConfig::WellTarget{1.0, Ewoms::Well::GuideRateTarget::COMB, 2.0};
}

Ewoms::DeckRecord getDeckRecord()
{
    Ewoms::DeckItem item1({1.0}, {2}, {"test3"}, {Ewoms::UDAValue(4)},
                       Ewoms::type_tag::string, "test5",
                       {Ewoms::value::status::deck_value},
                       true,
                       {Ewoms::Dimension(7.0, 8.0)},
                       {Ewoms::Dimension(10.0, 11.0)});

    Ewoms::DeckItem item2({1.0}, {2}, {"test3"}, {Ewoms::UDAValue(4)},
                       Ewoms::type_tag::string, "test6",
                       {Ewoms::value::status::deck_value},
                       true,
                       {Ewoms::Dimension(7.0, 8.0)},
                       {Ewoms::Dimension(10.0, 11.0)});

    return Ewoms::DeckRecord({item1, item2});
}

Ewoms::Tuning getTuning()
{
    return Ewoms::Tuning{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, true,
                       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
                       20.0, 21.0, false, 22.0, 3, 4, 5, 6, 7, 8, 9, 23.0, 24.0,
                       25.0, 26.0, true};
}

Ewoms::Action::Condition getCondition()
{
    Ewoms::Action::Quantity q;
    q.quantity = "test1";
    q.args = {"test2", "test3"};
    Ewoms::Action::Condition val1;
    val1.lhs = val1.rhs = q;
    val1.logic = Ewoms::Action::Condition::Logical::OR;
    val1.cmp = Ewoms::Action::Condition::Comparator::LESS;
    val1.cmp_string = "test";
    return val1;
}

Ewoms::Action::ActionX getActionX()
{
    std::shared_ptr<Ewoms::Action::ASTNode> node;
    node.reset(new Ewoms::Action::ASTNode(number, FuncType::field,
                                        "test1", {"test2"}, 1.0, {}));
    Ewoms::Action::AST ast(node);
    return Ewoms::Action::ActionX("test", 1, 2.0, 3,
                                {Ewoms::DeckKeyword("test", {"test",1},
                                                  {getDeckRecord(), getDeckRecord()},
                                                  true, false)},
                                ast, {getCondition()}, 4, 5);
}

Ewoms::AquiferCT getAquiferCT()
{
    Ewoms::AquiferCT::AQUCT_data data;
    data.aquiferID = 1;
    data.inftableID = 2;
    data.pvttableID = 3;
    data.phi_aq = 100;
    data.d0 = 1;
    data.C_t = 10;
    data.r_o = 1.5;
    data.k_a = 100;
    data.c1 = 0.78;
    data.h = 1;
    data.c2 = 45;
    data.p0 = std::make_pair(true, 98);
    data.td = {1,2,3};
    data.pi = {4,5,6};
    data.cell_id = {0,10,100};

    return Ewoms::AquiferCT( { data } );
}

Ewoms::Aquifetp getAquifetp()
{
    Ewoms::Aquifetp::AQUFETP_data data;

    data.aquiferID = 1;
    data.pvttableID = 3;
    data.C_t = 10;
    data.p0 = std::make_pair(true, 98);
    data.V0 = 0;
    data.d0 = 0;

    return Ewoms::Aquifetp( { data } );
}

Ewoms::Aquancon getAquancon()
{
    Ewoms::Aquancon::AquancCell cell(1, 100, std::make_pair(false, 0), 100, Ewoms::FaceDir::XPlus);
    return Ewoms::Aquancon( std::unordered_map<int, std::vector<Ewoms::Aquancon::AquancCell>>{{1, {cell}}});
}

#endif

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
#if HAVE_MPI
    Ewoms::data::Solution val1 = getSolution();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Solution)
#endif
}

BOOST_AUTO_TEST_CASE(Rates)
{
#if HAVE_MPI
    Ewoms::data::Rates val1 = getRates();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Rates)
#endif
}

BOOST_AUTO_TEST_CASE(dataConnection)
{
#if HAVE_MPI
    Ewoms::data::Connection val1 = getConnection();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Connection)
#endif
}

BOOST_AUTO_TEST_CASE(dataCurrentControl)
{
#if HAVE_MPI
    Ewoms::data::CurrentControl cur1 = getCurrentControl();
    auto cur2 = PackUnpack(cur1);
    BOOST_CHECK(std::get<1>(cur2) == std::get<2>(cur2));
    BOOST_CHECK(cur1 == std::get<0>(cur2));
#endif
}

BOOST_AUTO_TEST_CASE(dataSegment)
{
#if HAVE_MPI
    Ewoms::data::Segment val1 = getSegment();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Segment)
#endif
}

BOOST_AUTO_TEST_CASE(dataWell)
{
#if HAVE_MPI
    Ewoms::data::Well val1 = getWell();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Well)
#endif
}

BOOST_AUTO_TEST_CASE(WellRates)
{
#if HAVE_MPI
    Ewoms::data::WellRates val1;
    val1.insert({"test_well", getWell()});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::WellRates)
#endif
}

BOOST_AUTO_TEST_CASE(CellData)
{
#if HAVE_MPI
    Ewoms::data::CellData val1;
    val1.dim = Ewoms::UnitSystem::measure::length;
    val1.data = {1.0, 2.0, 3.0};
    val1.target = Ewoms::data::TargetType::RESTART_SOLUTION;
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::cellData)
#endif
}

BOOST_AUTO_TEST_CASE(RestartKey)
{
#if HAVE_MPI
    Ewoms::RestartKey val1("key", Ewoms::UnitSystem::measure::length, true);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RestartKey)
#endif
}

BOOST_AUTO_TEST_CASE(RestartValue)
{
#if HAVE_MPI
    Ewoms::data::WellRates wells1;
    wells1.insert({"test_well", getWell()});
    Ewoms::RestartValue val1(getSolution(), wells1);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RestartValue)
#endif
}

BOOST_AUTO_TEST_CASE(ThresholdPressure)
{
#if HAVE_MPI
    Ewoms::ThresholdPressure val1 = getThresholdPressure();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(ThresholdPressure)
#endif
}

BOOST_AUTO_TEST_CASE(RockConfig)
{
#if HAVE_MPI
    Ewoms::RockConfig val1 = getRockConfig();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(RockConfig)
#endif
}

BOOST_AUTO_TEST_CASE(EDITNNC)
{
#if HAVE_MPI
    Ewoms::EDITNNC val1({{1,2,1.0},{2,3,2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(EDITNNC)
#endif
}

BOOST_AUTO_TEST_CASE(NNC)
{
#if HAVE_MPI
    Ewoms::NNC val1({{1,2,1.0},{2,3,2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(NNC)
#endif
}

BOOST_AUTO_TEST_CASE(Rock2dTable)
{
#if HAVE_MPI
    Ewoms::Rock2dTable val1({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Rock2dTable)
#endif
}

BOOST_AUTO_TEST_CASE(Rock2dtrTable)
{
#if HAVE_MPI
    Ewoms::Rock2dtrTable val1({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Rock2dtrTable)
#endif
}

BOOST_AUTO_TEST_CASE(ColumnSchema)
{
#if HAVE_MPI
    Ewoms::ColumnSchema val1("test1", Ewoms::Table::INCREASING,
                           Ewoms::Table::DEFAULT_LINEAR);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(ColumnSchema)
    val1 = Ewoms::ColumnSchema("test2", Ewoms::Table::DECREASING, 1.0);
    val2 = PackUnpack2(val1);
    DO_CHECKS(ColumnSchema)
#endif
}

BOOST_AUTO_TEST_CASE(TableSchema)
{
#if HAVE_MPI
    Ewoms::TableSchema val1 = getTableSchema();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(TableSchema)
#endif
}

BOOST_AUTO_TEST_CASE(TableColumn)
{
#if HAVE_MPI
    Ewoms::TableColumn val1 = getTableColumn();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(TableColumn)
#endif
}

BOOST_AUTO_TEST_CASE(SimpleTable)
{
#if HAVE_MPI
    Ewoms::SimpleTable val1 = getSimpleTable();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SimpleTable)
#endif
}

BOOST_AUTO_TEST_CASE(TableContainer)
{
#if HAVE_MPI
    Ewoms::OrderedMap<std::string, Ewoms::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    Ewoms::SimpleTable tab1(getTableSchema(), data, true);
    Ewoms::TableContainer val1(2);
    val1.addTable(0, std::make_shared<Ewoms::SimpleTable>(tab1));
    val1.addTable(1, std::make_shared<Ewoms::SimpleTable>(tab1));
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(TableContainer)
#endif
}

BOOST_AUTO_TEST_CASE(EquilRecord)
{
#if HAVE_MPI
    Ewoms::EquilRecord val1 = getEquilRecord();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(EquilRecord)
#endif
}

BOOST_AUTO_TEST_CASE(Equil)
{
#if HAVE_MPI
    Ewoms::Equil val1({getEquilRecord(), getEquilRecord()});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Equil)
#endif
}

BOOST_AUTO_TEST_CASE(FoamData)
{
#if HAVE_MPI
    Ewoms::FoamData val1 = getFoamData();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(FoamData)
#endif
}

BOOST_AUTO_TEST_CASE(FoamConfig)
{
#if HAVE_MPI
    Ewoms::FoamConfig val1 = getFoamConfig();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(FoamConfig)
#endif
}

BOOST_AUTO_TEST_CASE(InitConfig)
{
#if HAVE_MPI
    Ewoms::InitConfig val1(Ewoms::Equil({getEquilRecord(), getEquilRecord()}),
                         getFoamConfig(),
                         true, true, true, 20, "test1");
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(InitConfig)
#endif
}

BOOST_AUTO_TEST_CASE(SimulationConfig)
{
#if HAVE_MPI
    Ewoms::SimulationConfig val1(getThresholdPressure(), getBCConfig(), getRockConfig(), false, true, false, true);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SimulationConfig)
#endif
}

BOOST_AUTO_TEST_CASE(BCConfig)
{
#if HAVE_MPI
    Ewoms::BCConfig val1({{10,11,12,13,14,15,Ewoms::BCType::RATE, Ewoms::FaceDir::XPlus, Ewoms::BCComponent::GAS, 100}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(BCConfig)
#endif
}

BOOST_AUTO_TEST_CASE(RestartSchedule)
{
#if HAVE_MPI
    Ewoms::RestartSchedule val1(1, 2, 3);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(RestartSchedule)
#endif
}

BOOST_AUTO_TEST_CASE(TimeMap)
{
#if HAVE_MPI
    Ewoms::TimeMap val1 = getTimeMap();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(TimeMap)
#endif
}

BOOST_AUTO_TEST_CASE(RestartConfig)
{
#if HAVE_MPI
    Ewoms::DynamicState<Ewoms::RestartSchedule> rsched({Ewoms::RestartSchedule(1, 2, 3)}, 2);
    Ewoms::DynamicState<std::map<std::string,int>> rkw({{{"test",3}}}, 3);
    Ewoms::IOConfig io(true, false, true, false, false, true, "test1", true,
                     "test2", true, "test3", false);
    Ewoms::RestartConfig val1(getTimeMap(), 1, true, rsched, rkw, {false, true});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(RestartConfig)
#endif
}

BOOST_AUTO_TEST_CASE(IOConfig)
{
#if HAVE_MPI
    Ewoms::IOConfig val1(true, false, true, false, false, true, "test1", true,
                       "test2", true, "test3", false);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(IOConfig)
#endif
}

BOOST_AUTO_TEST_CASE(Phases)
{
#if HAVE_MPI
    Ewoms::Phases val1(true, true, true, false, true, false, true, false);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Phases)
#endif
}

BOOST_AUTO_TEST_CASE(Tabdims)
{
#if HAVE_MPI
    Ewoms::Tabdims val1(1,2,3,4,5,6);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Tabdims)
#endif
}

BOOST_AUTO_TEST_CASE(EndpointScaling)
{
#if HAVE_MPI
    Ewoms::EndpointScaling val1(std::bitset<4>(13));
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(EndpointScaling)
#endif
}

BOOST_AUTO_TEST_CASE(Welldims)
{
#if HAVE_MPI
    Ewoms::Welldims val1(1,2,3,4);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Welldims)
#endif
}

BOOST_AUTO_TEST_CASE(WellSegmentDims)
{
#if HAVE_MPI
    Ewoms::WellSegmentDims val1(1,2,3);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellSegmentDims)
#endif
}

BOOST_AUTO_TEST_CASE(UDQParams)
{
#if HAVE_MPI
    Ewoms::UDQParams val1(true, 1, 2.0, 3.0, 4.0);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDQParams)
#endif
}

BOOST_AUTO_TEST_CASE(EclHysterConfig)
{
#if HAVE_MPI
    Ewoms::EclHysterConfig val1(true, 1, 2);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(EclHysterConfig)
#endif
}

BOOST_AUTO_TEST_CASE(Actdims)
{
#if HAVE_MPI
    Ewoms::Actdims val1(1,2,3,4);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Actdims)
#endif
}

BOOST_AUTO_TEST_CASE(Runspec)
{
#if HAVE_MPI
    Ewoms::Runspec val1(Ewoms::Phases(true, true, true, false, true, false, true, false),
                      Ewoms::Tabdims(1,2,3,4,5,6),
                      Ewoms::EndpointScaling(std::bitset<4>(13)),
                      Ewoms::Welldims(1,2,3,4),
                      Ewoms::WellSegmentDims(1,2,3),
                      Ewoms::UDQParams(true, 1, 2.0, 3.0, 4.0),
                      Ewoms::EclHysterConfig(true, 1, 2),
                      Ewoms::Actdims(1,2,3,4),
                      Ewoms::SatFuncControls(5.0e-7, Ewoms::SatFuncControls::ThreePhaseOilKrModel::Stone2));

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Runspec)
#endif
}

BOOST_AUTO_TEST_CASE(PvtgTable)
{
#if HAVE_MPI
    Ewoms::PvtgTable val1 = getPvtgTable();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(PvtgTable)
#endif
}

BOOST_AUTO_TEST_CASE(PvtoTable)
{
#if HAVE_MPI
    Ewoms::PvtoTable val1 = getPvtoTable();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(PvtoTable)
#endif
}

BOOST_AUTO_TEST_CASE(JFunc)
{
#if HAVE_MPI
    Ewoms::JFunc val1(Ewoms::JFunc::Flag::BOTH, 1.0, 2.0,
                    3.0, 4.0, Ewoms::JFunc::Direction::XY);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(JFunc)
#endif
}

BOOST_AUTO_TEST_CASE(PVTWRecord)
{
#if HAVE_MPI
    Ewoms::PVTWRecord val1{1.0, 2.0, 3.0, 4.0, 5.0};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(PVTWRecord)
#endif
}

BOOST_AUTO_TEST_CASE(PvtwTable)
{
#if HAVE_MPI
    Ewoms::PvtwTable val1({Ewoms::PVTWRecord{1.0, 2.0, 3.0, 4.0, 5.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(PvtwTable)
#endif
}

BOOST_AUTO_TEST_CASE(PVCDORecord)
{
#if HAVE_MPI
    Ewoms::PVTWRecord val1{1.0, 2.0, 3.0, 4.0, 5.0};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(PVTWRecord)
#endif
}

BOOST_AUTO_TEST_CASE(PvcdoTable)
{
#if HAVE_MPI
    Ewoms::PvcdoTable val1({Ewoms::PVCDORecord{1.0, 2.0, 3.0, 4.0, 5.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(PvcdoTable)
#endif
}

BOOST_AUTO_TEST_CASE(DENSITYRecord)
{
#if HAVE_MPI
    Ewoms::DENSITYRecord val1{1.0, 2.0, 3.0};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(DENSITYRecord)
#endif
}

BOOST_AUTO_TEST_CASE(DensityTable)
{
#if HAVE_MPI
    Ewoms::DensityTable val1({Ewoms::DENSITYRecord{1.0, 2.0, 3.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(DensityTable)
#endif
}

BOOST_AUTO_TEST_CASE(VISCREFRecord)
{
#if HAVE_MPI
    Ewoms::VISCREFRecord val1{1.0, 2.0};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(VISCREFRecord)
#endif
}

BOOST_AUTO_TEST_CASE(ViscrefTable)
{
#if HAVE_MPI
    Ewoms::ViscrefTable val1({Ewoms::VISCREFRecord{1.0, 2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(ViscrefTable)
#endif
}

BOOST_AUTO_TEST_CASE(WATDENTRecord)
{
#if HAVE_MPI
    Ewoms::WATDENTRecord val1{1.0, 2.0, 3.0};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WATDENTRecord)
#endif
}

BOOST_AUTO_TEST_CASE(WatdentTable)
{
#if HAVE_MPI
    Ewoms::WatdentTable val1({Ewoms::WATDENTRecord{1.0, 2.0, 3.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WatdentTable)
#endif
}

BOOST_AUTO_TEST_CASE(PlymwinjTable)
{
#if HAVE_MPI
    Ewoms::PlymwinjTable val1({1.0}, {2.0}, 1, {{1.0}, {2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(PlymwinjTable)
#endif
}

BOOST_AUTO_TEST_CASE(SkprpolyTable)
{
#if HAVE_MPI
    Ewoms::SkprpolyTable val1({1.0}, {2.0}, 1, {{1.0}, {2.0}}, 3.0);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SkprpolyTable)
#endif
}

BOOST_AUTO_TEST_CASE(SkprwatTable)
{
#if HAVE_MPI
    Ewoms::SkprwatTable val1({1.0}, {2.0}, 1, {{1.0}, {2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SkprwatTable)
#endif
}

BOOST_AUTO_TEST_CASE(Regdims)
{
#if HAVE_MPI
    Ewoms::Regdims val1(1,2,3,4,5);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Regdims)
#endif
}

BOOST_AUTO_TEST_CASE(Eqldims)
{
#if HAVE_MPI
    Ewoms::Eqldims val1(1,2,3,4,5);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Eqldims)
#endif
}

BOOST_AUTO_TEST_CASE(Aqudims)
{
#if HAVE_MPI
    Ewoms::Aqudims val1(1,2,3,4,5,6,7,8);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Aqudims)
#endif
}

BOOST_AUTO_TEST_CASE(ROCKRecord)
{
#if HAVE_MPI
    Ewoms::ROCKRecord val1{1.0,2.0};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(ROCKRecord)
#endif
}

BOOST_AUTO_TEST_CASE(RockTable)
{
#if HAVE_MPI
    Ewoms::RockTable val1({Ewoms::ROCKRecord{1.0,2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(RockTable)
#endif
}

BOOST_AUTO_TEST_CASE(TableManager)
{
#if HAVE_MPI
    auto jfunc = std::make_shared<Ewoms::JFunc>(Ewoms::JFunc::Flag::BOTH,
                                              1.0, 2.0, 3.0, 4.0,
                                              Ewoms::JFunc::Direction::XY);
    Ewoms::TableManager val1({{"test", getTableContainer()}},
                           {getPvtgTable()},
                           {getPvtoTable()},
                           {Ewoms::Rock2dTable({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0})},
                           {Ewoms::Rock2dtrTable({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0})},
                           Ewoms::PvtwTable({Ewoms::PVTWRecord{1.0, 2.0, 3.0, 4.0, 5.0}}),
                           Ewoms::PvcdoTable({Ewoms::PVCDORecord{1.0, 2.0, 3.0, 4.0, 5.0}}),
                           Ewoms::DensityTable({Ewoms::DENSITYRecord{1.0, 2.0, 3.0}}),
                           Ewoms::PlyvmhTable({Ewoms::PlyvmhRecord{1.0, 2.0, 3.0, 4.0}}),
                           Ewoms::RockTable({Ewoms::ROCKRecord{1.0,2.0}}),
                           Ewoms::PlmixparTable({Ewoms::PlmixparRecord{1.0}}),
                           Ewoms::ShrateTable({Ewoms::ShrateRecord{1.0}}),
                           Ewoms::Stone1exTable({Ewoms::Stone1exRecord{1.0}}),
                           Ewoms::TlmixparTable({Ewoms::TlmixparRecord{1.0, 2.0}}),
                           Ewoms::ViscrefTable({Ewoms::VISCREFRecord{1.0, 2.0}}),
                           Ewoms::WatdentTable({Ewoms::WATDENTRecord{1.0, 2.0, 3.0}}),
                           {{1.0, 2.0, {1.0, 2.0, 3.0}}},
                           {{{1.0, 2.0, 3.0}}},
                           {{{4.0, 5.0, 6.0}}},
                           {{1, Ewoms::PlymwinjTable({1.0}, {2.0}, 1, {{1.0}, {2.0}})}},
                           {{2, Ewoms::SkprwatTable({1.0}, {2.0}, 1, {{1.0}, {2.0}})}},
                           {{3, Ewoms::SkprpolyTable({1.0}, {2.0}, 1, {{1.0}, {2.0}}, 3.0)}},
                           Ewoms::Tabdims(1,2,3,4,5,6),
                           Ewoms::Regdims(1,2,3,4,5),
                           Ewoms::Eqldims(1,2,3,4,5),
                           Ewoms::Aqudims(1,2,3,4,5,6,7,8),
                           true,
                           true,
                           true,
                           true,
                           jfunc,
                           getDenT(),
                           getDenT(),
                           getDenT(),
                           {7.0, 8.0},
                           77,
                           1.0);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(TableManager)
#endif
}

BOOST_AUTO_TEST_CASE(OilVaporizationProperties)
{
#ifdef HAVE_MPI
    using VapType = Ewoms::OilVaporizationProperties::OilVaporization;
    Ewoms::OilVaporizationProperties val1(VapType::VAPPARS,
                                        1.0, 2.0, {5.0, 6.0},
                                        {false, true}, {7.0, 8.0});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(OilVaporizationProperties)
    val1 = Ewoms::OilVaporizationProperties(VapType::DRDT,
                                          1.0, 2.0, {5.0, 6.0},
                                          {false, true}, {7.0, 8.0});
    val2 = PackUnpack2(val1);
    DO_CHECKS(OilVaporizationProperties)
#endif
}

BOOST_AUTO_TEST_CASE(Events)
{
#ifdef HAVE_MPI
    Ewoms::Events val1(Ewoms::DynamicVector<uint64_t>({1,2,3,4,5}));
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Events)
#endif
}

BOOST_AUTO_TEST_CASE(MLimits)
{
#ifdef HAVE_MPI
    Ewoms::MLimits val1{1,2,3,4,5,6,7,8,9,10,11,12};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(MLimits)
#endif
}

BOOST_AUTO_TEST_CASE(MessageLimits)
{
#ifdef HAVE_MPI
    std::vector<Ewoms::MLimits> limits{Ewoms::MLimits{1,2,3,4,5,6,7,8,9,10,11,12}};
    Ewoms::MessageLimits val1(Ewoms::DynamicState<Ewoms::MLimits>(limits,2));
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(MessageLimits)
#endif
}

BOOST_AUTO_TEST_CASE(VFPInjTable)
{
#ifdef HAVE_MPI
    Ewoms::VFPInjTable val1 = getVFPInjTable();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(VFPInjTable)
#endif
}

BOOST_AUTO_TEST_CASE(VFPProdTable)
{
#ifdef HAVE_MPI
    Ewoms::VFPProdTable val1 = getVFPProdTable();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(VFPProdTable)
#endif
}

BOOST_AUTO_TEST_CASE(WTESTWell)
{
#ifdef HAVE_MPI
    Ewoms::WellTestConfig::WTESTWell val1{"test", Ewoms::WellTestConfig::ECONOMIC,
                                         1.0, 2, 3.0, 4};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellTestConfig::WTESTWell)
#endif
}

BOOST_AUTO_TEST_CASE(WellTestConfig)
{
#ifdef HAVE_MPI
    Ewoms::WellTestConfig::WTESTWell tw{"test", Ewoms::WellTestConfig::ECONOMIC,
                                         1.0, 2, 3.0, 4};
    Ewoms::WellTestConfig val1({tw, tw, tw});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellTestConfig)
#endif
}

BOOST_AUTO_TEST_CASE(WellPolymerProperties)
{
#ifdef HAVE_MPI
    Ewoms::WellPolymerProperties val1{1.0, 2.0, 3, 4, 5};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellPolymerProperties)
#endif
}

BOOST_AUTO_TEST_CASE(WellFoamProperties)
{
#ifdef HAVE_MPI
    Ewoms::WellFoamProperties val1{1.0};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellFoamProperties)
#endif
}

BOOST_AUTO_TEST_CASE(WellTracerProperties)
{
#ifdef HAVE_MPI
    Ewoms::WellTracerProperties val1({{"test", 1.0}, {"test2", 2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellTracerProperties)
#endif
}

BOOST_AUTO_TEST_CASE(UDAValue)
{
#ifdef HAVE_MPI
    Ewoms::UDAValue val1("test");
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDAValue)
    val1 = Ewoms::UDAValue(1.0);
    val2 = PackUnpack2(val1);
    DO_CHECKS(UDAValue)
#endif
}

BOOST_AUTO_TEST_CASE(Connection)
{
#ifdef HAVE_MPI
    Ewoms::Connection val1(Ewoms::Connection::Direction::Y,
                         1.0, Ewoms::Connection::State::SHUT,
                         2, 3, 4.0, 5.0, 6.0, 7.0, 8.0,
                         {9, 10, 11}, 12345, Ewoms::Connection::CTFKind::Defaulted,
                         12, 13.0, 14.0, true,
                         15, 16, 17.0);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Connection)
#endif
}

BOOST_AUTO_TEST_CASE(WellInjectionProperties)
{
#ifdef HAVE_MPI
    Ewoms::Well::WellInjectionProperties val1("test",
                                            Ewoms::UDAValue(1.0),
                                            Ewoms::UDAValue("test"),
                                            Ewoms::UDAValue(2.0),
                                            Ewoms::UDAValue(3.0),
                                            2.0, 3.0,
                                            4.0, 5.0, 6.0,
                                            7,
                                            true,
                                            8,
                                            Ewoms::InjectorType::OIL,
                                            Ewoms::Well::InjectorCMode::BHP);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Well::WellInjectionProperties)
#endif
}

BOOST_AUTO_TEST_CASE(WellEconProductionLimits)
{
#ifdef HAVE_MPI
    Ewoms::WellEconProductionLimits val1(1.0, 2.0, 3.0, 4.0, 5.0,
                                       Ewoms::WellEconProductionLimits::EconWorkover::CONP,
                                       true, "test",
                                       Ewoms::WellEconProductionLimits::QuantityLimit::POTN,
                                       6.0,
                                       Ewoms::WellEconProductionLimits::EconWorkover::WELL,
                                       7.0, 8.0, 9.0, 10.0);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellEconProductionLimits)
#endif
}

BOOST_AUTO_TEST_CASE(WellGuideRate)
{
#ifdef HAVE_MPI
    Ewoms::Well::WellGuideRate val1{true, 1.0, Ewoms::Well::GuideRateTarget::COMB, 2.0};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Well::WellGuideRate)
#endif
}

BOOST_AUTO_TEST_CASE(WellConnections)
{
#ifdef HAVE_MPI
    Ewoms::Connection conn(Ewoms::Connection::Direction::Y,
                         1.0, Ewoms::Connection::State::SHUT,
                         2, 3, 4.0, 5.0, 6.0, 7.0, 8.0,
                         {9, 10, 11}, 12345, Ewoms::Connection::CTFKind::Defaulted,
                         12, 13.0, 14.0, true,
                         15, 16, 17.0);
    Ewoms::WellConnections val1(Ewoms::Connection::Order::TRACK, 1, 2, {conn, conn});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellConnections)
#endif
}

BOOST_AUTO_TEST_CASE(WellProductionProperties)
{
#ifdef HAVE_MPI
    Ewoms::Well::WellProductionProperties val1("test",
                                             Ewoms::UDAValue(1.0),
                                             Ewoms::UDAValue("test"),
                                             Ewoms::UDAValue(2.0),
                                             Ewoms::UDAValue(3.0),
                                             Ewoms::UDAValue(4.0),
                                             Ewoms::UDAValue(5.0),
                                             Ewoms::UDAValue(6.0),
                                             5.0, 6.0,
                                             7.0, 8.0,
                                             9,
                                             10.0,
                                             true,
                                             Ewoms::Well::ProducerCMode::CRAT,
                                             Ewoms::Well::ProducerCMode::BHP, 11);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Well::WellProductionProperties)
#endif
}

BOOST_AUTO_TEST_CASE(SpiralICD)
{
#ifdef HAVE_MPI
    Ewoms::SpiralICD val1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8, 9.0,
                        Ewoms::ICDStatus::OPEN, 10.0);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SpiralICD)
#endif
}

BOOST_AUTO_TEST_CASE(Valve)
{
#ifdef HAVE_MPI
    Ewoms::Valve val1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, Ewoms::ICDStatus::OPEN);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Valve)
#endif
}

BOOST_AUTO_TEST_CASE(Segment)
{
#ifdef HAVE_MPI
    Ewoms::Segment val1(1, 2, 3, {1, 2}, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, false,
                      Ewoms::Segment::SegmentType::SICD,
                      std::make_shared<Ewoms::SpiralICD>(),
                      std::make_shared<Ewoms::Valve>());
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Segment)
#endif
}

BOOST_AUTO_TEST_CASE(Dimension)
{
#ifdef HAVE_MPI
    Ewoms::Dimension val1(1.0, 2.0);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Dimension)
#endif
}

BOOST_AUTO_TEST_CASE(UnitSystem)
{
#ifdef HAVE_MPI
    Ewoms::UnitSystem val1(Ewoms::UnitSystem::UnitType::UNIT_TYPE_METRIC);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UnitSystem)
#endif
}

BOOST_AUTO_TEST_CASE(WellSegments)
{
#ifdef HAVE_MPI
    Ewoms::Segment seg(1, 2, 3, {1, 2}, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, false,
                     Ewoms::Segment::SegmentType::SICD,
                     std::make_shared<Ewoms::SpiralICD>(),
                     std::make_shared<Ewoms::Valve>());
    Ewoms::WellSegments val1(Ewoms::WellSegments::CompPressureDrop::HF_,
                           {seg, seg});

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellSegments)
#endif
}

BOOST_AUTO_TEST_CASE(Well)
{
#ifdef HAVE_MPI
    Ewoms::Well val1 = getFullWell();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Well)
#endif
}

BOOST_AUTO_TEST_CASE(GroupInjectionProperties)
{
#ifdef HAVE_MPI
    Ewoms::Group::GroupInjectionProperties val1{Ewoms::Phase::WATER,
                                              Ewoms::Group::InjectionCMode::REIN,
                                              Ewoms::UDAValue(1.0),
                                              Ewoms::UDAValue(2.0),
                                              Ewoms::UDAValue(3.0),
                                              Ewoms::UDAValue(4.0),
                                              "test1", "test2", 5};

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Group::GroupInjectionProperties)
#endif
}

BOOST_AUTO_TEST_CASE(GroupProductionProperties)
{
#ifdef HAVE_MPI
    Ewoms::Group::GroupProductionProperties val1{Ewoms::Group::ProductionCMode::PRBL,
                                               Ewoms::Group::ExceedAction::WELL,
                                               Ewoms::UDAValue(1.0),
                                               Ewoms::UDAValue(2.0),
                                               Ewoms::UDAValue(3.0),
                                               Ewoms::UDAValue(4.0),
                                               5.0, Ewoms::Group::GuideRateTarget::COMB,
                                               6.0, 7};

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Group::GroupProductionProperties)
#endif
}

BOOST_AUTO_TEST_CASE(Group)
{
#ifdef HAVE_MPI
    Ewoms::UnitSystem unitSystem;

    std::map<Ewoms::Phase, Ewoms::Group::GroupInjectionProperties> injection;
    Ewoms::Group val1("test1", 1, 2, 3.0, unitSystem,
                    Ewoms::Group::GroupType::PRODUCTION,
                    4.0, true, false, 5, "test2",
                    Ewoms::IOrderSet<std::string>({"test3", "test4"}, {"test3","test4"}),
                    Ewoms::IOrderSet<std::string>({"test5", "test6"}, {"test5","test6"}),
                    injection,
                    Ewoms::Group::GroupProductionProperties());

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Group)
#endif
}

BOOST_AUTO_TEST_CASE(WList)
{
#ifdef HAVE_MPI
    Ewoms::WList val1({"test1", "test2", "test3"});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WList)
#endif
}

BOOST_AUTO_TEST_CASE(WListManager)
{
#ifdef HAVE_MPI
    Ewoms::WList wl({"test1", "test2", "test3"});
    std::map<std::string,Ewoms::WList> data{{"test", wl}, {"test2", wl}};
    Ewoms::WListManager val1(data);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WListManager)
#endif
}

BOOST_AUTO_TEST_CASE(UDQASTNode)
{
#ifdef HAVE_MPI
  std::shared_ptr<Ewoms::UDQASTNode> n0;
  std::shared_ptr<Ewoms::UDQASTNode> n1(new Ewoms::UDQASTNode(Ewoms::UDQVarType::NONE,
                                                          Ewoms::UDQTokenType::error,
                                                          "test1", 1.0, {"test2"},
                                                          n0, n0));
    Ewoms::UDQASTNode val1(Ewoms::UDQVarType::NONE,
                         Ewoms::UDQTokenType::error,
                         "test", 1.0, {"test3"}, n1, n1);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDQASTNode)
#endif
}

BOOST_AUTO_TEST_CASE(UDQDefine)
{
#ifdef HAVE_MPI
    std::shared_ptr<Ewoms::UDQASTNode> n0;
    Ewoms::UDQASTNode n1(Ewoms::UDQVarType::NONE,
                       Ewoms::UDQTokenType::error,
                       "test", 1.0, {"test1", "test2"}, n0, n0);
    Ewoms::UDQDefine val1("test", std::make_shared<Ewoms::UDQASTNode>(n1),
                        Ewoms::UDQVarType::NONE, "test2");
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDQDefine)
#endif
}

BOOST_AUTO_TEST_CASE(UDQAssign)
{
#ifdef HAVE_MPI
    Ewoms::UDQAssign val1("test", Ewoms::UDQVarType::NONE,
                        {Ewoms::UDQAssign::AssignRecord{{"test1"}, 1.0},
                         Ewoms::UDQAssign::AssignRecord{{"test2"}, 2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDQAssign)
#endif
}

BOOST_AUTO_TEST_CASE(UDQIndex)
{
#ifdef HAVE_MPI
    Ewoms::UDQIndex val1(1, 2, Ewoms::UDQAction::ASSIGN, Ewoms::UDQVarType::WELL_VAR);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDQIndex)
#endif
}

BOOST_AUTO_TEST_CASE(UDQConfig)
{
#ifdef HAVE_MPI
    Ewoms::UDQConfig val1 = getUDQConfig();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDQConfig)
#endif
}

BOOST_AUTO_TEST_CASE(UDQActiveInputRecord)
{
#ifdef HAVE_MPI
    Ewoms::UDQActive::InputRecord val1(1, "test1", "test2",
                                     Ewoms::UDAControl::WCONPROD_ORAT);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDQActive::InputRecord)
#endif
}

BOOST_AUTO_TEST_CASE(UDQActiveRecord)
{
#ifdef HAVE_MPI
    Ewoms::UDQActive::Record val1("test1", 1, 2, "test2",
                                Ewoms::UDAControl::WCONPROD_ORAT);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDQActive::Record)
#endif
}

BOOST_AUTO_TEST_CASE(UDQActive)
{
#ifdef HAVE_MPI
    Ewoms::UDQActive val1({Ewoms::UDQActive::InputRecord(1, "test1", "test2",
                                                     Ewoms::UDAControl::WCONPROD_ORAT)},
                        {Ewoms::UDQActive::Record("test1", 1, 2, "test2",
                                                  Ewoms::UDAControl::WCONPROD_ORAT)},
                        {{"test1", 1}}, {{"test2", 2}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(UDQActive)
#endif
}

BOOST_AUTO_TEST_CASE(AquiferCT)
{
#ifdef HAVE_MPI
    Ewoms::AquiferCT val1 = getAquiferCT();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(AquiferCT);
#endif
}

BOOST_AUTO_TEST_CASE(Aquifetp)
{
#ifdef HAVE_MPI
    Ewoms::Aquifetp val1 = getAquifetp();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Aquifetp);
#endif
}

BOOST_AUTO_TEST_CASE(Aquancon)
{
#ifdef HAVE_MPI
    Ewoms::Aquancon val1 = getAquancon();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Aquancon);
#endif
}

BOOST_AUTO_TEST_CASE(AquferConfig)
{
#ifdef HAVE_MPI
    Ewoms::Aquifetp fetp = getAquifetp();
    Ewoms::AquiferCT ct = getAquiferCT();
    Ewoms::Aquancon conn = getAquancon();
    Ewoms::AquiferConfig val1(fetp, ct, conn);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(AquiferConfig);
#endif
}

BOOST_AUTO_TEST_CASE(GuideRateModel)
{
#ifdef HAVE_MPI
    Ewoms::GuideRateModel val1 = getGuideRateModel();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GuideRateModel)
#endif
}

BOOST_AUTO_TEST_CASE(GuideRateConfigGroup)
{
#ifdef HAVE_MPI
    Ewoms::GuideRateConfig::GroupTarget val1 = getGuideRateConfigGroup();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GuideRateConfig::GroupTarget)
#endif
}

BOOST_AUTO_TEST_CASE(GuideRateConfigWell)
{
#ifdef HAVE_MPI
    Ewoms::GuideRateConfig::WellTarget val1 = getGuideRateConfigWell();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GuideRateConfig::WellTarget)
#endif
}

BOOST_AUTO_TEST_CASE(GuideRateConfig)
{
#ifdef HAVE_MPI
    auto model = std::make_shared<Ewoms::GuideRateModel>(getGuideRateModel());
    Ewoms::GuideRateConfig val1(model,
                              {{"test1", getGuideRateConfigWell()}},
                              {{"test2", getGuideRateConfigGroup()}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GuideRateConfig)
#endif
}

BOOST_AUTO_TEST_CASE(GConSaleGroup)
{
#ifdef HAVE_MPI
    Ewoms::GConSale::GCONSALEGroup val1{Ewoms::UDAValue(1.0),
                                      Ewoms::UDAValue(2.0),
                                      Ewoms::UDAValue(3.0),
                                      Ewoms::GConSale::MaxProcedure::PLUG,
                                      4.0, Ewoms::UnitSystem()};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GConSale::GCONSALEGroup)
#endif
}

BOOST_AUTO_TEST_CASE(GConSale)
{
#ifdef HAVE_MPI
    Ewoms::GConSale::GCONSALEGroup group{Ewoms::UDAValue(1.0),
                                       Ewoms::UDAValue(2.0),
                                       Ewoms::UDAValue(3.0),
                                       Ewoms::GConSale::MaxProcedure::PLUG,
                                       4.0, Ewoms::UnitSystem()};
    Ewoms::GConSale val1({{"test1", group}, {"test2", group}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GConSale)
#endif
}

BOOST_AUTO_TEST_CASE(GConSumpGroup)
{
#ifdef HAVE_MPI
    Ewoms::GConSump::GCONSUMPGroup val1{Ewoms::UDAValue(1.0),
                                      Ewoms::UDAValue(2.0),
                                      "test",
                                      3.0, Ewoms::UnitSystem()};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GConSump::GCONSUMPGroup)
#endif
}

BOOST_AUTO_TEST_CASE(GConSump)
{
#ifdef HAVE_MPI
    Ewoms::GConSump::GCONSUMPGroup group{Ewoms::UDAValue(1.0),
                                       Ewoms::UDAValue(2.0),
                                       "test",
                                       3.0, Ewoms::UnitSystem()};
    Ewoms::GConSump val1({{"test1", group}, {"test2", group}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GConSump)
#endif
}

BOOST_AUTO_TEST_CASE(RFTConfig)
{
#ifdef HAVE_MPI
    Ewoms::RFTConfig val1(getTimeMap(),
                        std::size_t{1729},
                        {true, 1},
                        {{"test1", 2}, {"test2", 3}},
                        {{"test3", 2}},
                        {{"test1", {{{Ewoms::RFTConfig::RFT::TIMESTEP, 3}}, 4}}},
                        {{"test2", {{{Ewoms::RFTConfig::PLT::REPT, 5}}, 6}}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(RFTConfig)
#endif
}

BOOST_AUTO_TEST_CASE(DeckItem)
{
#ifdef HAVE_MPI
    Ewoms::DeckItem val1({1.0}, {2}, {"test3"}, {Ewoms::UDAValue(4)},
                       Ewoms::type_tag::string, "test5",
                       {Ewoms::value::status::deck_value},
                       true,
                       {Ewoms::Dimension(7.0, 8.0)},
                       {Ewoms::Dimension(10.0, 11.0)});

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(DeckItem)
#endif
}

BOOST_AUTO_TEST_CASE(DeckRecord)
{
#ifdef HAVE_MPI
    Ewoms::DeckRecord val1 = getDeckRecord();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(DeckRecord)
#endif
}

BOOST_AUTO_TEST_CASE(Location)
{
#ifdef HAVE_MPI
    Ewoms::Location val1{"test", 1};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Location)
#endif
}

BOOST_AUTO_TEST_CASE(DeckKeyword)
{
#ifdef HAVE_MPI
    Ewoms::DeckKeyword val1("test", {"test",1},
                          {getDeckRecord(), getDeckRecord()}, true, false);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(DeckKeyword)
#endif
}

BOOST_AUTO_TEST_CASE(Deck)
{
#ifdef HAVE_MPI
    std::unique_ptr<Ewoms::UnitSystem> unitSys(new Ewoms::UnitSystem);
    Ewoms::Deck val1({Ewoms::DeckKeyword("test", {"test",1},
                                     {getDeckRecord(), getDeckRecord()}, true, false)},
                   Ewoms::UnitSystem(), unitSys.get(),
                   "test2", "test3", 2);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Deck)
#endif
}

BOOST_AUTO_TEST_CASE(Tuning)
{
#ifdef HAVE_MPI
    Ewoms::Tuning val1 = getTuning();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Tuning)
#endif
}

BOOST_AUTO_TEST_CASE(ASTNode)
{
#ifdef HAVE_MPI
    Ewoms::Action::ASTNode child(number, FuncType::field, "test3", {"test2"}, 2.0, {});
    Ewoms::Action::ASTNode val1(number, FuncType::field, "test1", {"test2"}, 1.0, {child});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Action::ASTNode)
#endif
}

BOOST_AUTO_TEST_CASE(AST)
{
#ifdef HAVE_MPI
    std::shared_ptr<Ewoms::Action::ASTNode> node;
    node.reset(new Ewoms::Action::ASTNode(number, FuncType::field,
                                        "test1", {"test2"}, 1.0, {}));
    Ewoms::Action::AST val1(node);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Action::AST)
#endif
}

BOOST_AUTO_TEST_CASE(Quantity)
{
#ifdef HAVE_MPI
    Ewoms::Action::Quantity val1;
    val1.quantity = "test1";
    val1.args = {"test2", "test3"};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Action::Quantity)
#endif
}

BOOST_AUTO_TEST_CASE(Condition)
{
#ifdef HAVE_MPI
    Ewoms::Action::Condition val1 = getCondition();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Action::Condition)
#endif
}

BOOST_AUTO_TEST_CASE(ActionX)
{
#ifdef HAVE_MPI
    Ewoms::Action::ActionX val1 = getActionX();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Action::ActionX)
#endif
}

BOOST_AUTO_TEST_CASE(PyAction)
{
#ifdef HAVE_MPI
    Ewoms::Action::PyAction val1("name", Ewoms::Action::PyAction::RunCount::single, "import opm");
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Action::PyAction)
#endif
}

BOOST_AUTO_TEST_CASE(Actions)
{
#ifdef HAVE_MPI
    Ewoms::Action::Actions val1({getActionX()}, {{"name", Ewoms::Action::PyAction::RunCount::unlimited, "import numpy"}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Action::Actions)
#endif
}

BOOST_AUTO_TEST_CASE(Schedule)
{
#ifdef HAVE_MPI
    Ewoms::UnitSystem unitSystem;
    Ewoms::Schedule::WellMap wells;
    wells.insert({"test", {{std::make_shared<Ewoms::Well>(getFullWell())},1}});
    Ewoms::Schedule::GroupMap groups;
    std::map<Ewoms::Phase, Ewoms::Group::GroupInjectionProperties> injection;
    groups.insert({"test", {{std::make_shared<Ewoms::Group>("test1", 1, 2, 3.0, unitSystem,
                                                          Ewoms::Group::GroupType::PRODUCTION,
                                                          4.0, true, false, 5, "test2",
                                                          Ewoms::IOrderSet<std::string>({"test3", "test4"}, {"test3","test4"}),
                                                          Ewoms::IOrderSet<std::string>({"test5", "test6"}, {"test5","test6"}),
                                                          injection,
                                                          Ewoms::Group::GroupProductionProperties())},1}});
    using VapType = Ewoms::OilVaporizationProperties::OilVaporization;
    Ewoms::DynamicState<Ewoms::OilVaporizationProperties> oilvap{{Ewoms::OilVaporizationProperties(VapType::VAPPARS,
                                                                                   1.0, 2.0, {5.0, 6.0},
                                                                                   {false, true}, {7.0, 8.0})},1};
    Ewoms::Events events(Ewoms::DynamicVector<uint64_t>({1,2,3,4,5}));
    std::unique_ptr<Ewoms::UnitSystem> unitSys(new Ewoms::UnitSystem);
    Ewoms::Deck sdeck({Ewoms::DeckKeyword("test", {"test",1},
                                     {getDeckRecord(), getDeckRecord()}, true, false)},
                   Ewoms::UnitSystem(), unitSys.get(),
                   "test2", "test3", 2);
    Ewoms::DynamicVector<Ewoms::Deck> modifierDeck({sdeck});
    std::vector<Ewoms::MLimits> limits{Ewoms::MLimits{1,2,3,4,5,6,7,8,9,10,11,12}};
    Ewoms::MessageLimits messageLimits(Ewoms::DynamicState<Ewoms::MLimits>(limits,2));
    Ewoms::Runspec runspec(Ewoms::Phases(true, true, true, false, true, false, true, false),
                         Ewoms::Tabdims(1,2,3,4,5,6),
                         Ewoms::EndpointScaling(std::bitset<4>(13)),
                         Ewoms::Welldims(1,2,3,4),
                         Ewoms::WellSegmentDims(1,2,3),
                         Ewoms::UDQParams(true, 1, 2.0, 3.0, 4.0),
                         Ewoms::EclHysterConfig(true, 1, 2),
                         Ewoms::Actdims(1,2,3,4),
                         Ewoms::SatFuncControls(5.6e-7, Ewoms::SatFuncControls::ThreePhaseOilKrModel::Stone1));
    Ewoms::Schedule::VFPProdMap vfpProd {{1, {{std::make_shared<Ewoms::VFPProdTable>(getVFPProdTable())},1}}};
    Ewoms::Schedule::VFPInjMap vfpIn{{1, {{std::make_shared<Ewoms::VFPInjTable>(getVFPInjTable())},1}}};
    Ewoms::WellTestConfig::WTESTWell tw{"test", Ewoms::WellTestConfig::ECONOMIC,
                                         1.0, 2, 3.0, 4};
    Ewoms::WellTestConfig wtc({tw, tw, tw});

    Ewoms::WList wl({"test1", "test2", "test3"});
    std::map<std::string,Ewoms::WList> data{{"test", wl}, {"test2", wl}};
    Ewoms::WListManager wlm(data);

    Ewoms::UDQActive udqa({Ewoms::UDQActive::InputRecord(1, "test1", "test2",
                                                     Ewoms::UDAControl::WCONPROD_ORAT)},
                        {Ewoms::UDQActive::Record("test1", 1, 2, "test2",
                                                  Ewoms::UDAControl::WCONPROD_ORAT)},
                        {{"test1", 1}}, {{"test2", 2}});

    auto model = std::make_shared<Ewoms::GuideRateModel>(getGuideRateModel());
    Ewoms::GuideRateConfig grc(model,
                             {{"test1", getGuideRateConfigWell()}},
                             {{"test2", getGuideRateConfigGroup()}});

    Ewoms::GConSale::GCONSALEGroup group{Ewoms::UDAValue(1.0),
                                       Ewoms::UDAValue(2.0),
                                       Ewoms::UDAValue(3.0),
                                       Ewoms::GConSale::MaxProcedure::PLUG,
                                       4.0, Ewoms::UnitSystem()};
    Ewoms::GConSale gcs({{"test1", group}, {"test2", group}});

    Ewoms::GConSump::GCONSUMPGroup grp{Ewoms::UDAValue(1.0),
                                     Ewoms::UDAValue(2.0),
                                     "test",
                                     3.0, Ewoms::UnitSystem()};
    Ewoms::GConSump gcm({{"test1", grp}, {"test2", grp}});

    Ewoms::Action::Actions actions({getActionX()}, {{"pyaction", Ewoms::Action::PyAction::RunCount::single, "import os"}});

    Ewoms::RFTConfig rftc(getTimeMap(),
                        std::size_t{1729},
                        {true, 1},
                        {{"test1", 2}, {"test2", 3}},
                        {{"test3", 2}},
                        {{"test1", {{{Ewoms::RFTConfig::RFT::TIMESTEP, 3}}, 4}}},
                        {{"test2", {{{Ewoms::RFTConfig::PLT::REPT, 5}}, 6}}});

    Ewoms::Schedule val1(getTimeMap(),
                       wells,
                       groups,
                       oilvap,
                       events,
                       modifierDeck,
                       Ewoms::DynamicState<Ewoms::Tuning>({getTuning()}, 1),
                       messageLimits,
                       runspec,
                       vfpProd,
                       vfpIn,
                       {{std::make_shared<Ewoms::WellTestConfig>(wtc)}, 1},
                       {{std::make_shared<Ewoms::WListManager>(wlm)}, 1},
                       {{std::make_shared<Ewoms::UDQConfig>(getUDQConfig())}, 1},
                       {{std::make_shared<Ewoms::UDQActive>(udqa)}, 1},
                       {{std::make_shared<Ewoms::GuideRateConfig>(grc)}, 1},
                       {{std::make_shared<Ewoms::GConSale>(gcs)}, 1},
                       {{std::make_shared<Ewoms::GConSump>(gcm)}, 1},
                       {{Ewoms::Well::ProducerCMode::CRAT}, 1},
                       {{std::make_shared<Ewoms::Action::Actions>(actions)}, 1},
                       rftc,
                       {std::vector<int>{1}, 1},
                       getRestartConfig(),
                       {{"test", events}});

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Schedule)
#endif
}

BOOST_AUTO_TEST_CASE(BrineDensityTable)
{
#ifdef HAVE_MPI
    Ewoms::BrineDensityTable val1({1.0, 2.0, 3.0});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(BrineDensityTable)
#endif
}

BOOST_AUTO_TEST_CASE(SummaryConfigNode)
{
#ifdef HAVE_MPI
    auto val1 = Ewoms::SummaryConfigNode{"test1", Ewoms::SummaryConfigNode::Category::Region,
                                 Ewoms::Location{"test2", 1}}
                                 .parameterType(Ewoms::SummaryConfigNode::Type::Pressure)
                                 .namedEntity("test3")
                                 .number(2)
                                 .isUserDefined(true);

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SummaryConfigNode)
#endif
}

BOOST_AUTO_TEST_CASE(SummaryConfig)
{
#ifdef HAVE_MPI
    auto node = Ewoms::SummaryConfigNode{"test1", Ewoms::SummaryConfigNode::Category::Region,
                                 Ewoms::Location{"test2", 1}}
                                 .parameterType(Ewoms::SummaryConfigNode::Type::Pressure)
                                 .namedEntity("test3")
                                 .number(2)
                                 .isUserDefined(true);
    Ewoms::SummaryConfig val1({node}, {"test1", "test2"}, {"test3", "test4"});

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SummaryConfig)
#endif
}

BOOST_AUTO_TEST_CASE(PvtwsaltTable)
{
#ifdef HAVE_MPI
    Ewoms::PvtwsaltTable val1(1.0, 2.0, {3.0, 4.0, 5.0});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(PvtwsaltTable)
#endif
}

BOOST_AUTO_TEST_CASE(WellBrineProperties)
{
#ifdef HAVE_MPI
    Ewoms::WellBrineProperties val1{1.0};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellBrineProperties)
#endif
}

BOOST_AUTO_TEST_CASE(MULTREGTRecord)
{
#ifdef HAVE_MPI
    Ewoms::MULTREGTRecord val1{1, 2, 3.0, 4, Ewoms::MULTREGT::ALL, "test"};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(MULTREGTRecord)
#endif
}

BOOST_AUTO_TEST_CASE(MULTREGTScanner)
{
#ifdef HAVE_MPI
    std::vector<Ewoms::MULTREGTRecord> records{{1, 2, 3.0, 4, Ewoms::MULTREGT::ALL, "test1"}};
    std::map<std::pair<int, int>, int> searchRecord{{{5,6},0}};
    Ewoms::MULTREGTScanner::ExternalSearchMap searchMap;
    searchMap.insert({"test2", searchRecord});
    Ewoms::MULTREGTScanner val1({1, 2, 3},
                              records,
                              searchMap,
                              {{"test3", {7,8}}},
                              "test4");

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(MULTREGTScanner)
#endif
}

BOOST_AUTO_TEST_CASE(EclipseConfig)
{
#ifdef HAVE_MPI
    Ewoms::IOConfig io(true, false, true, false, false, true, "test1", true,
                     "test2", true, "test3", false);
    Ewoms::InitConfig init(Ewoms::Equil({getEquilRecord(), getEquilRecord()}),
                         getFoamConfig(),
                         true, true, true, 20, "test1");
    Ewoms::EclipseConfig val1{init, io};

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(EclipseConfig)
#endif
}

BOOST_AUTO_TEST_CASE(TransMult)
{
#ifdef HAVE_MPI
    std::vector<Ewoms::MULTREGTRecord> records{{1, 2, 3.0, 4, Ewoms::MULTREGT::ALL, "test1"}};
    std::map<std::pair<int, int>, int> searchRecord{{{5,6},0}};
    Ewoms::MULTREGTScanner::ExternalSearchMap searchMap;
    searchMap.insert({"test2", searchRecord});
    Ewoms::MULTREGTScanner scanner({1, 2, 3},
                                 records,
                                 searchMap,
                                 {{"test3", {7,8}}},
                                 "test4");

    Ewoms::TransMult val1({1, 2, 3},
                        {{Ewoms::FaceDir::YPlus, {4.0, 5.0}}},
                        {{Ewoms::FaceDir::ZPlus, "test1"}},
                        scanner);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(TransMult)
#endif
}

BOOST_AUTO_TEST_CASE(FaultFace)
{
#ifdef HAVE_MPI
    Ewoms::FaultFace val1({1,2,3,4,5,6}, Ewoms::FaceDir::YPlus);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(FaultFace)
#endif
}

BOOST_AUTO_TEST_CASE(Fault)
{
#ifdef HAVE_MPI
    Ewoms::Fault val1("test", 1.0, {{{1,2,3,4,5,6}, Ewoms::FaceDir::YPlus}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Fault)
#endif
}

BOOST_AUTO_TEST_CASE(WellType)
{
#ifdef HAVE_MPI
    Ewoms::WellType val1(true, Ewoms::Phase::OIL);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(WellType)
#endif
}

BOOST_AUTO_TEST_CASE(DenT)
{
#ifdef HAVE_MPI
    Ewoms::DenT val1 = getDenT();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(DenT)
#endif
}

BOOST_AUTO_TEST_CASE(FaultCollection)
{
#ifdef HAVE_MPI
    Ewoms::Fault fault("test", 1.0, {{{1,2,3,4,5,6}, Ewoms::FaceDir::YPlus}});
    Ewoms::OrderedMap<std::string, Ewoms::Fault> faults;
    faults.insert({"test2", fault});
    Ewoms::FaultCollection val1(faults);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(FaultCollection)
#endif
}

BOOST_AUTO_TEST_CASE(SolventDensityTable)
{
#ifdef HAVE_MPI
    Ewoms::SolventDensityTable val1({1.0, 2.0, 3.0});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SolventDensityTable)
#endif
}

BOOST_AUTO_TEST_CASE(GridDims)
{
#ifdef HAVE_MPI
    Ewoms::GridDims val1{ 1,  2,  3};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GridDims)
#endif
}

BOOST_AUTO_TEST_CASE(PlyshlogTable)
{
#ifdef HAVE_MPI
    Ewoms::OrderedMap<std::string, Ewoms::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    Ewoms::PlyshlogTable val1(getTableSchema(), data, true, 1.0, 2.0, 3.0, true, true);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(PlyshlogTable)
#endif
}

BOOST_AUTO_TEST_CASE(RocktabTable)
{
#ifdef HAVE_MPI
    Ewoms::OrderedMap<std::string, Ewoms::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    Ewoms::RocktabTable val1(getTableSchema(), data, true, true);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(RocktabTable)
#endif
}

bool init_unit_test_func()
{
    return true;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
