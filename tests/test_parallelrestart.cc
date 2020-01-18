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
#include <ewoms/material/fluidsystems/blackoilpvt/drygaspvt.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/solventpvt.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/wetgaspvt.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/deck/deckitem.hh>
#include <ewoms/eclio/parser/eclipsestate/runspec.hh>
#include <ewoms/eclio/parser/eclipsestate/edit/editnnc.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/nnc.hh>
#include <ewoms/eclio/parser/eclipsestate/initconfig/equil.hh>
#include <ewoms/eclio/parser/eclipsestate/initconfig/foamconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/initconfig/initconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/ioconfig/ioconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/ioconfig/restartconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actionast.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actions.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actionx.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/astnode.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/condition.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/events.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/gconsale.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/group.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/guideratemodel.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/messagelimits.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/spiralicd.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/valve.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/oilvaporizationproperties.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/rftconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
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
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/simulationconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/thresholdpressure.hh>
#include <ewoms/eclio/parser/eclipsestate/summaryconfig/summaryconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/aqudims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/columnschema.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/eqldims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/flattable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/jfunc.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/plymwinjtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/pvtgtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/pvtotable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/regdims.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rock2dtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rock2dtrtable.hh>
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
    return well1;
}

Ewoms::ThresholdPressure getThresholdPressure()
{
    return Ewoms::ThresholdPressure(false, true, {{true, 1.0}, {false, 2.0}},
                                  {{{1,2},{false,3.0}},{{2,3},{true,4.0}}});
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

Ewoms::TimeMap getTimeMap()
{
    return Ewoms::TimeMap({123},
                        {{1, Ewoms::TimeStampUTC(123)}},
                        {{2, Ewoms::TimeStampUTC(456)}});
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
    result.addTable(0, std::make_shared<const Ewoms::SimpleTable>(tab1));
    result.addTable(1, std::make_shared<const Ewoms::SimpleTable>(tab1));
    return result;
}

Ewoms::Well getFullWell()
{
    Ewoms::UnitSystem unitSystem;
    return Ewoms::Well("test1", "test2", 1, 2, 3, 4, 5.0,
                     Ewoms::Phase::WATER, Ewoms::Connection::Order::DEPTH,
                     unitSystem, 6.0, Ewoms::Well::Status::SHUT,
                     7.0, true, true, false,
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
    Ewoms::VFPInjTable::extents shape;
    shape[0] = 3;
    shape[1] = 2;
    table.resize(shape);
    double foo = 1.0;
    for (size_t i = 0; i < table.num_elements(); ++i)
        *(table.data() + i) = foo++;
    return Ewoms::VFPInjTable(1, 2.0, Ewoms::VFPInjTable::FLO_WAT, {1.0, 2.0},
                            {3.0, 4.0, 5.0}, table);
}

Ewoms::VFPProdTable getVFPProdTable()
{
    Ewoms::VFPProdTable::array_type table;
    Ewoms::VFPProdTable::extents shape;
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 3;
    shape[3] = 4;
    shape[4] = 5;
    table.resize(shape);
    double foo = 1.0;
    for (size_t i = 0; i < table.num_elements(); ++i)
        *(table.data() + i) = foo++;
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
                       {Ewoms::Dimension("DimensionLess", 7.0, 8.0)},
                       {Ewoms::Dimension("Metric", 10.0, 11.0)});

    Ewoms::DeckItem item2({1.0}, {2}, {"test3"}, {Ewoms::UDAValue(4)},
                       Ewoms::type_tag::string, "test6",
                       {Ewoms::value::status::deck_value},
                       true,
                       {Ewoms::Dimension("DimensionLess", 7.0, 8.0)},
                       {Ewoms::Dimension("Metric", 10.0, 11.0)});

    return Ewoms::DeckRecord({item1, item2});
}

Ewoms::Tuning getTuning()
{
    return Ewoms::Tuning(Ewoms::DynamicState<double>(std::vector<double>{1.0}, 1),  //TSINIT
                       Ewoms::DynamicState<double>(std::vector<double>{2.0}, 1),  //TSMAXZ
                       Ewoms::DynamicState<double>(std::vector<double>{3.0}, 1),  //TSMINZ
                       Ewoms::DynamicState<double>(std::vector<double>{4.0}, 1),  //TSMCHP
                       Ewoms::DynamicState<double>(std::vector<double>{5.0}, 1),  //TSFMAX
                       Ewoms::DynamicState<double>(std::vector<double>{6.0}, 1),  //TSFMIN
                       Ewoms::DynamicState<double>(std::vector<double>{7.0}, 1),  //TSFCNV
                       Ewoms::DynamicState<double>(std::vector<double>{8.0}, 1),  //TFDIFF
                       Ewoms::DynamicState<double>(std::vector<double>{9.0}, 1),  //THRUPT
                       Ewoms::DynamicState<double>(std::vector<double>{10.0}, 1), //TMAXWC
                       Ewoms::DynamicState<int>(std::vector<int>{1}, 1),       //TMAXWC_has_value
                       Ewoms::DynamicState<double>(std::vector<double>{11.0}, 1), //TRGTTE
                       Ewoms::DynamicState<double>(std::vector<double>{12.0}, 1), //TRGCNV
                       Ewoms::DynamicState<double>(std::vector<double>{13.0}, 1), //TRGMBE
                       Ewoms::DynamicState<double>(std::vector<double>{14.0}, 1), //TRGLCV
                       Ewoms::DynamicState<double>(std::vector<double>{15.0}, 1), //XXXTTE
                       Ewoms::DynamicState<double>(std::vector<double>{16.0}, 1), //XXXCNV
                       Ewoms::DynamicState<double>(std::vector<double>{17.0}, 1), //XXXMBE
                       Ewoms::DynamicState<double>(std::vector<double>{18.0}, 1), //XXXLCV
                       Ewoms::DynamicState<double>(std::vector<double>{19.0}, 1), //XXXWFL
                       Ewoms::DynamicState<double>(std::vector<double>{20.0}, 1), ///TRGFIP
                       Ewoms::DynamicState<double>(std::vector<double>{21.0}, 1), //TRGSFT
                       Ewoms::DynamicState<int>(std::vector<int>{2}, 1),       //TRGSFT_has_value
                       Ewoms::DynamicState<double>(std::vector<double>{22.0}, 1), // THIONX
                       Ewoms::DynamicState<int>(std::vector<int>{3}, 1),       //TRWGHT
                       Ewoms::DynamicState<int>(std::vector<int>{4}, 1),       //NEWTMX
                       Ewoms::DynamicState<int>(std::vector<int>{5}, 1),       //NEWTMN
                       Ewoms::DynamicState<int>(std::vector<int>{6}, 1),       //LITMAX
                       Ewoms::DynamicState<int>(std::vector<int>{7}, 1),       //LITMIN
                       Ewoms::DynamicState<int>(std::vector<int>{8}, 1),       //MXWSIT
                       Ewoms::DynamicState<int>(std::vector<int>{9}, 1),       //MXWPIT
                       Ewoms::DynamicState<double>(std::vector<double>{23.0}, 1), //DDPLIM
                       Ewoms::DynamicState<double>(std::vector<double>{24.0}, 1), //DDSLIM
                       Ewoms::DynamicState<double>(std::vector<double>{25.0}, 1), //TGRDPR
                       Ewoms::DynamicState<double>(std::vector<double>{26.0}, 1), //XXXDPR
                       Ewoms::DynamicState<int>(std::vector<int>{10}, 1),      //XXDPR_has_value
                       std::map<std::string,bool>{{"test", false}}); // resetValue
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

BOOST_AUTO_TEST_CASE(Solution)
{
#if HAVE_MPI
    Ewoms::data::Solution sol1 = getSolution();
    auto sol2 = PackUnpack(sol1);
    BOOST_CHECK(std::get<1>(sol2) == std::get<2>(sol2));
    BOOST_CHECK(sol1 == std::get<0>(sol2));
#endif
}

BOOST_AUTO_TEST_CASE(Rates)
{
#if HAVE_MPI
    Ewoms::data::Rates rat1 = getRates();
    auto rat2 = PackUnpack(rat1);
    BOOST_CHECK(std::get<1>(rat2) == std::get<2>(rat2));
    BOOST_CHECK(rat1 == std::get<0>(rat2));
#endif
}

BOOST_AUTO_TEST_CASE(dataConnection)
{
#if HAVE_MPI
    Ewoms::data::Connection con1 = getConnection();
    auto con2 = PackUnpack(con1);
    BOOST_CHECK(std::get<1>(con2) == std::get<2>(con2));
    BOOST_CHECK(con1 == std::get<0>(con2));
#endif
}

BOOST_AUTO_TEST_CASE(dataSegment)
{
#if HAVE_MPI
    Ewoms::data::Segment seg1 = getSegment();
    auto seg2 = PackUnpack(seg1);
    BOOST_CHECK(std::get<1>(seg2) == std::get<2>(seg2));
    BOOST_CHECK(seg1 == std::get<0>(seg2));
#endif
}

BOOST_AUTO_TEST_CASE(dataWell)
{
#if HAVE_MPI
    Ewoms::data::Well well1 = getWell();
    auto well2 = PackUnpack(well1);
    BOOST_CHECK(std::get<1>(well2) == std::get<2>(well2));
    BOOST_CHECK(well1 == std::get<0>(well2));
#endif
}

BOOST_AUTO_TEST_CASE(WellRates)
{
#if HAVE_MPI
    Ewoms::data::WellRates wells1;
    wells1.insert({"test_well", getWell()});
    auto wells2 = PackUnpack(wells1);
    BOOST_CHECK(std::get<1>(wells2) == std::get<2>(wells2));
    BOOST_CHECK(wells1 == std::get<0>(wells2));
#endif
}

BOOST_AUTO_TEST_CASE(CellData)
{
#if HAVE_MPI
    Ewoms::data::CellData data1;
    data1.dim = Ewoms::UnitSystem::measure::length;
    data1.data = {1.0, 2.0, 3.0};
    data1.target = Ewoms::data::TargetType::RESTART_SOLUTION;
    auto data2 = PackUnpack(data1);
    BOOST_CHECK(std::get<1>(data2) == std::get<2>(data2));
    BOOST_CHECK(data1 == std::get<0>(data2));
#endif
}

BOOST_AUTO_TEST_CASE(RestartKey)
{
#if HAVE_MPI
    Ewoms::RestartKey key1("key", Ewoms::UnitSystem::measure::length, true);
    auto key2 = PackUnpack(key1);
    BOOST_CHECK(std::get<1>(key2) == std::get<2>(key2));
    BOOST_CHECK(key1 == std::get<0>(key2));
#endif
}

BOOST_AUTO_TEST_CASE(RestartValue)
{
#if HAVE_MPI
    Ewoms::data::WellRates wells1;
    wells1.insert({"test_well", getWell()});
    Ewoms::RestartValue val1(getSolution(), wells1);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(ThresholdPressure)
{
#if HAVE_MPI
    Ewoms::ThresholdPressure val1 = getThresholdPressure();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(EDITNNC)
{
#if HAVE_MPI
    Ewoms::EDITNNC val1({{1,2,1.0},{2,3,2.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(NNC)
{
#if HAVE_MPI
    Ewoms::NNC val1({{1,2,1.0},{2,3,2.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Rock2dTable)
{
#if HAVE_MPI
    Ewoms::Rock2dTable val1({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Rock2dtrTable)
{
#if HAVE_MPI
    Ewoms::Rock2dtrTable val1({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(ColumnSchema)
{
#if HAVE_MPI
    Ewoms::ColumnSchema val1("test1", Ewoms::Table::INCREASING,
                           Ewoms::Table::DEFAULT_LINEAR);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
    Ewoms::ColumnSchema val3("test2", Ewoms::Table::DECREASING, 1.0);
    val2 = PackUnpack(val3);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val3 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(TableSchema)
{
#if HAVE_MPI
    Ewoms::TableSchema val1 = getTableSchema();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(TableColumn)
{
#if HAVE_MPI
    Ewoms::TableColumn val1 = getTableColumn();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(SimpleTable)
{
#if HAVE_MPI
    Ewoms::SimpleTable val1 = getSimpleTable();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(TableContainer)
{
#if HAVE_MPI
    Ewoms::OrderedMap<std::string, Ewoms::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    Ewoms::SimpleTable tab1(getTableSchema(), data, true);
    Ewoms::TableContainer val1(2);
    val1.addTable(0, std::make_shared<const Ewoms::SimpleTable>(tab1));
    val1.addTable(1, std::make_shared<const Ewoms::SimpleTable>(tab1));
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(EquilRecord)
{
#if HAVE_MPI
    Ewoms::EquilRecord val1 = getEquilRecord();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Equil)
{
#if HAVE_MPI
    Ewoms::Equil val1({getEquilRecord(), getEquilRecord()});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(FoamData)
{
#if HAVE_MPI
    Ewoms::FoamData val1 = getFoamData();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(FoamConfig)
{
#if HAVE_MPI
    Ewoms::FoamConfig val1({getFoamData(), getFoamData()});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(InitConfig)
{
#if HAVE_MPI
    Ewoms::InitConfig val1(Ewoms::Equil({getEquilRecord(), getEquilRecord()}),
                         Ewoms::FoamConfig({getFoamData(), getFoamData()}),
                         true, true, 20, "test1");
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(SimulationConfig)
{
#if HAVE_MPI
    Ewoms::SimulationConfig val1(getThresholdPressure(), false, true, false, true);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(RestartSchedule)
{
#if HAVE_MPI
    Ewoms::RestartSchedule val1(1, 2, 3);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(StepData)
{
#if HAVE_MPI
    Ewoms::TimeMap::StepData val1{1, Ewoms::TimeStampUTC(123456)};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(TimeMap)
{
#if HAVE_MPI
    Ewoms::TimeMap val1 = getTimeMap();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(RestartConfig)
{
#if HAVE_MPI
    Ewoms::DynamicState<Ewoms::RestartSchedule> rsched({Ewoms::RestartSchedule(1, 2, 3)}, 2);
    Ewoms::DynamicState<std::map<std::string,int>> rkw({{{"test",3}}}, 3);
    Ewoms::RestartConfig val1(getTimeMap(), 1, true, rsched, rkw, {false, true});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(IOConfig)
{
#if HAVE_MPI
    Ewoms::IOConfig val1(true, false, true, false, false, true, 1, "test1", true,
                       "test2", true, "test3", false);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Phases)
{
#if HAVE_MPI
    Ewoms::Phases val1(true, true, true, false, true, false, true, false);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Tabdims)
{
#if HAVE_MPI
    Ewoms::Tabdims val1(1,2,3,4,5,6);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(EndpointScaling)
{
#if HAVE_MPI
    Ewoms::EndpointScaling val1(std::bitset<4>(13));
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Welldims)
{
#if HAVE_MPI
    Ewoms::Welldims val1(1,2,3,4);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WellSegmentDims)
{
#if HAVE_MPI
    Ewoms::WellSegmentDims val1(1,2,3);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UDQParams)
{
#if HAVE_MPI
    Ewoms::UDQParams val1(true, 1, 2.0, 3.0, 4.0);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(EclHysterConfig)
{
#if HAVE_MPI
    Ewoms::EclHysterConfig val1(true, 1, 2);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Actdims)
{
#if HAVE_MPI
    Ewoms::Actdims val1(1,2,3,4);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
                      Ewoms::Actdims(1,2,3,4));

    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(PvtgTable)
{
#if HAVE_MPI
    Ewoms::PvtgTable val1 = getPvtgTable();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(PvtoTable)
{
#if HAVE_MPI
    Ewoms::PvtoTable val1 = getPvtoTable();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(JFunc)
{
#if HAVE_MPI
    Ewoms::JFunc val1(Ewoms::JFunc::Flag::BOTH, 1.0, 2.0,
                    3.0, 4.0, Ewoms::JFunc::Direction::XY);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(PVTWRecord)
{
#if HAVE_MPI
    Ewoms::PVTWRecord val1{1.0, 2.0, 3.0, 4.0, 5.0};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(PvtwTable)
{
#if HAVE_MPI
    Ewoms::PvtwTable val1({Ewoms::PVTWRecord{1.0, 2.0, 3.0, 4.0, 5.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(PVCDORecord)
{
#if HAVE_MPI
    Ewoms::PVTWRecord val1{1.0, 2.0, 3.0, 4.0, 5.0};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(PvcdoTable)
{
#if HAVE_MPI
    Ewoms::PvcdoTable val1({Ewoms::PVCDORecord{1.0, 2.0, 3.0, 4.0, 5.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(DENSITYRecord)
{
#if HAVE_MPI
    Ewoms::DENSITYRecord val1{1.0, 2.0, 3.0};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(DensityTable)
{
#if HAVE_MPI
    Ewoms::DensityTable val1({Ewoms::DENSITYRecord{1.0, 2.0, 3.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(VISCREFRecord)
{
#if HAVE_MPI
    Ewoms::VISCREFRecord val1{1.0, 2.0};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(ViscrefTable)
{
#if HAVE_MPI
    Ewoms::ViscrefTable val1({Ewoms::VISCREFRecord{1.0, 2.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WATDENTRecord)
{
#if HAVE_MPI
    Ewoms::WATDENTRecord val1{1.0, 2.0, 3.0};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WatdentTable)
{
#if HAVE_MPI
    Ewoms::WatdentTable val1({Ewoms::WATDENTRecord{1.0, 2.0, 3.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(PlymwinjTable)
{
#if HAVE_MPI
    Ewoms::PlymwinjTable val1({1.0}, {2.0}, 1, {{1.0}, {2.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(SkprpolyTable)
{
#if HAVE_MPI
    Ewoms::SkprpolyTable val1({1.0}, {2.0}, 1, {{1.0}, {2.0}}, 3.0);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(SkprwatTable)
{
#if HAVE_MPI
    Ewoms::SkprwatTable val1({1.0}, {2.0}, 1, {{1.0}, {2.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Regdims)
{
#if HAVE_MPI
    Ewoms::Regdims val1(1,2,3,4,5);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Eqldims)
{
#if HAVE_MPI
    Ewoms::Eqldims val1(1,2,3,4,5);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Aqudims)
{
#if HAVE_MPI
    Ewoms::Aqudims val1(1,2,3,4,5,6,7,8);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(ROCKRecord)
{
#if HAVE_MPI
    Ewoms::ROCKRecord val1{1.0,2.0};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(RockTable)
{
#if HAVE_MPI
    Ewoms::RockTable val1({Ewoms::ROCKRecord{1.0,2.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
                           Ewoms::RockTable({Ewoms::ROCKRecord{1.0,2.0}}),
                           Ewoms::ViscrefTable({Ewoms::VISCREFRecord{1.0, 2.0}}),
                           Ewoms::WatdentTable({Ewoms::WATDENTRecord{1.0, 2.0, 3.0}}),
                           {{1.0, 2.0, {1.0, 2.0, 3.0}}},
                           {{{1.0, 2.0, 3.0}}},
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
                           jfunc,
                           1.0);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(TabulatedOneDFunction)
{
#ifdef HAVE_MPI
    Ewoms::Tabulated1DFunction<double> val1(2, std::vector<double>{1.0, 2.0},
                                             std::vector<double>{3.0, 4.0});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(IntervalTabulatedTwoDFunction)
{
#ifdef HAVE_MPI
    std::vector<double> xPos{1.0, 2.0};
    std::vector<double> yPos{3.0, 4.0};
    std::vector<std::vector<double>> samples{{1.0, 2.0}, {3.0, 4.0}};
    Ewoms::IntervalTabulated2DFunction<double> val1(xPos, yPos, samples, true, true);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UniformXTabulatedTwoDFunction)
{
#ifdef HAVE_MPI
    std::vector<double> xPos{1.0, 2.0};
    std::vector<double> yPos{3.0, 4.0};
    using SampleType = std::vector<std::vector<std::tuple<double,double,double>>>;
    SampleType samples{{std::make_tuple(1.0, 2.0, 3.0)},
                       {std::make_tuple(4.0, 5.0, 6.0)}};
    using FFuncType = Ewoms::UniformXTabulated2DFunction<double>;
    FFuncType val1(xPos, yPos, samples, FFuncType::Vertical);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(SolventPvt)
{
#ifdef HAVE_MPI
    Ewoms::Tabulated1DFunction<double> func(2, std::vector<double>{1.0, 2.0},
                                             std::vector<double>{3.0, 4.0});
    Ewoms::SolventPvt<double> val1({1.0, 2.0}, {func}, {func}, {func});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(DryGasPvt)
{
#ifdef HAVE_MPI
    Ewoms::Tabulated1DFunction<double> func(2, std::vector<double>{1.0, 2.0},
                                             std::vector<double>{3.0, 4.0});
    Ewoms::DryGasPvt<double> val1({1.0, 2.0}, {func}, {func}, {func});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(GasPvtThermal)
{
#ifdef HAVE_MPI
    Ewoms::Tabulated1DFunction<double> func(2, std::vector<double>{1.0, 2.0},
                                             std::vector<double>{3.0, 4.0});
    Ewoms::GasPvtThermal<double>::IsothermalPvt* pvt = new Ewoms::GasPvtThermal<double>::IsothermalPvt;
    Ewoms::GasPvtThermal<double> val1(pvt, {func}, {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0},
                                    {func}, true, true, false);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WetGasPvt)
{
#ifdef HAVE_MPI
    Ewoms::Tabulated1DFunction<double> func(2, std::vector<double>{1.0, 2.0},
                                             std::vector<double>{3.0, 4.0});
    std::vector<double> xPos{1.0, 2.0};
    std::vector<double> yPos{3.0, 4.0};
    using FFuncType = Ewoms::UniformXTabulated2DFunction<double>;
    using Samples = std::vector<std::vector<FFuncType::SamplePoint>>;
    Samples samples({{std::make_tuple(1.0, 2.0, 3.0), std::make_tuple(3.0, 4.0, 5.0)}});
    FFuncType func2(xPos, yPos, samples, FFuncType::Vertical);
    Ewoms::WetGasPvt<double> val1({1.0, 2.0}, {3.0, 4.0},
                                {func2}, {func}, {func2},
                                {func2}, {func}, {func}, {func}, 5.0);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(ConstantCompressibilityOilPvt)
{
#ifdef HAVE_MPI
    Ewoms::ConstantCompressibilityOilPvt<double> val1({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0},
                                                    {7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(DeadOilPvt)
{
#ifdef HAVE_MPI
    Ewoms::Tabulated1DFunction<double> func(2, std::vector<double>{1.0, 2.0},
                                             std::vector<double>{3.0, 4.0});
    Ewoms::DeadOilPvt<double> val1({1.0, 2.0}, {func}, {func}, {func});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(LiveOilPvt)
{
#ifdef HAVE_MPI
    Ewoms::Tabulated1DFunction<double> func(2, std::vector<double>{1.0, 2.0},
                                             std::vector<double>{3.0, 4.0});
    std::vector<double> xPos{1.0, 2.0};
    std::vector<double> yPos{3.0, 4.0};
    using FFuncType = Ewoms::UniformXTabulated2DFunction<double>;
    using Samples = std::vector<std::vector<FFuncType::SamplePoint>>;
    Samples samples({{std::make_tuple(1.0, 2.0, 3.0), std::make_tuple(3.0, 4.0, 5.0)}});
    FFuncType func2(xPos, yPos, samples, FFuncType::Vertical);
    Ewoms::LiveOilPvt<double> val1({1.0, 2.0}, {3.0, 4.0},
                                 {func2}, {func2}, {func2},
                                 {func}, {func}, {func}, {func}, {func}, 5.0);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(OilPvtThermal)
{
#ifdef HAVE_MPI
    Ewoms::Tabulated1DFunction<double> func(2, std::vector<double>{1.0, 2.0},
                                             std::vector<double>{3.0, 4.0});
    Ewoms::OilPvtThermal<double>::IsothermalPvt* pvt = new Ewoms::OilPvtThermal<double>::IsothermalPvt;
    Ewoms::OilPvtThermal<double> val1(pvt, {func}, {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0},
                                    {7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0},
                                    {func}, true, true, false);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(ConstantCompressibilityWaterPvt)
{
#ifdef HAVE_MPI
    Ewoms::ConstantCompressibilityWaterPvt<double> val1({1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0},
                                                      {7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WaterPvtThermal)
{
#ifdef HAVE_MPI
    Ewoms::Tabulated1DFunction<double> func(2, std::vector<double>{1.0, 2.0},
                                             std::vector<double>{3.0, 4.0});
    Ewoms::WaterPvtThermal<double>::IsothermalPvt* pvt = new Ewoms::WaterPvtThermal<double>::IsothermalPvt;
    Ewoms::WaterPvtThermal<double> val1(pvt, {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0},
                                      {7.0, 8.0}, {9.0, 10.0}, {11.0, 12.0},
                                      {13.0, 14.0}, {15.0, 16.0}, {17.0, 18.0},
                                      {func}, {func}, true, true, false);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(OilVaporizationProperties)
{
#ifdef HAVE_MPI
    using VapType = Ewoms::OilVaporizationProperties::OilVaporization;
    Ewoms::OilVaporizationProperties val1(VapType::VAPPARS,
                                        {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0},
                                        {false, true}, {7.0, 8.0});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
    val1 = Ewoms::OilVaporizationProperties(VapType::DRDT,
                                          {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0},
                                          {false, true}, {7.0, 8.0});
    val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Events)
{
#ifdef HAVE_MPI
    Ewoms::Events val1(Ewoms::DynamicVector<uint64_t>({1,2,3,4,5}));
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(MLimits)
{
#ifdef HAVE_MPI
    Ewoms::MLimits val1{1,2,3,4,5,6,7,8,9,10,11,12};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(MessageLimits)
{
#ifdef HAVE_MPI
    std::vector<Ewoms::MLimits> limits{Ewoms::MLimits{1,2,3,4,5,6,7,8,9,10,11,12}};
    Ewoms::MessageLimits val1(Ewoms::DynamicState<Ewoms::MLimits>(limits,2));
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(VFPInjTable)
{
#ifdef HAVE_MPI
    Ewoms::VFPInjTable val1 = getVFPInjTable();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(VFPProdTable)
{
#ifdef HAVE_MPI
    Ewoms::VFPProdTable val1 = getVFPProdTable();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WTESTWell)
{
#ifdef HAVE_MPI
    Ewoms::WellTestConfig::WTESTWell val1{"test", Ewoms::WellTestConfig::ECONOMIC,
                                         1.0, 2, 3.0, 4};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WellTestConfig)
{
#ifdef HAVE_MPI
    Ewoms::WellTestConfig::WTESTWell tw{"test", Ewoms::WellTestConfig::ECONOMIC,
                                         1.0, 2, 3.0, 4};
    Ewoms::WellTestConfig val1({tw, tw, tw});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WellPolymerProperties)
{
#ifdef HAVE_MPI
    Ewoms::WellPolymerProperties val1{1.0, 2.0, 3, 4, 5};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WellFoamProperties)
{
#ifdef HAVE_MPI
    Ewoms::WellFoamProperties val1{1.0};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WellTracerProperties)
{
#ifdef HAVE_MPI
    Ewoms::WellTracerProperties val1({{"test", 1.0}, {"test2", 2.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UDAValue)
{
#ifdef HAVE_MPI
    Ewoms::UDAValue val1("test");
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
    val1 = Ewoms::UDAValue(1.0);
    auto val22 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val22) == std::get<2>(val22));
    BOOST_CHECK(val1 == std::get<0>(val22));
#endif
}

BOOST_AUTO_TEST_CASE(Connection)
{
#ifdef HAVE_MPI
    Ewoms::Connection val1(Ewoms::Connection::Direction::Y,
                         1.0, Ewoms::Connection::State::SHUT,
                         2, 3, 4.0, 5.0, 6.0, 7.0, 8.0,
                         {9, 10, 11}, Ewoms::Connection::CTFKind::Defaulted,
                         12, 13.0, 14.0, true,
                         15, 16, 17.0);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
                                            Ewoms::Well::InjectorType::OIL,
                                            Ewoms::Well::InjectorCMode::BHP);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WellGuideRate)
{
#ifdef HAVE_MPI
    Ewoms::Well::WellGuideRate val1{true, 1.0, Ewoms::Well::GuideRateTarget::COMB, 2.0};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WellConnections)
{
#ifdef HAVE_MPI
    Ewoms::Connection conn(Ewoms::Connection::Direction::Y,
                         1.0, Ewoms::Connection::State::SHUT,
                         2, 3, 4.0, 5.0, 6.0, 7.0, 8.0,
                         {9, 10, 11}, Ewoms::Connection::CTFKind::Defaulted,
                         12, 13.0, 14.0, true,
                         15, 16, 17.0);
    Ewoms::WellConnections val1(1, 2, 3, {conn, conn});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(SpiralICD)
{
#ifdef HAVE_MPI
    Ewoms::SpiralICD val1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8, 9.0,
                        Ewoms::SpiralICD::Status::OPEN, 10.0);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Valve)
{
#ifdef HAVE_MPI
    Ewoms::Valve val1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, Ewoms::Valve::Status::OPEN);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Segment)
{
#ifdef HAVE_MPI
    Ewoms::Segment val1(1, 2, 3, {1, 2}, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, false,
                      Ewoms::Segment::SegmentType::SICD,
                      std::make_shared<Ewoms::SpiralICD>(),
                      std::make_shared<Ewoms::Valve>());
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Dimension)
{
#ifdef HAVE_MPI
    Ewoms::Dimension val1("test", 1.0, 2.0);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UnitSystem)
{
#ifdef HAVE_MPI
    Ewoms::UnitSystem val1(Ewoms::UnitSystem::UnitType::UNIT_TYPE_METRIC);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WellSegments)
{
#ifdef HAVE_MPI
    Ewoms::Segment seg(1, 2, 3, {1, 2}, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, false,
                     Ewoms::Segment::SegmentType::SICD,
                     std::make_shared<Ewoms::SpiralICD>(),
                     std::make_shared<Ewoms::Valve>());
    Ewoms::WellSegments val1("test", 1.0, 2.0, 3.0,
                           Ewoms::WellSegments::LengthDepth::ABS,
                           Ewoms::WellSegments::CompPressureDrop::HF_,
                           Ewoms::WellSegments::MultiPhaseModel::DF,
                           {seg, seg}, {{1,2},{3,4}});

    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Well)
{
#ifdef HAVE_MPI
    Ewoms::Well val1 = getFullWell();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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

    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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

    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Group)
{
#ifdef HAVE_MPI
    Ewoms::UnitSystem unitSystem;
    Ewoms::Group val1("test1", 1, 2, 3.0, unitSystem,
                    Ewoms::Group::GroupType::PRODUCTION,
                    4.0, true, 5, "test2",
                    Ewoms::IOrderSet<std::string>({"test3", "test4"}, {"test3","test4"}),
                    Ewoms::IOrderSet<std::string>({"test5", "test6"}, {"test5","test6"}),
                    Ewoms::Group::GroupInjectionProperties(),
                    Ewoms::Group::GroupProductionProperties());

    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WList)
{
#ifdef HAVE_MPI
    Ewoms::WList val1({"test1", "test2", "test3"});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WListManager)
{
#ifdef HAVE_MPI
    Ewoms::WList wl({"test1", "test2", "test3"});
    std::map<std::string,Ewoms::WList> data{{"test", wl}, {"test2", wl}};
    Ewoms::WListManager val1(data);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UDQFunction)
{
#ifdef HAVE_MPI
    Ewoms::UDQFunction val1("test", Ewoms::UDQTokenType::binary_op_add);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UDQFunctionTable)
{
#ifdef HAVE_MPI
    Ewoms::UDQFunctionTable::FunctionMap map{{"test",
                                            std::make_shared<Ewoms::UDQFunction>()}};
    Ewoms::UDQFunctionTable val1(Ewoms::UDQParams(true, 1, 2.0, 3.0, 4.0), map);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UDQAssign)
{
#ifdef HAVE_MPI
    Ewoms::UDQAssign val1("test", Ewoms::UDQVarType::NONE,
                        {Ewoms::UDQAssign::AssignRecord{{"test1"}, 1.0},
                         Ewoms::UDQAssign::AssignRecord{{"test2"}, 2.0}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UDQIndex)
{
#ifdef HAVE_MPI
    Ewoms::UDQIndex val1(1, 2, Ewoms::UDQAction::ASSIGN, Ewoms::UDQVarType::WELL_VAR);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UDQConfig)
{
#ifdef HAVE_MPI
    Ewoms::UDQConfig val1 = getUDQConfig();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UDQActiveInputRecord)
{
#ifdef HAVE_MPI
    Ewoms::UDQActive::InputRecord val1(1, "test1", "test2",
                                     Ewoms::UDAControl::WCONPROD_ORAT);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(UDQActiveRecord)
{
#ifdef HAVE_MPI
    Ewoms::UDQActive::Record val1("test1", 1, 2, "test2",
                                Ewoms::UDAControl::WCONPROD_ORAT);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(GuideRateModel)
{
#ifdef HAVE_MPI
    Ewoms::GuideRateModel val1 = getGuideRateModel();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(GuideRateConfigGroup)
{
#ifdef HAVE_MPI
    Ewoms::GuideRateConfig::GroupTarget val1 = getGuideRateConfigGroup();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(GuideRateConfigWell)
{
#ifdef HAVE_MPI
    Ewoms::GuideRateConfig::WellTarget val1 = getGuideRateConfigWell();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(GuideRateConfig)
{
#ifdef HAVE_MPI
    auto model = std::make_shared<Ewoms::GuideRateModel>(getGuideRateModel());
    Ewoms::GuideRateConfig val1(model,
                              {{"test1", getGuideRateConfigWell()}},
                              {{"test2", getGuideRateConfigGroup()}});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(GConSumpGroup)
{
#ifdef HAVE_MPI
    Ewoms::GConSump::GCONSUMPGroup val1{Ewoms::UDAValue(1.0),
                                      Ewoms::UDAValue(2.0),
                                      "test",
                                      3.0, Ewoms::UnitSystem()};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(DeckItem)
{
#ifdef HAVE_MPI
    Ewoms::DeckItem val1({1.0}, {2}, {"test3"}, {Ewoms::UDAValue(4)},
                       Ewoms::type_tag::string, "test5",
                       {Ewoms::value::status::deck_value},
                       true,
                       {Ewoms::Dimension("DimensionLess", 7.0, 8.0)},
                       {Ewoms::Dimension("Metric", 10.0, 11.0)});

    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(DeckRecord)
{
#ifdef HAVE_MPI
    Ewoms::DeckRecord val1 = getDeckRecord();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Location)
{
#ifdef HAVE_MPI
    Ewoms::Location val1{"test", 1};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(DeckKeyword)
{
#ifdef HAVE_MPI
    Ewoms::DeckKeyword val1("test", {"test",1},
                          {getDeckRecord(), getDeckRecord()}, true, false);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Tuning)
{
#ifdef HAVE_MPI
    Ewoms::Tuning val1 = getTuning();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(ASTNode)
{
#ifdef HAVE_MPI
    Ewoms::Action::ASTNode child(number, FuncType::field, "test3", {"test2"}, 2.0, {});
    Ewoms::Action::ASTNode val1(number, FuncType::field, "test1", {"test2"}, 1.0, {child});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(AST)
{
#ifdef HAVE_MPI
    std::shared_ptr<Ewoms::Action::ASTNode> node;
    node.reset(new Ewoms::Action::ASTNode(number, FuncType::field,
                                        "test1", {"test2"}, 1.0, {}));
    Ewoms::Action::AST val1(node);
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Quantity)
{
#ifdef HAVE_MPI
    Ewoms::Action::Quantity val1;
    val1.quantity = "test1";
    val1.args = {"test2", "test3"};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Condition)
{
#ifdef HAVE_MPI
    Ewoms::Action::Condition val1 = getCondition();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(ActionX)
{
#ifdef HAVE_MPI
    Ewoms::Action::ActionX val1 = getActionX();
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Actions)
{
#ifdef HAVE_MPI
    Ewoms::Action::Actions val1({getActionX()});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(Schedule)
{
#ifdef HAVE_MPI
    Ewoms::UnitSystem unitSystem;
    Ewoms::Schedule::WellMap wells;
    wells.insert({"test", {{std::make_shared<Ewoms::Well>(getFullWell())},1}});
    Ewoms::Schedule::GroupMap groups;
    groups.insert({"test", {{std::make_shared<Ewoms::Group>("test1", 1, 2, 3.0, unitSystem,
                                                          Ewoms::Group::GroupType::PRODUCTION,
                                                          4.0, true, 5, "test2",
                                                          Ewoms::IOrderSet<std::string>({"test3", "test4"}, {"test3","test4"}),
                                                          Ewoms::IOrderSet<std::string>({"test5", "test6"}, {"test5","test6"}),
                                                          Ewoms::Group::GroupInjectionProperties(),
                                                          Ewoms::Group::GroupProductionProperties())},1}});
    using VapType = Ewoms::OilVaporizationProperties::OilVaporization;
    Ewoms::DynamicState<Ewoms::OilVaporizationProperties> oilvap{{Ewoms::OilVaporizationProperties(VapType::VAPPARS,
                                                                                   {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0},
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
                         Ewoms::Actdims(1,2,3,4));
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

    Ewoms::Action::Actions acnts({getActionX()});

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
                       getTuning(),
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
                       {{std::make_shared<Ewoms::Action::Actions>(acnts)}, 1},
                       rftc,
                       {std::vector<int>{1}, 1},
                       {{"test", events}});

    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(BrineDensityTable)
{
#ifdef HAVE_MPI
    Ewoms::BrineDensityTable val1({1.0, 2.0, 3.0});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(SummaryNode)
{
#ifdef HAVE_MPI
    auto val1 = Ewoms::SummaryNode{"test1", Ewoms::SummaryNode::Category::Region,
                                 Ewoms::Location{"test2", 1}}
                                 .parameterType(Ewoms::SummaryNode::Type::Pressure)
                                 .namedEntity("test3")
                                 .number(2)
                                 .isUserDefined(true);

    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(SummaryConfig)
{
#ifdef HAVE_MPI
    auto node = Ewoms::SummaryNode{"test1", Ewoms::SummaryNode::Category::Region,
                                 Ewoms::Location{"test2", 1}}
                                 .parameterType(Ewoms::SummaryNode::Type::Pressure)
                                 .namedEntity("test3")
                                 .number(2)
                                 .isUserDefined(true);
    Ewoms::SummaryConfig val1({node}, {"test1", "test2"}, {"test3", "test4"});

    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(PvtwsaltTable)
{
#ifdef HAVE_MPI
    Ewoms::PvtwsaltTable val1(1.0, 2.0, {3.0, 4.0, 5.0});
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
#endif
}

BOOST_AUTO_TEST_CASE(WellBrineProperties)
{
#ifdef HAVE_MPI
    Ewoms::WellBrineProperties val1{1.0};
    auto val2 = PackUnpack(val1);
    BOOST_CHECK(std::get<1>(val2) == std::get<2>(val2));
    BOOST_CHECK(val1 == std::get<0>(val2));
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
