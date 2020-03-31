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

#define BOOST_TEST_MODULE WellStateFIBOTest

#include <ewoms/eclsimulators/wells/wellstatefullyimplicitblackoil.hh>

#include <boost/test/unit_test.hpp>

#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/parsecontext.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/units/units.hh>

#include <ewoms/eclgrids/gridhelpers.hh>

#include <ewoms/eclsimulators/deprecated/props/blackoilphases.hh>
#include <ewoms/eclsimulators/deprecated/props/phaseusagefromdeck.hh>

#include <ewoms/eclgrids/gridmanager.hh>

#include <chrono>
#include <cstddef>
#include <string>

struct Setup
{
    Setup(const std::string& filename)
        : Setup(Ewoms::Parser{}.parseFile(filename))
    {}

    Setup(const Ewoms::Deck& deck)
        : es   (deck)
        , pu   (Ewoms::phaseUsageFromDeck(es))
        , grid (es.getInputGrid())
        , sched(deck, es)
        , st(std::chrono::system_clock::from_time_t(sched.getStartTime()))
    {
        initWellPerfData();
    }

    void initWellPerfData()
    {
        const auto& wells = sched.getWells(0);
        const auto& cartDims = Ewoms::UgGridHelpers::cartDims(*grid.c_grid());
        const int* compressed_to_cartesian = Ewoms::UgGridHelpers::globalCell(*grid.c_grid());
        std::vector<int> cartesian_to_compressed(cartDims[0] * cartDims[1] * cartDims[2], -1);
        for (int ii = 0; ii < Ewoms::UgGridHelpers::numCells(*grid.c_grid()); ++ii) {
            cartesian_to_compressed[compressed_to_cartesian[ii]] = ii;
        }
        well_perf_data.resize(wells.size());
        int well_index = 0;
        for (const auto& well : wells) {
            well_perf_data[well_index].clear();
            well_perf_data[well_index].reserve(well.getConnections().size());
            for (const auto& completion : well.getConnections()) {
                if (completion.state() == Ewoms::Connection::State::OPEN) {
                    const int i = completion.getI();
                    const int j = completion.getJ();
                    const int k = completion.getK();
                    const int cart_grid_indx = i + cartDims[0] * (j + cartDims[1] * k);
                    const int active_index = cartesian_to_compressed[cart_grid_indx];
                    if (active_index < 0) {
                        const std::string msg
                            = ("Cell with i,j,k indices " + std::to_string(i) + " " + std::to_string(j) + " "
                               + std::to_string(k) + " not found in grid (well = " + well.name() + ").");
                        EWOMS_THROW(std::runtime_error, msg);
                    } else {
                        Ewoms::PerforationData pd;
                        pd.cell_index = active_index;
                        pd.connection_transmissibility_factor = completion.CF();
                        pd.satnum_id = completion.satTableId();
                        well_perf_data[well_index].push_back(pd);
                    }
                } else {
                    if (completion.state() != Ewoms::Connection::State::SHUT) {
                        EWOMS_THROW(std::runtime_error,
                                  "Completion state: " << Ewoms::Connection::State2String(completion.state()) << " not handled");
                    }
                }
            }
            ++well_index;
        }
    }

    Ewoms::EclipseState es;
    Ewoms::PhaseUsage   pu;
    Ewoms::GridManager  grid;
    Ewoms::Schedule     sched;
    Ewoms::SummaryState st;
    std::vector<std::vector<Ewoms::PerforationData>> well_perf_data;
};

namespace {
    Ewoms::WellStateFullyImplicitBlackoil
    buildWellState(const Setup& setup, const std::size_t timeStep)
    {
        auto state  = Ewoms::WellStateFullyImplicitBlackoil{};

        const auto cpress =
            std::vector<double>(setup.grid.c_grid()->number_of_cells,
                                100.0*Ewoms::unit::barsa);

        state.init(cpress, setup.sched,
                   setup.sched.getWells(timeStep),
                   timeStep, nullptr, setup.pu, setup.well_perf_data, setup.st, setup.sched.getWells(timeStep).size());

        state.initWellStateMSWell(setup.sched.getWells(timeStep),
                                  setup.pu, nullptr);

        return state;
    }

    void setSegPress(const std::vector<Ewoms::Well>& wells,
                     Ewoms::WellStateFullyImplicitBlackoil& wstate)
    {
        const auto nWell = wells.size();

        auto& segPress = wstate.segPress();

        for (auto wellID = 0*nWell; wellID < nWell; ++wellID) {
            const auto& well     = wells[wellID];
            const auto  topSegIx = wstate.topSegmentIndex(wellID);
            const auto  pressTop = 100.0 * wellID;

            auto* press = &segPress[topSegIx];

            press[0] = pressTop;

            if (! well.isMultiSegment()) {
                continue;
            }

            const auto& segSet = well.getSegments();
            const auto  nSeg   = segSet.size();

            for (auto segID = 0*nSeg + 1; segID < nSeg; ++segID) {
                // One-based numbering scheme for segments.
                const auto segNo = segSet[segID].segmentNumber();
                press[segNo - 1] = pressTop + 1.0*(segNo - 1);
            }
        }
    }

  void setSegRates(const std::vector<Ewoms::Well>& wells,
                     const Ewoms::PhaseUsage&               pu,
                     Ewoms::WellStateFullyImplicitBlackoil& wstate)
    {
        const auto wat = pu.phase_used[Ewoms::BlackoilPhases::Aqua];
        const auto iw  = wat ? pu.phase_pos[Ewoms::BlackoilPhases::Aqua] : -1;

        const auto oil = pu.phase_used[Ewoms::BlackoilPhases::Liquid];
        const auto io  = oil ? pu.phase_pos[Ewoms::BlackoilPhases::Liquid] : -1;

        const auto gas = pu.phase_used[Ewoms::BlackoilPhases::Vapour];
        const auto ig  = gas ? pu.phase_pos[Ewoms::BlackoilPhases::Vapour] : -1;

        const auto np = wstate.numPhases();

        const auto nWell = wells.size();

        auto& segRates = wstate.segRates();

        for (auto wellID = 0*nWell; wellID < nWell; ++wellID) {
            const auto& well     = wells[wellID];
            const auto  topSegIx = wstate.topSegmentIndex(wellID);
            const auto  rateTop  = 1000.0 * wellID;

            if (wat) { segRates[np*topSegIx + iw] = rateTop; }
            if (oil) { segRates[np*topSegIx + io] = rateTop; }
            if (gas) { segRates[np*topSegIx + ig] = rateTop; }

            if (! well.isMultiSegment()) {
                continue;
            }

            const auto& segSet = well.getSegments();
            const auto  nSeg   = segSet.size();

            for (auto segID = 0*nSeg + 1; segID < nSeg; ++segID) {
                // One-based numbering scheme for segments.
                const auto segNo = segSet[segID].segmentNumber();

                auto* rates = &segRates[(topSegIx + segNo - 1) * np];

                if (wat) { rates[iw] = rateTop + 100.0*(segNo - 1); }
                if (oil) { rates[io] = rateTop + 200.0*(segNo - 1); }
                if (gas) { rates[ig] = rateTop + 400.0*(segNo - 1); }
            }
        }
    }
} // Anonymous

BOOST_AUTO_TEST_SUITE(Segment)

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Linearisation)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    const auto wstate = buildWellState(setup, tstep);

    BOOST_CHECK_EQUAL(wstate.numSegment(), 6 + 1);

    const auto& wells = setup.sched.getWellsatEnd();
    BOOST_CHECK_EQUAL(wells.size(), 2);

    const auto prod01_first = wells[0].name() == "PROD01";

    BOOST_CHECK_EQUAL(wstate.topSegmentIndex(0), 0);
    BOOST_CHECK_EQUAL(wstate.topSegmentIndex(1),
                      prod01_first ? 6 : 1);
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Pressure)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    auto wstate = buildWellState(setup, tstep);

    const auto& wells = setup.sched.getWells(tstep);
    const auto prod01_first = wells[0].name() == "PROD01";

    setSegPress(wells, wstate);

    const auto rpt = wstate.report(setup.pu, setup.grid.c_grid()->global_cell);

    {
        const auto& xw = rpt.at("INJE01");

        BOOST_CHECK_EQUAL(xw.segments.size(), 1); // Top Segment

        const auto& xseg = xw.segments.at(1);

        BOOST_CHECK_EQUAL(xseg.segNumber, 1);
        BOOST_CHECK_CLOSE(xseg.pressure, prod01_first ? 100.0 : 0.0, 1.0e-10);
    }

    {
        const auto expect_nSeg = 6;
        const auto& xw = rpt.at("PROD01");

        BOOST_CHECK_EQUAL(xw.segments.size(), expect_nSeg);

        const auto pressTop = prod01_first ? 0.0 : 100.0;

        for (auto segID = 0; segID < expect_nSeg; ++segID) {
            const auto& xseg = xw.segments.at(segID + 1);

            BOOST_CHECK_EQUAL(xseg.segNumber, segID + 1);
            BOOST_CHECK_CLOSE(xseg.pressure, pressTop + 1.0*segID, 1.0e-10);
        }
    }
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Rates)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    auto wstate = buildWellState(setup, tstep);

    const auto wells = setup.sched.getWells(tstep);
    const auto prod01_first = wells[0].name() == "PROD01";

    const auto& pu = setup.pu;

    setSegRates(wells, pu, wstate);

    const auto rpt = wstate.report(pu, setup.grid.c_grid()->global_cell);

    const auto wat = pu.phase_used[Ewoms::BlackoilPhases::Aqua];
    const auto oil = pu.phase_used[Ewoms::BlackoilPhases::Liquid];
    const auto gas = pu.phase_used[Ewoms::BlackoilPhases::Vapour];

    BOOST_CHECK(wat && oil && gas);

    {
        const auto rateTop = prod01_first ? 1000.0 : 0.0;

        const auto& xw = rpt.at("INJE01");

        BOOST_CHECK_EQUAL(xw.segments.size(), 1); // Top Segment

        const auto& xseg = xw.segments.at(1);

        BOOST_CHECK_EQUAL(xseg.segNumber, 1);
        BOOST_CHECK_CLOSE(xseg.rates.get(Ewoms::data::Rates::opt::wat),
                          rateTop, 1.0e-10);

        BOOST_CHECK_CLOSE(xseg.rates.get(Ewoms::data::Rates::opt::oil),
                          rateTop, 1.0e-10);

        BOOST_CHECK_CLOSE(xseg.rates.get(Ewoms::data::Rates::opt::gas),
                          rateTop, 1.0e-10);
    }

    {
        const auto expect_nSeg = 6;
        const auto& xw = rpt.at("PROD01");

        BOOST_CHECK_EQUAL(xw.segments.size(), expect_nSeg);

        const auto rateTop = prod01_first ? 0.0 : 1000.0;

        for (auto segNum = 1; segNum <= expect_nSeg; ++segNum) {
            const auto& xseg = xw.segments.at(segNum);

            BOOST_CHECK_EQUAL(xseg.segNumber, segNum);

            BOOST_CHECK_CLOSE(xseg.rates.get(Ewoms::data::Rates::opt::wat),
                              rateTop + 100.0*(segNum - 1), 1.0e-10);

            BOOST_CHECK_CLOSE(xseg.rates.get(Ewoms::data::Rates::opt::oil),
                              rateTop + 200.0*(segNum - 1), 1.0e-10);

            BOOST_CHECK_CLOSE(xseg.rates.get(Ewoms::data::Rates::opt::gas),
                              rateTop + 400.0*(segNum - 1), 1.0e-10);
        }
    }
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(STOP_well)
{
    /*
      This test verifies that the perforation pressures is correctly initialized
      also for wells in the STOP state.
    */
    const Setup setup{ "wells_manager_data_wellSTOP.data" };
    auto wstate = buildWellState(setup, 0);
    for (const auto& p : wstate.perfPress())
        BOOST_CHECK(p > 0);
}

BOOST_AUTO_TEST_SUITE_END()
