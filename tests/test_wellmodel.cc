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

#define BOOST_TEST_MODULE WellModelTest

#include <chrono>

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>

#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablemanager.hh>

#include <ewoms/eclgrids/gridmanager.hh>
#include <ewoms/eclio/parser/units/units.hh>

#include <ewoms/material/fluidmatrixinteractions/eclmateriallawmanager.hh>
#include <ewoms/eclgrids/gridhelpers.hh>
#include <ewoms/eclsimulators/eflow/eflowmain.hh>
#include <ewoms/eclsimulators/eflow/blackoilmodel.hh>

#include <eebos/eclproblem.hh>
#include <ewoms/numerics/utils/start.hh>

#include <ewoms/eclsimulators/wells/standardwell.hh>
#include <ewoms/eclsimulators/wells/blackoilwellmodel.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

using StandardWell = Ewoms::StandardWell<TTAG(EclEFlowProblem)>;

struct SetupTest {

    using Grid = UnstructuredGrid;

    SetupTest ()
    {
        Ewoms::Parser parser;
        auto deck = parser.parseFile("TESTWELLMODEL.DATA");
        ecl_state.reset(new Ewoms::EclipseState(deck) );
        {
          const Ewoms::TableManager table ( deck );
          const Ewoms::Eclipse3DProperties eclipseProperties ( deck , table, ecl_state->getInputGrid());
          const Ewoms::Runspec runspec (deck);
          schedule.reset( new Ewoms::Schedule(deck, *ecl_state));
          summaryState.reset( new Ewoms::SummaryState(std::chrono::system_clock::from_time_t(schedule->getStartTime())));
        }
        current_timestep = 0;
    };

    std::unique_ptr<const Ewoms::EclipseState> ecl_state;
    std::unique_ptr<const Ewoms::Schedule> schedule;
    std::unique_ptr<Ewoms::SummaryState> summaryState;
    std::vector<std::vector<Ewoms::PerforationData>> well_perf_data;
    int current_timestep;
};

struct GlobalFixture {
    GlobalFixture()
    {
        int argcDummy = 1;
        const char *tmp[] = {"test_wellmodel"};
        char **argvDummy = const_cast<char**>(tmp);

        // MPI setup.
#if HAVE_DUNE_FEM
        Dune::Fem::MPIManager::initialize(argcDummy, argvDummy);
#else
        Dune::MPIHelper::instance(argcDummy, argvDummy);
#endif

        Ewoms::EFlowMain<TTAG(EclEFlowProblem)>::setupParameters_(argcDummy, argvDummy);
    }
};

BOOST_GLOBAL_FIXTURE(GlobalFixture);

BOOST_AUTO_TEST_CASE(TestStandardWellInput) {
    const SetupTest setup_test;
    const auto& wells_ecl = setup_test.schedule->getWells(setup_test.current_timestep);
    BOOST_CHECK_EQUAL( wells_ecl.size(), 2);
    const Ewoms::Well& well = wells_ecl[1];
    const Ewoms::BlackoilModelParameters<TTAG(EclEFlowProblem) > param;

    // For the conversion between the surface volume rate and resrevoir voidage rate
    typedef Ewoms::BlackOilFluidSystem<double> FluidSystem;
    using RateConverterType = Ewoms::RateConverter::
        SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;
    // Compute reservoir volumes for RESV controls.
    Ewoms::PhaseUsage phaseUsage;
    std::unique_ptr<RateConverterType> rateConverter;
    // Compute reservoir volumes for RESV controls.
    rateConverter.reset(new RateConverterType (phaseUsage,
                                     std::vector<int>(10, 0)));

    Ewoms::PerforationData dummy;
    std::vector<Ewoms::PerforationData> pdata(well.getConnections().size(), dummy);

    BOOST_CHECK_THROW( StandardWell( well, -1, param, *rateConverter, 0, 3, 3, 0, 0, pdata), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestBehavoir) {
    const SetupTest setup_test;
    const auto& wells_ecl = setup_test.schedule->getWells(setup_test.current_timestep);
    const int current_timestep = setup_test.current_timestep;
    std::vector<std::unique_ptr<const StandardWell> >  wells;

    {
        const int nw = wells_ecl.size();
        const Ewoms::BlackoilModelParameters<TTAG(EclEFlowProblem)> param;

        for (int w = 0; w < nw; ++w) {
            // For the conversion between the surface volume rate and resrevoir voidage rate
            typedef Ewoms::BlackOilFluidSystem<double> FluidSystem;
            using RateConverterType = Ewoms::RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;
            // Compute reservoir volumes for RESV controls.
            // TODO: not sure why for this class the initlizer list does not work
            // otherwise we should make a meaningful const PhaseUsage here.
            Ewoms::PhaseUsage phaseUsage;
            std::unique_ptr<RateConverterType> rateConverter;
            // Compute reservoir volumes for RESV controls.
            rateConverter.reset(new RateConverterType (phaseUsage,
                                             std::vector<int>(10, 0)));
            Ewoms::PerforationData dummy;
            std::vector<Ewoms::PerforationData> pdata(wells_ecl[w].getConnections().size(), dummy);
            wells.emplace_back(new StandardWell(wells_ecl[w], current_timestep, param, *rateConverter, 0, 3, 3, w, 0, pdata) );
        }
    }

    // first well, it is a production well from the deck
    {
        const auto& well = wells[0];
        BOOST_CHECK_EQUAL(well->name(), "PROD1");
        BOOST_CHECK(well->isProducer());
        BOOST_CHECK(well->numEq == 3);
        BOOST_CHECK(well->numStaticWellEq== 4);
    }

    // second well, it is the injection well from the deck
    {
        const auto& well = wells[1];
        BOOST_CHECK_EQUAL(well->name(), "INJE1");
        BOOST_CHECK(well->isInjector());
        BOOST_CHECK(well->numEq == 3);
        BOOST_CHECK(well->numStaticWellEq== 4);
    }
}
