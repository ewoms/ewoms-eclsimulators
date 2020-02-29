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
#include "config.h"

#if HAVE_DUNE_FEM
#warning "eflow is incompatible with dune-fem. Disabling."
#undef HAVE_DUNE_FEM
#endif // HAVE_DUNE_FEM

// Define making clear that the simulator supports AMG
#define EFLOW_SUPPORT_AMG 1

#include <eflow/eflow_oilwater.hh>

#include <ewoms/common/resetlocale.hh>
#include <ewoms/numerics/models/blackoil/blackoiltwophaseindices.hh>

#include <ewoms/eclgrids/cpgrid.hh>
#include <ewoms/eclsimulators/eflow/simulatorfullyimplicitblackoil.hh>
#include <ewoms/eclsimulators/eflow/eflowmain.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

namespace Ewoms {
namespace Properties {
NEW_TYPE_TAG(EclEFlowOilWaterProblem, INHERITS_FROM(EclEFlowProblem));

//! The indices required by the model
SET_PROP(EclEFlowOilWaterProblem, Indices)
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    typedef TTAG(EclEFlowProblem) BaseTypeTag;
    typedef typename GET_PROP_TYPE(BaseTypeTag, FluidSystem) FluidSystem;

public:
    typedef Ewoms::BlackOilTwoPhaseIndices<GET_PROP_VALUE(TypeTag, EnableSolvent),
                                         GET_PROP_VALUE(TypeTag, EnablePolymer),
                                         GET_PROP_VALUE(TypeTag, EnableEnergy),
                                         GET_PROP_VALUE(TypeTag, EnableFoam),
                                         GET_PROP_VALUE(TypeTag, EnableBrine),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::gasCompIdx> type;
};
}}

namespace Ewoms {
void eflowOilWaterSetDeck(double setupTime, Deck *deck, EclipseState& eclState, Schedule& schedule, SummaryConfig& summaryConfig)
{
    typedef TTAG(EclEFlowOilWaterProblem) TypeTag;
    typedef GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;

    Vanguard::setExternalSetupTime(setupTime);
    Vanguard::setExternalDeck(deck);
    Vanguard::setExternalEclState(&eclState);
    Vanguard::setExternalSchedule(&schedule);
    Vanguard::setExternalSummaryConfig(&summaryConfig);
}

// ----------------- Main program -----------------
int eflowOilWaterMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    Ewoms::resetLocale();

#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    Ewoms::EFlowMain<TTAG(EclEFlowOilWaterProblem)> mainfunc;
    return mainfunc.execute(argc, argv, outputCout, outputFiles);
}

}
