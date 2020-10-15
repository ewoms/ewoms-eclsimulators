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

#include <eflow/eflow_oilwater_polymer_injectivity.hh>

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

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclEFlowOilWaterPolymerInjectivityProblem, INHERITS_FROM(EclEFlowProblem));
SET_BOOL_PROP(EclEFlowOilWaterPolymerInjectivityProblem, EnablePolymer, true);
SET_BOOL_PROP(EclEFlowOilWaterPolymerInjectivityProblem, EnablePolymerMW, true);
//! The indices required by the model
// For this case, there will be two primary variables introduced for the polymer
// polymer concentration and polymer molecular weight
SET_PROP(EclEFlowOilWaterPolymerInjectivityProblem, Indices)
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using BaseTypeTag = TTAG(EclEFlowProblem);
    using FluidSystem = GET_PROP_TYPE(BaseTypeTag, FluidSystem);

public:
    typedef Ewoms::BlackOilTwoPhaseIndices<0,
                                         2,
                                         0,
                                         GET_PROP_VALUE(TypeTag, EnableFoam),
                                         GET_PROP_VALUE(TypeTag, EnableBrine),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::gasCompIdx> type;
};

END_PROPERTIES

namespace Ewoms {
/* void eflowOilWaterPolymerInjectivitySetDeck(Deck& deck, EclipseState& eclState)
{
    using TypeTag = TTAG(EclEFlowOilWaterPolymerInjectivityProblem);
    using Vanguard = GET_PROP_TYPE(TypeTag, Vanguard);

    Vanguard::setExternalDeck(&deck, &eclState);
} */

// ----------------- Main program -----------------
int eflowOilWaterPolymerInjectivityMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    Ewoms::resetLocale();

#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    Ewoms::EFlowMain<TTAG(EclEFlowOilWaterPolymerInjectivityProblem)>
        mainfunc {argc, argv, outputCout, outputFiles};
    return mainfunc.execute();
}

}
