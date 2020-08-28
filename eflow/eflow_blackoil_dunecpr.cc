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

#include <ewoms/eclsimulators/eflow/main.hh>
#include  <ewoms/eclsimulators/linalg/istlsolverflexible.hh>

namespace Ewoms {
  namespace Properties {

    NEW_TYPE_TAG(EclEFlowProblemSimple, INHERITS_FROM(EclEFlowProblem));

    SET_BOOL_PROP(EclEFlowProblemSimple, MatrixAddWellContributions, true);
    SET_INT_PROP(EclEFlowProblemSimple, LinearSolverVerbosity,0);
    SET_SCALAR_PROP(EclEFlowProblemSimple, LinearSolverReduction, 1e-2);
    SET_INT_PROP(EclEFlowProblemSimple, LinearSolverMaxIter, 100);
    SET_BOOL_PROP(EclEFlowProblemSimple, UseAmg, true);//probably not used
    SET_BOOL_PROP(EclEFlowProblemSimple, UseCpr, true);
    SET_INT_PROP(EclEFlowProblemSimple, CprMaxEllIter, 1);
    SET_INT_PROP(EclEFlowProblemSimple, CprEllSolvetype, 3);
    SET_INT_PROP(EclEFlowProblemSimple, CprReuseSetup, 3);
    SET_INT_PROP(EclEFlowProblemSimple, CprSolverVerbose, 0);
    SET_STRING_PROP(EclEFlowProblemSimple, LinearSolverConfiguration, "ilu0");
    SET_STRING_PROP(EclEFlowProblemSimple, SystemStrategy, "quasiimpes");

    SET_PROP(EclEFlowProblemSimple, FluidSystem)
    {
    private:
      using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
      using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);

    public:
        typedef Ewoms::BlackOilFluidSystem<Scalar> type;
    };
    //NEW_TYPE_TAG(EclEFlowProblem, INHERITS_FROM(BlackOilModel, EclBaseProblem));
    SET_TYPE_PROP(EclEFlowProblemSimple, IntensiveQuantities, Ewoms::BlackOilIntensiveQuantities<TypeTag>);
    //SET_TYPE_PROP(EclEFlowProblemSimple, LinearSolverBackend, Ewoms::ISTLSolver<TypeTag>);
    //SET_TAG_PROP(EclEFlowProblemSimple, LinearSolverSplice, ParallelBiCGStabLinearSolver);
    //SET_TYPE_PROP(EclEFlowProblemSimple, LinearSolverBackend, Ewoms::Linear::ParallelBiCGStabSolverBackend<TypeTag>);//not work
    //SET_TYPE_PROP(EclEFlowProblemSimple, LinearSolverBackend, Ewoms::Linear::SuperLUBackend<TypeTag>)//not work
    //SET_TAG_PROP(EclEFlowProblem, FluidState, Ewoms::BlackOilFluidState);
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
    SET_TYPE_PROP(EclEFlowProblemSimple, LinearSolverBackend, Ewoms::ISTLSolverFlexible<TypeTag>);
#else
    SET_TYPE_PROP(EclEFlowProblemSimple, LinearSolverBackend, Ewoms::ISTLSolverCpr<TypeTag>);
#endif
    SET_BOOL_PROP(EclEFlowProblemSimple, EnableStorageCache, true);
    SET_BOOL_PROP(EclEFlowProblemSimple, EnableIntensiveQuantityCache, true);

    //SET_INT_PROP(EclEFlowProblemSimple, NumWellAdjoint, 1);
    //SET_BOOL_PROP(EclEFlowProblem, EnableStorageCache, true);
    //SET_BOOL_PROP(EclEFlowProblem, EnableIntensiveQuantityCache, true);

END_PROPERTIES

int main(int argc, char** argv)
{
    using TypeTag = TTAG(EclEFlowProblemSimple);
    auto mainObject = Ewoms::EFlowNihMain(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
