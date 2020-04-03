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
#include "eflow/eflow_tag.hh"

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 6)
#include  <ewoms/eclsimulators/linalg/istlsolverflexible.hh>
#else
#include  <ewoms/eclsimulators/linalg/istlsolvercpr.hh>
#endif

BEGIN_PROPERTIES
NEW_TYPE_TAG(EclEFlowProblemSimple, INHERITS_FROM(EclEFlowProblem));
NEW_PROP_TAG(FluidState);
//SET_TYPE_PROP(EclBaseProblem, Problem, Ewoms::EclProblem<TypeTag>);
SET_PROP(EclEFlowProblemSimple, FluidState)
    {
    private:
      typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
      typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
      enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
      enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
      enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
      enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
      static const bool compositionSwitchEnabled = Indices::gasEnabled;

    public:
//typedef Ewoms::BlackOilFluidSystemSimple<Scalar> type;
       typedef Ewoms::BlackOilFluidState<Evaluation, FluidSystem, enableTemperature, enableEnergy, compositionSwitchEnabled,  Indices::numPhases > type;
};

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
END_PROPERTIES

namespace Ewoms {
  namespace Properties {

    SET_PROP(EclEFlowProblemSimple, FluidSystem)
    {
    private:
      //typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
      typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;

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
  }
}

int main(int argc, char** argv)
{
    typedef TTAG(EclEFlowProblemSimple) TypeTag;
    return mainEFlow<TypeTag>(argc, argv);
}
