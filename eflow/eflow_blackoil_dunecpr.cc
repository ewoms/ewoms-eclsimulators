// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#include <dune/common/version.hh>

#if !DUNE_VERSION_NEWER(DUNE_GRID, 2,6)

# warning "The CPR variant of eflow requires Dune 2.6 or newer"

#include <cstdlib>

int main(int, char**)
{
  return EXIT_FAILURE;
}

#else

#include <ewoms/eclsimulators/eflow/main.hh>
#include  <ewoms/eclsimulators/linalg/istlsolverflexible.hh>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclEFlowProblemSimple, INHERITS_FROM(EclEFlowProblem));

SET_BOOL_PROP(EclEFlowProblemSimple, MatrixAddWellContributions, true);
SET_INT_PROP(EclEFlowProblemSimple, LinearSolverVerbosity, 0);
SET_SCALAR_PROP(EclEFlowProblemSimple, LinearSolverReduction, 1e-2);
SET_INT_PROP(EclEFlowProblemSimple, LinearSolverMaxIter, 100);
SET_INT_PROP(EclEFlowProblemSimple, CprMaxEllIter, 1);
SET_INT_PROP(EclEFlowProblemSimple, CprEllSolvetype, 3);
SET_INT_PROP(EclEFlowProblemSimple, CprReuseSetup, 3);
SET_STRING_PROP(EclEFlowProblemSimple, Linsolver, "ilu0");

SET_PROP(EclEFlowProblemSimple, FluidSystem)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);

public:
    typedef Ewoms::BlackOilFluidSystem<Scalar> type;
};

SET_TYPE_PROP(EclEFlowProblemSimple,
              IntensiveQuantities,
              Ewoms::BlackOilIntensiveQuantities<TypeTag>);

SET_TYPE_PROP(EclEFlowProblemSimple, LinearSolverBackend, Ewoms::ISTLSolverFlexible<TypeTag>);
SET_BOOL_PROP(EclEFlowProblemSimple, EnableStorageCache, true);
SET_BOOL_PROP(EclEFlowProblemSimple, EnableIntensiveQuantityCache, true);

END_PROPERTIES

namespace Ewoms {
namespace CO2DefaultTables {
#include <ewoms/material/components/co2tables.inc.cc>
}}

int main(int argc, char** argv)
{
    using TypeTag = TTAG(EclEFlowProblemSimple);
    auto mainObject = Ewoms::EFlowNihMain(argc, argv);
    return mainObject.runStatic<TypeTag>();
}

#endif
