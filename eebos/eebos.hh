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
/*!
 * \file
 *
 * \brief The common settings for all eebos variants.
 */
#ifndef EEBOS_HH
#define EEBOS_HH

#include "eclproblem.hh"
#include "eclbicgstabbackend.hh"

#include <ewoms/eclsimulators/wells/blackoilwellmodel.hh>
#include <ewoms/eclsimulators/aquifers/blackoilaquifermodel.hh>

#include <ewoms/numerics/utils/start.hh>

namespace Ewoms {
template <class TypeTag>
class EebosProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EebosTypeTag, INHERITS_FROM(BlackOilModel, EclBaseProblem, EFlowModelParameters));

// Set the problem class
SET_TYPE_PROP(EebosTypeTag, Problem, Ewoms::EebosProblem<TypeTag>);

// Enable the more experimental features for eebos
SET_BOOL_PROP(EebosTypeTag, EnableExperiments, true);

// use eflow's well model for now
SET_TYPE_PROP(EebosTypeTag, EclWellModel, Ewoms::BlackoilWellModel<TypeTag>);

// currently, eebos uses the non-multisegment well model by default to avoid
// regressions. the --use-multisegment-well=true|false command line parameter is still
// available in eebos, but hidden from view.
SET_BOOL_PROP(EebosTypeTag, UseMultisegmentWell, false);

// set some properties that are only required by the well model
SET_BOOL_PROP(EebosTypeTag, MatrixAddWellContributions, true);

SET_BOOL_PROP(EebosTypeTag, EnableTerminalOutput, false);

// eflow's well model only works with surface volumes
SET_BOOL_PROP(EebosTypeTag, BlackoilConserveSurfaceVolume, true);

// the values for the residual are for the whole cell instead of for a cubic meter of the cell
SET_BOOL_PROP(EebosTypeTag, UseVolumetricResidual, false);

// by default use eflow's aquifer model for now
SET_TYPE_PROP(EebosTypeTag, EclAquiferModel, Ewoms::BlackoilAquiferModel<TypeTag>);

// use the eebos optimized linear solver backend
SET_TAG_PROP(EebosTypeTag, LinearSolverSplice, EclBiCGStabSolverBackend);

// reducing the residual in the linear solver by a factor of 20 is enough for us
SET_SCALAR_PROP(EebosTypeTag, LinearSolverTolerance, 0.05);

// the default for the allowed volumetric error for oil per second. note that the "main"
// convergence criterium usually is the sum tolerance (specified below)
SET_SCALAR_PROP(EebosTypeTag, NewtonTolerance, 0.1);

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
SET_SCALAR_PROP(EebosTypeTag, EclNewtonRelaxedVolumeFraction, 0.05);

// the maximum volumetric error of a cell in the relaxed region
SET_SCALAR_PROP(EebosTypeTag, EclNewtonRelaxedTolerance, 1e10);

// the tolerated amount of "incorrect" amount of oil per second for the complete
// reservoir. this is scaled by the pore volume of the reservoir, i.e., larger reservoirs
// will tolerate larger residuals.
SET_SCALAR_PROP(EebosTypeTag, EclNewtonSumTolerance, 1.0);

// make all Newton iterations strict, i.e., the volumetric Newton tolerance must be
// always be upheld in the majority of the spatial domain. In this context, "majority"
// means 1 - EclNewtonRelaxedVolumeFraction.
SET_INT_PROP(EebosTypeTag, EclNewtonStrictIterations, 100);

// set the maximum number of Newton iterations to 8 so that we fail quickly (albeit
// relatively often)
SET_INT_PROP(EebosTypeTag, NewtonMaxIterations, 8);

// By default, eebos accepts the result of the time integration unconditionally if the
// smallest time step size is reached.
SET_BOOL_PROP(EebosTypeTag, ContinueOnConvergenceError, true);

END_PROPERTIES

namespace Ewoms {
template <class TypeTag>
class EebosProblem : public EclProblem<TypeTag>
{
    typedef EclProblem<TypeTag> ParentType;

public:
    static void registerParameters()
    {
        ParentType::registerParameters();

        Ewoms::BlackoilModelParameters<TypeTag>::registerParameters();
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTerminalOutput, "Do *NOT* use!");
        EWOMS_HIDE_PARAM(TypeTag, DbhpMaxRel);
        EWOMS_HIDE_PARAM(TypeTag, DwellFractionMax);
        EWOMS_HIDE_PARAM(TypeTag, MaxResidualAllowed);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceMb);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceCnv);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceCnvRelaxed);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceWells);
        EWOMS_HIDE_PARAM(TypeTag, ToleranceWellControl);
        EWOMS_HIDE_PARAM(TypeTag, MaxWelleqIter);
        EWOMS_HIDE_PARAM(TypeTag, UseMultisegmentWell);
        EWOMS_HIDE_PARAM(TypeTag, TolerancePressureMsWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxPressureChangeMsWells);
        EWOMS_HIDE_PARAM(TypeTag, UseInnerIterationsMsWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxInnerIterMsWells);
        EWOMS_HIDE_PARAM(TypeTag, UseInnerIterationsWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxInnerIterWells);
        EWOMS_HIDE_PARAM(TypeTag, MaxSinglePrecisionDays);
        EWOMS_HIDE_PARAM(TypeTag, MaxStrictIter);
        EWOMS_HIDE_PARAM(TypeTag, SolveWelleqInitially);
        EWOMS_HIDE_PARAM(TypeTag, UpdateEquationsScaling);
        EWOMS_HIDE_PARAM(TypeTag, UseUpdateStabilization);
        EWOMS_HIDE_PARAM(TypeTag, MatrixAddWellContributions);
        EWOMS_HIDE_PARAM(TypeTag, EnableTerminalOutput);
    }

    // inherit the constructors
    using ParentType::EclProblem;
};
}

#endif // EEBOS_HH
