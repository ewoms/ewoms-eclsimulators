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
#ifndef EWOMS_BLACKOILMODELPARAMETERS_EEBOS_HH
#define EWOMS_BLACKOILMODELPARAMETERS_EEBOS_HH

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <string>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EFlowModelParameters);

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(EclDeckFileName);

NEW_PROP_TAG(DbhpMaxRel);
NEW_PROP_TAG(DwellFractionMax);
NEW_PROP_TAG(MaxResidualAllowed);
NEW_PROP_TAG(RelaxedMaxPvFraction);
NEW_PROP_TAG(ToleranceMb);
NEW_PROP_TAG(ToleranceCnv);
NEW_PROP_TAG(ToleranceCnvRelaxed);
NEW_PROP_TAG(ToleranceWells);
NEW_PROP_TAG(ToleranceWellControl);
NEW_PROP_TAG(MaxWelleqIter);
NEW_PROP_TAG(UseMultisegmentWell);
NEW_PROP_TAG(MaxSinglePrecisionDays);
NEW_PROP_TAG(MaxStrictIter);
NEW_PROP_TAG(SolveWelleqInitially);
NEW_PROP_TAG(UpdateEquationsScaling);
NEW_PROP_TAG(UseUpdateStabilization);
NEW_PROP_TAG(MatrixAddWellContributions);
NEW_PROP_TAG(EnableWellOperabilityCheck);

// parameters for multisegment wells
NEW_PROP_TAG(TolerancePressureMsWells);
NEW_PROP_TAG(MaxPressureChangeMsWells);
NEW_PROP_TAG(UseInnerIterationsMsWells);
NEW_PROP_TAG(MaxInnerIterMsWells);
NEW_PROP_TAG(StrictInnerIterMsWells);
NEW_PROP_TAG(RelaxedEFlowTolInnerIterMsw);
NEW_PROP_TAG(RelaxedPressureTolInnerIterMsw);
NEW_PROP_TAG(RegularizationFactorMsw);
NEW_PROP_TAG(UseInnerIterationsWells);
NEW_PROP_TAG(MaxInnerIterWells);
NEW_PROP_TAG(AlternativeWellRateInit);

SET_SCALAR_PROP(EFlowModelParameters, DbhpMaxRel, 1.0);
SET_SCALAR_PROP(EFlowModelParameters, DwellFractionMax, 0.2);
SET_SCALAR_PROP(EFlowModelParameters, MaxResidualAllowed, 1e7);
SET_SCALAR_PROP(EFlowModelParameters, RelaxedMaxPvFraction, 0.03);
SET_SCALAR_PROP(EFlowModelParameters, ToleranceMb, 1e-6);
SET_SCALAR_PROP(EFlowModelParameters, ToleranceCnv, 1e-2);
SET_SCALAR_PROP(EFlowModelParameters, ToleranceCnvRelaxed, 1);
SET_SCALAR_PROP(EFlowModelParameters, ToleranceWells, 1e-4);
SET_SCALAR_PROP(EFlowModelParameters, ToleranceWellControl, 1e-7);
SET_INT_PROP(EFlowModelParameters, MaxWelleqIter, 30);
SET_BOOL_PROP(EFlowModelParameters, UseMultisegmentWell, true);
SET_SCALAR_PROP(EFlowModelParameters, MaxSinglePrecisionDays, 20.0);
SET_INT_PROP(EFlowModelParameters, MaxStrictIter, 0);
SET_BOOL_PROP(EFlowModelParameters, SolveWelleqInitially, true);
SET_BOOL_PROP(EFlowModelParameters, UpdateEquationsScaling, false);
SET_BOOL_PROP(EFlowModelParameters, UseUpdateStabilization, true);
SET_BOOL_PROP(EFlowModelParameters, MatrixAddWellContributions, false);
SET_SCALAR_PROP(EFlowModelParameters, TolerancePressureMsWells, 0.01*1e5);
SET_SCALAR_PROP(EFlowModelParameters, MaxPressureChangeMsWells, 10*1e5);
SET_BOOL_PROP(EFlowModelParameters, UseInnerIterationsMsWells, true);
SET_INT_PROP(EFlowModelParameters, MaxInnerIterMsWells, 100);
SET_BOOL_PROP(EFlowModelParameters, UseInnerIterationsWells, false);
SET_INT_PROP(EFlowModelParameters, MaxInnerIterWells, 50);
SET_BOOL_PROP(EFlowModelParameters, AlternativeWellRateInit, true);
SET_INT_PROP(EFlowModelParameters, StrictInnerIterMsWells, 40);
SET_SCALAR_PROP(EFlowModelParameters, RegularizationFactorMsw, 1);
SET_BOOL_PROP(EFlowModelParameters, EnableWellOperabilityCheck, true);
SET_SCALAR_PROP(EFlowModelParameters, RelaxedEFlowTolInnerIterMsw, 1);
SET_SCALAR_PROP(EFlowModelParameters, RelaxedPressureTolInnerIterMsw, 0.5e5);

END_PROPERTIES

namespace Ewoms
{
    /// Solver parameters for the BlackoilModel.
    template <class TypeTag>
    struct BlackoilModelParameters
    {
    private:
        using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

    public:
        /// Max relative change in bhp in single iteration.
        double dbhp_max_rel_;
        /// Max absolute change in well volume fraction in single iteration.
        double dwell_fraction_max_;
        /// Absolute max limit for residuals.
        double max_residual_allowed_;
        //// Max allowed pore volume faction where CNV is violated. Below the
        //// relaxed tolerance tolerance_cnv_relaxed_ is used.
        double relaxed_max_pv_fraction_;
        /// Relative mass balance tolerance (total mass balance error).
        double tolerance_mb_;
        /// Local convergence tolerance (max of local saturation errors).
        double tolerance_cnv_;
        /// Relaxed local convergence tolerance (used when iter >= max_strict_iter_).
        double tolerance_cnv_relaxed_;
        /// Well convergence tolerance.
        double tolerance_wells_;
        /// Tolerance for the well control equations
        //  TODO: it might need to distinguish between rate control and pressure control later
        double tolerance_well_control_;
        /// Tolerance for the pressure equations for multisegment wells
        double tolerance_pressure_ms_wells_;

        /// Relaxed tolerance for the inner iteration for the MSW flow solution
        double relaxed_inner_tolerance_flow_ms_well_;
        /// Relaxed tolerance for the inner iteration for the MSW pressure solution
        double relaxed_inner_tolerance_pressure_ms_well_;

        /// Maximum pressure change over an iteratio for ms wells
        double max_pressure_change_ms_wells_;

        /// Whether to use inner iterations for ms wells
        bool use_inner_iterations_ms_wells_;

        /// Maximum inner iteration number for ms wells
        int max_inner_iter_ms_wells_;

        /// Strict inner iteration number for ms wells
        int strict_inner_iter_ms_wells_;

        /// Regularization factor for ms wells
        int regularization_factor_ms_wells_;

        /// Whether to use inner iterations for standard wells
        bool use_inner_iterations_wells_;

        /// Maximum inner iteration number for standard wells
        int max_inner_iter_wells_;

        /// Maximum iteration number of the well equation solution
        int max_welleq_iter_;

        /// Tolerance for time step in seconds where single precision can be used
        /// for solving for the Jacobian
        double maxSinglePrecisionTimeStep_;

        /// Maximum number of Newton iterations before we give up on the CNV convergence criterion
        int max_strict_iter_;

        /// Solve well equation initially
        bool solve_welleq_initially_;

        /// Update scaling factors for mass balance equations
        bool update_equations_scaling_;

        /// Try to detect oscillation or stagnation.
        bool use_update_stabilization_;

        /// Whether to use MultisegmentWell to handle multisegment wells
        /// it is something temporary before the multisegment well model is considered to be
        /// well developed and tested.
        /// if it is false, we will handle multisegment wells as standard wells, which will be
        /// the default behavoir for the moment. Later, we might set it to be true by default if necessary
        bool use_multisegment_well_;

        /// The file name of the deck
        std::string deck_file_name_;

        // Whether to add influences of wells between cells to the matrix and preconditioner matrix
        bool matrix_add_well_contributions_;

        /// Construct from user parameters or defaults.
        BlackoilModelParameters()
        {
            dbhp_max_rel_=  EWOMS_GET_PARAM(TypeTag, Scalar, DbhpMaxRel);
            dwell_fraction_max_ = EWOMS_GET_PARAM(TypeTag, Scalar, DwellFractionMax);
            max_residual_allowed_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxResidualAllowed);
            relaxed_max_pv_fraction_ = EWOMS_GET_PARAM(TypeTag, Scalar, RelaxedMaxPvFraction);
            tolerance_mb_ = EWOMS_GET_PARAM(TypeTag, Scalar, ToleranceMb);
            tolerance_cnv_ = EWOMS_GET_PARAM(TypeTag, Scalar, ToleranceCnv);
            tolerance_cnv_relaxed_ = EWOMS_GET_PARAM(TypeTag, Scalar, ToleranceCnvRelaxed);
            tolerance_wells_ = EWOMS_GET_PARAM(TypeTag, Scalar, ToleranceWells);
            tolerance_well_control_ = EWOMS_GET_PARAM(TypeTag, Scalar, ToleranceWellControl);
            max_welleq_iter_ = EWOMS_GET_PARAM(TypeTag, int, MaxWelleqIter);
            use_multisegment_well_ = EWOMS_GET_PARAM(TypeTag, bool, UseMultisegmentWell);
            tolerance_pressure_ms_wells_ = EWOMS_GET_PARAM(TypeTag, Scalar, TolerancePressureMsWells);
            relaxed_inner_tolerance_flow_ms_well_ = EWOMS_GET_PARAM(TypeTag, Scalar, RelaxedEFlowTolInnerIterMsw);
            relaxed_inner_tolerance_pressure_ms_well_ = EWOMS_GET_PARAM(TypeTag, Scalar, RelaxedPressureTolInnerIterMsw);
            max_pressure_change_ms_wells_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxPressureChangeMsWells);
            use_inner_iterations_ms_wells_ = EWOMS_GET_PARAM(TypeTag, bool, UseInnerIterationsMsWells);
            max_inner_iter_ms_wells_ = EWOMS_GET_PARAM(TypeTag, int, MaxInnerIterMsWells);
            strict_inner_iter_ms_wells_ = EWOMS_GET_PARAM(TypeTag, int, StrictInnerIterMsWells);
            regularization_factor_ms_wells_ = EWOMS_GET_PARAM(TypeTag, Scalar, RegularizationFactorMsw);
            use_inner_iterations_wells_ = EWOMS_GET_PARAM(TypeTag, bool, UseInnerIterationsWells);
            max_inner_iter_wells_ = EWOMS_GET_PARAM(TypeTag, int, MaxInnerIterWells);
            maxSinglePrecisionTimeStep_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxSinglePrecisionDays) *24*60*60;
            max_strict_iter_ = EWOMS_GET_PARAM(TypeTag, int, MaxStrictIter);
            solve_welleq_initially_ = EWOMS_GET_PARAM(TypeTag, bool, SolveWelleqInitially);
            update_equations_scaling_ = EWOMS_GET_PARAM(TypeTag, bool, UpdateEquationsScaling);
            use_update_stabilization_ = EWOMS_GET_PARAM(TypeTag, bool, UseUpdateStabilization);
            matrix_add_well_contributions_ = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);

            deck_file_name_ = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
        }

        static void registerParameters()
        {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DbhpMaxRel, "Maximum relative change of the bottom-hole pressure in a single iteration");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DwellFractionMax, "Maximum absolute change of a well's volume fraction in a single iteration");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxResidualAllowed, "Absolute maximum tolerated for residuals without cutting the time step size");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, RelaxedMaxPvFraction, "The fraction of the pore volume of the reservoir "
                                 "where the volumetric error (CNV) may be voilated during strict Newton iterations.");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, ToleranceMb, "Tolerated mass balance error relative to total mass present");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, ToleranceCnv, "Local convergence tolerance (Maximum of local saturation errors)");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, ToleranceCnvRelaxed, "Relaxed local convergence tolerance that applies for iterations after the iterations with the strict tolerance");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, ToleranceWells, "Well convergence tolerance");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, ToleranceWellControl, "Tolerance for the well control equations");
            EWOMS_REGISTER_PARAM(TypeTag, int, MaxWelleqIter, "Maximum number of iterations to determine solution the  well equations");
            EWOMS_REGISTER_PARAM(TypeTag, bool, UseMultisegmentWell, "Use the well model for multi-segment wells instead of the one for single-segment wells");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, TolerancePressureMsWells, "Tolerance for the pressure equations for multi-segment wells");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, RelaxedEFlowTolInnerIterMsw, "Relaxed tolerance for the inner iteration for the MSW flow solution");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, RelaxedPressureTolInnerIterMsw, "Relaxed tolerance for the inner iteration for the MSW pressure solution");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxPressureChangeMsWells, "Maximum relative pressure change for a single iteration of the multi-segment well model");
            EWOMS_REGISTER_PARAM(TypeTag, bool, UseInnerIterationsMsWells, "Use nested iterations for multi-segment wells");
            EWOMS_REGISTER_PARAM(TypeTag, int, MaxInnerIterMsWells, "Maximum number of inner iterations for multi-segment wells");
            EWOMS_REGISTER_PARAM(TypeTag, int, StrictInnerIterMsWells, "Number of inner iterations for multi-segment wells with strict tolerance");
            EWOMS_REGISTER_PARAM(TypeTag, bool, UseInnerIterationsWells, "Use nested iterations for standard wells");
            EWOMS_REGISTER_PARAM(TypeTag, int, MaxInnerIterWells, "Maximum number of inner iterations for standard wells");
            EWOMS_REGISTER_PARAM(TypeTag, bool, AlternativeWellRateInit, "Use alternative well rate initialization procedure");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, RegularizationFactorMsw, "Regularization factor for ms wells");
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxSinglePrecisionDays, "Maximum time step size where single precision floating point arithmetic can be used solving for the linear systems of equations");
            EWOMS_REGISTER_PARAM(TypeTag, int, MaxStrictIter, "Maximum number of Newton iterations before relaxed tolerances are used for the CNV convergence criterion");
            EWOMS_REGISTER_PARAM(TypeTag, bool, SolveWelleqInitially, "Fully solve the well equations before each iteration of the reservoir model");
            EWOMS_REGISTER_PARAM(TypeTag, bool, UpdateEquationsScaling, "Update scaling factors for mass balance equations during the run");
            EWOMS_REGISTER_PARAM(TypeTag, bool, UseUpdateStabilization, "Try to detect and correct oscillations or stagnation during the Newton method");
            EWOMS_REGISTER_PARAM(TypeTag, bool, MatrixAddWellContributions, "Explicitly specify the influences of wells between cells in the Jacobian and preconditioner matrices");
            EWOMS_REGISTER_PARAM(TypeTag, bool, EnableWellOperabilityCheck, "Enable the well operability checking");
        }
    };
} // namespace Ewoms

#endif // EWOMS_BLACKOILMODELPARAMETERS_EEBOS_HH
