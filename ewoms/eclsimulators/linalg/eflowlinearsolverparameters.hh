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

#ifndef EWOMS_EFLOWLINEARSOLVERPARAMETERS_HH
#define EWOMS_EFLOWLINEARSOLVERPARAMETERS_HH

#include <ewoms/eclio/utility/parameters/parametergroup.hh>
#include <ewoms/eclsimulators/linalg/paralleloverlappingilu0.hh>

#include <ewoms/common/parametersystem.hh>

#include <array>
#include <memory>

namespace Ewoms {
template <class TypeTag>
class ISTLSolver;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EFlowIstlSolverParams);

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(LinearSolverReduction);
NEW_PROP_TAG(IluRelaxation);
NEW_PROP_TAG(LinearSolverMaxIter);
NEW_PROP_TAG(LinearSolverRestart);
NEW_PROP_TAG(EflowLinearSolverVerbosity);
NEW_PROP_TAG(IluFillinLevel);
NEW_PROP_TAG(MiluVariant);
NEW_PROP_TAG(IluRedblack);
NEW_PROP_TAG(IluReorderSpheres);
NEW_PROP_TAG(UseGmres);
NEW_PROP_TAG(LinearSolverRequireFullSparsityPattern);
NEW_PROP_TAG(LinearSolverIgnoreConvergenceFailure);
NEW_PROP_TAG(PreconditionerAddWellContributions);
NEW_PROP_TAG(ScaleLinearSystem);
NEW_PROP_TAG(CprMaxEllIter);
NEW_PROP_TAG(CprEllSolvetype);
NEW_PROP_TAG(CprReuseSetup);
NEW_PROP_TAG(Linsolver);
NEW_PROP_TAG(GpuMode);
NEW_PROP_TAG(BdaDeviceId);
NEW_PROP_TAG(OpenclPlatformId);
NEW_PROP_TAG(LinearSolverBackend);

SET_SCALAR_PROP(EFlowIstlSolverParams, LinearSolverReduction, 1e-2);
SET_SCALAR_PROP(EFlowIstlSolverParams, IluRelaxation, 0.9);
SET_INT_PROP(EFlowIstlSolverParams, LinearSolverMaxIter, 200);
SET_INT_PROP(EFlowIstlSolverParams, LinearSolverRestart, 40);
SET_INT_PROP(EFlowIstlSolverParams, EflowLinearSolverVerbosity, 0);
SET_INT_PROP(EFlowIstlSolverParams, IluFillinLevel, 0);
SET_STRING_PROP(EFlowIstlSolverParams, MiluVariant, "ILU");
SET_BOOL_PROP(EFlowIstlSolverParams, IluRedblack, false);
SET_BOOL_PROP(EFlowIstlSolverParams, IluReorderSpheres, false);
SET_BOOL_PROP(EFlowIstlSolverParams, UseGmres, false);
SET_BOOL_PROP(EFlowIstlSolverParams, LinearSolverRequireFullSparsityPattern, false);
SET_BOOL_PROP(EFlowIstlSolverParams, LinearSolverIgnoreConvergenceFailure, false);
SET_TYPE_PROP(EFlowIstlSolverParams, LinearSolverBackend, Ewoms::ISTLSolver<TypeTag>);
SET_BOOL_PROP(EFlowIstlSolverParams, PreconditionerAddWellContributions, false);
SET_BOOL_PROP(EFlowIstlSolverParams, ScaleLinearSystem, false);
SET_INT_PROP(EFlowIstlSolverParams, CprMaxEllIter, 20);
SET_INT_PROP(EFlowIstlSolverParams, CprEllSolvetype, 0);
SET_INT_PROP(EFlowIstlSolverParams, CprReuseSetup, 3);
SET_STRING_PROP(EFlowIstlSolverParams, Linsolver, "ilu0");
SET_STRING_PROP(EFlowIstlSolverParams, GpuMode, "none");
SET_INT_PROP(EFlowIstlSolverParams, BdaDeviceId, 0);
SET_INT_PROP(EFlowIstlSolverParams, OpenclPlatformId, 0);

END_PROPERTIES

namespace Ewoms
{

    /// This class carries all parameters for the NewtonIterationBlackoilInterleaved class.
    struct EFlowLinearSolverParameters
    {
        double linear_solver_reduction_;
        double ilu_relaxation_;
        int    linear_solver_maxiter_;
        int    linear_solver_restart_;
        int    linear_solver_verbosity_;
        int    ilu_fillin_level_;
        Ewoms::MILU_VARIANT   ilu_milu_;
        bool   ilu_redblack_;
        bool   ilu_reorder_sphere_;
        bool   newton_use_gmres_;
        bool   require_full_sparsity_pattern_;
        bool   ignoreConvergenceFailure_;
        bool scale_linear_system_;
        std::string linsolver_;
        std::string gpu_mode_;
        int bda_device_id_;
        int opencl_platform_id_;
        int cpr_max_ell_iter_ = 20;
        int cpr_reuse_setup_ = 0;
        bool use_gpu_;

        template <class TypeTag>
        void init()
        {
            // TODO: these parameters have undocumented non-trivial dependencies
            linear_solver_reduction_ = EWOMS_GET_PARAM(TypeTag, double, LinearSolverReduction);
            ilu_relaxation_ = EWOMS_GET_PARAM(TypeTag, double, IluRelaxation);
            linear_solver_maxiter_ = EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIter);
            linear_solver_restart_ = EWOMS_GET_PARAM(TypeTag, int, LinearSolverRestart);
            linear_solver_verbosity_ = EWOMS_GET_PARAM(TypeTag, int, EflowLinearSolverVerbosity);
            ilu_fillin_level_ = EWOMS_GET_PARAM(TypeTag, int, IluFillinLevel);
            ilu_milu_ = convertString2Milu(EWOMS_GET_PARAM(TypeTag, std::string, MiluVariant));
            ilu_redblack_ = EWOMS_GET_PARAM(TypeTag, bool, IluRedblack);
            ilu_reorder_sphere_ = EWOMS_GET_PARAM(TypeTag, bool, IluReorderSpheres);
            newton_use_gmres_ = EWOMS_GET_PARAM(TypeTag, bool, UseGmres);
            require_full_sparsity_pattern_ = EWOMS_GET_PARAM(TypeTag, bool, LinearSolverRequireFullSparsityPattern);
            ignoreConvergenceFailure_ = EWOMS_GET_PARAM(TypeTag, bool, LinearSolverIgnoreConvergenceFailure);
            scale_linear_system_ = EWOMS_GET_PARAM(TypeTag, bool, ScaleLinearSystem);
            cpr_max_ell_iter_  =  EWOMS_GET_PARAM(TypeTag, int, CprMaxEllIter);
            cpr_reuse_setup_  =  EWOMS_GET_PARAM(TypeTag, int, CprReuseSetup);
            linsolver_ = EWOMS_GET_PARAM(TypeTag, std::string, Linsolver);
            gpu_mode_ = EWOMS_GET_PARAM(TypeTag, std::string, GpuMode);
            bda_device_id_ = EWOMS_GET_PARAM(TypeTag, int, BdaDeviceId);
            opencl_platform_id_ = EWOMS_GET_PARAM(TypeTag, int, OpenclPlatformId);
        }

        template <class TypeTag>
        static void registerParameters()
        {
            EWOMS_REGISTER_PARAM(TypeTag, double, LinearSolverReduction, "The minimum reduction of the residual which the linear solver must achieve");
            EWOMS_REGISTER_PARAM(TypeTag, double, IluRelaxation, "The relaxation factor of the linear solver's ILU preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverMaxIter, "The maximum number of iterations of the linear solver");
            EWOMS_REGISTER_PARAM(TypeTag, int, LinearSolverRestart, "The number of iterations after which GMRES is restarted");
            EWOMS_REGISTER_PARAM(TypeTag, int, EflowLinearSolverVerbosity, "The verbosity level of the linear solver (0: off, 2: all)");
            EWOMS_REGISTER_PARAM(TypeTag, int, IluFillinLevel, "The fill-in level of the linear solver's ILU preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, MiluVariant, "Specify which variant of the modified-ILU preconditioner ought to be used. Possible variants are: ILU (default, plain ILU), MILU_1 (lump diagonal with dropped row entries), MILU_2 (lump diagonal with the sum of the absolute values of the dropped row  entries), MILU_3 (if diagonal is positive add sum of dropped row entrires. Otherwise substract them), MILU_4 (if diagonal is positive add sum of dropped row entrires. Otherwise do nothing");
            EWOMS_REGISTER_PARAM(TypeTag, bool, IluRedblack, "Use red-black partioning for the ILU preconditioner");
            EWOMS_REGISTER_PARAM(TypeTag, bool, IluReorderSpheres, "Whether to reorder the entries of the matrix in the red-black ILU preconditioner in spheres starting at an edge. If false the original ordering is preserved in each color. Otherwise why try to ensure D4 ordering (in a 2D structured grid, the diagonal elements are consecutive).");
            EWOMS_REGISTER_PARAM(TypeTag, bool, UseGmres, "Use GMRES as the linear solver");
            EWOMS_REGISTER_PARAM(TypeTag, bool, LinearSolverRequireFullSparsityPattern, "Produce the full sparsity pattern for the linear solver");
            EWOMS_REGISTER_PARAM(TypeTag, bool, LinearSolverIgnoreConvergenceFailure, "Continue with the simulation like nothing happened after the linear solver did not converge");
            EWOMS_REGISTER_PARAM(TypeTag, bool, ScaleLinearSystem, "Scale linear system according to equation scale and primary variable types");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprMaxEllIter, "MaxIterations of the elliptic pressure part of the cpr solver");
            EWOMS_REGISTER_PARAM(TypeTag, int, CprReuseSetup, "Reuse preconditioner setup. Valid options are 0: recreate the preconditioner for every linear solve, 1: recreate once every timestep, 2: recreate if last linear solve took more than 10 iterations, 3: never recreate");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, Linsolver, "Configuration of solver. Valid options are: ilu0 (default), cpr (an alias for cpr_trueimpes), cpr_quasiimpes, cpr_trueimpes or amg. Alternatively, you can request a configuration to be read from a JSON file by giving the filename here, ending with '.json.'");
            EWOMS_REGISTER_PARAM(TypeTag, std::string, GpuMode, "Use GPU cusparseSolver or openclSolver as the linear solver, usage: '--gpu-mode=[none|cusparse|opencl]'");
            EWOMS_REGISTER_PARAM(TypeTag, int, BdaDeviceId, "Choose device ID for cusparseSolver or openclSolver, use 'nvidia-smi' or 'clinfo' to determine valid IDs");
            EWOMS_REGISTER_PARAM(TypeTag, int, OpenclPlatformId, "Choose platform ID for openclSolver, use 'clinfo' to determine valid platform IDs");
        }

        EFlowLinearSolverParameters() { reset(); }

        // set default values
        void reset()
        {
            newton_use_gmres_        = false;
            linear_solver_reduction_ = 1e-2;
            linear_solver_maxiter_   = 150;
            linear_solver_restart_   = 40;
            linear_solver_verbosity_ = 0;
            require_full_sparsity_pattern_ = false;
            ignoreConvergenceFailure_ = false;
            ilu_fillin_level_         = 0;
            ilu_relaxation_           = 0.9;
            ilu_milu_                 = MILU_VARIANT::ILU;
            ilu_redblack_             = false;
            ilu_reorder_sphere_       = true;
            gpu_mode_                 = "none";
            bda_device_id_            = 0;
            opencl_platform_id_       = 0;
        }
    };

} // namespace Ewoms

#endif // EWOMS_EFLOWLINEARSOLVERPARAMETERS_HH
