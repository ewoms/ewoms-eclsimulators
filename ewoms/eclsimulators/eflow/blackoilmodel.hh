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

#ifndef EWOMS_EFLOWBLACKOILMODEL_HH
#define EWOMS_EFLOWBLACKOILMODEL_HH

#include <eebos/eclproblem.hh>
#include <ewoms/numerics/utils/start.hh>

#include <ewoms/eclsimulators/timestepping/adaptivetimestepper.hh>

#include <ewoms/eclsimulators/eflow/nonlinearsolver.hh>
#include <ewoms/eclsimulators/eflow/blackoilmodelparameters.hh>
#include <ewoms/eclsimulators/wells/blackoilwellmodel.hh>
#include <ewoms/eclsimulators/aquifers/blackoilaquifermodel.hh>
#include <ewoms/eclsimulators/wells/wellconnectionauxiliarymodule.hh>
#include <ewoms/eclsimulators/eflow/countglobalcells.hh>

#include <ewoms/eclgrids/unstructuredgrid.h>
#include <ewoms/eclsimulators/timestepping/simulatorreport.hh>
#include <ewoms/eclsimulators/linalg/parallelistlinformation.hh>
#include <ewoms/eclsimulators/deprecated/props/phaseusagefromdeck.hh>
#include <ewoms/eclio/errormacros.hh>
#include <ewoms/eclio/exceptions.hh>
#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/parser/units/units.hh>
#include <ewoms/eclsimulators/timestepping/simulatortimer.hh>
#include <ewoms/eclio/utility/parameters/parametergroup.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/tablemanager.hh>

#include <ewoms/eclsimulators/linalg/istlsolver.hh>
#include <ewoms/eclio/data/simulationdatacontainer.hh>

#include <dune/istl/owneroverlapcopy.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
#include <dune/common/parallel/communication.hh>
#else
#include <dune/common/parallel/collectivecommunication.hh>
#endif
#include <dune/common/timer.hh>
#include <dune/common/unused.hh>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <algorithm>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclEFlowProblem, INHERITS_FROM(BlackOilModel, EclBaseProblem, EFlowNonLinearSolver, EFlowModelParameters, EFlowTimeSteppingParameters));
SET_STRING_PROP(EclEFlowProblem, OutputDir, "");
SET_BOOL_PROP(EclEFlowProblem, EnableDebuggingChecks, false);
// default in eflow is to formulate the equations in surface volumes
SET_BOOL_PROP(EclEFlowProblem, BlackoilConserveSurfaceVolume, true);
SET_BOOL_PROP(EclEFlowProblem, UseVolumetricResidual, false);

SET_TYPE_PROP(EclEFlowProblem, EclAquiferModel, Ewoms::BlackoilAquiferModel<TypeTag>);

// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
SET_BOOL_PROP(EclEFlowProblem, EnablePolymer, false);
SET_BOOL_PROP(EclEFlowProblem, EnableSolvent, false);
SET_BOOL_PROP(EclEFlowProblem, EnableTemperature, true);
SET_BOOL_PROP(EclEFlowProblem, EnableEnergy, false);
SET_BOOL_PROP(EclEFlowProblem, EnableFoam, false);
SET_BOOL_PROP(EclEFlowProblem, EnableBrine, false);

SET_TYPE_PROP(EclEFlowProblem, EclWellModel, Ewoms::BlackoilWellModel<TypeTag>);
SET_TAG_PROP(EclEFlowProblem, LinearSolverSplice, EFlowIstlSolver);

END_PROPERTIES

namespace Ewoms {
    /// A model implementation for three-phase black oil.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil and vice versa. It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    template <class TypeTag>
    class BlackoilModel
    {
    public:
        // ---------  Types and enums  ---------
        typedef WellStateFullyImplicitBlackoil WellState;
        typedef BlackoilModelParameters<TypeTag> ModelParameters;

        using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
        using Grid = GET_PROP_TYPE(TypeTag, Grid);
        using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
        using SparseMatrixAdapter = GET_PROP_TYPE(TypeTag, SparseMatrixAdapter);
        using SolutionVector = GET_PROP_TYPE(TypeTag, SolutionVector);
        using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
        using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
        using Indices = GET_PROP_TYPE(TypeTag, Indices);
        using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
        using MaterialLawParams = GET_PROP_TYPE(TypeTag, MaterialLawParams);

        typedef double Scalar;
        static const int numEq = Indices::numEq;
        static const int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static const int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
        static const int contiEnergyEqIdx = Indices::contiEnergyEqIdx;
        static const int contiPolymerMWEqIdx = Indices::contiPolymerMWEqIdx;
        static const int contiFoamEqIdx = Indices::contiFoamEqIdx;
        static const int contiBrineEqIdx = Indices::contiBrineEqIdx;
        static const int solventSaturationIdx = Indices::solventSaturationIdx;
        static const int polymerConcentrationIdx = Indices::polymerConcentrationIdx;
        static const int polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;
        static const int temperatureIdx = Indices::temperatureIdx;
        static const int foamConcentrationIdx = Indices::foamConcentrationIdx;
        static const int saltConcentrationIdx = Indices::saltConcentrationIdx;

        typedef Dune::FieldVector<Scalar, numEq >        VectorBlockType;
        typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
        typedef typename SparseMatrixAdapter::IstlMatrix Mat;
        typedef Dune::BlockVector<VectorBlockType>      BVector;

        typedef ISTLSolver<TypeTag> ISTLSolverType;
        //typedef typename SolutionVector :: value_type            PrimaryVariables ;

        // ---------  Public methods  ---------

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] wells            well structure
        /// \param[in] vfp_properties   Vertical flow performance tables
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilModel(Simulator& eebosSimulator,
                          const ModelParameters& param,
                          BlackoilWellModel<TypeTag>& well_model,
                          const bool terminal_output)
        : eebosSimulator_(eebosSimulator)
        , grid_(eebosSimulator_.vanguard().grid())
        , phaseUsage_(phaseUsageFromDeck(eclState()))
        , has_disgas_(FluidSystem::enableDissolvedGas())
        , has_vapoil_(FluidSystem::enableVaporizedOil())
        , has_solvent_(GET_PROP_VALUE(TypeTag, EnableSolvent))
        , has_polymer_(GET_PROP_VALUE(TypeTag, EnablePolymer))
        , has_polymermw_(GET_PROP_VALUE(TypeTag, EnablePolymerMW))
        , has_energy_(GET_PROP_VALUE(TypeTag, EnableEnergy))
        , has_foam_(GET_PROP_VALUE(TypeTag, EnableFoam))
        , has_brine_(GET_PROP_VALUE(TypeTag, EnableBrine))
        , param_( param )
        , well_model_ (well_model)
        , terminal_output_ (terminal_output)
        , current_relaxation_(1.0)
        , dx_old_(UgGridHelpers::numCells(grid_))
        {
            // compute global sum of number of cells
            global_nc_ = detail::countGlobalCells(grid_);
            convergence_reports_.reserve(300); // Often insufficient, but avoids frequent moves.
        }

        bool isParallel() const
        { return  grid_.comm().size() > 1; }

        const EclipseState& eclState() const
        { return eebosSimulator_.vanguard().eclState(); }

        /// Called once before each time step.
        /// \param[in] timer                  simulation timer
        void prepareStep(const SimulatorTimerInterface& timer)
        {
            // update the solution variables in eebos
            if ( timer.lastStepFailed() ) {
                eebosSimulator_.model().updateFailed();
            } else {
                eebosSimulator_.model().advanceTimeLevel();
            }

            // set the timestep size and episode index for eebos explicitly. eebos needs to
            // know the report step/episode index because of timing dependend data
            // despite the fact that eflow uses its own time stepper. (The length of the
            // episode does not matter, though.)
            eebosSimulator_.setTime(timer.simulationTimeElapsed());
            eebosSimulator_.setTimeStepSize(timer.currentStepLength());
            eebosSimulator_.problem().beginTimeStep();

            unsigned numDof = eebosSimulator_.model().numGridDof();
            wasSwitched_.resize(numDof);
            std::fill(wasSwitched_.begin(), wasSwitched_.end(), false);

            if (param_.update_equations_scaling_) {
                std::cout << "equation scaling not suported yet" << std::endl;
                //updateEquationsScaling();
            }
        }

        /// Called once per nonlinear iteration.
        /// This model will perform a Newton-Raphson update, changing reservoir_state
        /// and well_state. It will also use the nonlinear_solver to do relaxation of
        /// updates if necessary.
        /// \param[in] iteration              should be 0 for the first call of a new timestep
        /// \param[in] timer                  simulation timer
        /// \param[in] nonlinear_solver       nonlinear solver used (for oscillation/relaxation control)
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIteration(const int iteration,
                                                 const SimulatorTimerInterface& timer,
                                                 NonlinearSolverType& nonlinear_solver)
        {
            SimulatorReportSingle report;
            failureReport_ = SimulatorReportSingle();
            Dune::Timer perfTimer;

            perfTimer.start();
            if (iteration == 0) {
                // For each iteration we store in a vector the norms of the residual of
                // the mass balance for each active phase, the well flux and the well equations.
                residual_norms_history_.clear();
                current_relaxation_ = 1.0;
                dx_old_ = 0.0;
                convergence_reports_.push_back({timer.reportStepNum(), timer.currentStepNum(), {}});
                convergence_reports_.back().report.reserve(11);
            }

            report.total_linearizations = 1;

            try {
                report += assembleReservoir(timer, iteration);
                report.assemble_time += perfTimer.stop();
            }
            catch (...) {
                report.assemble_time += perfTimer.stop();
                failureReport_ += report;
                // todo (?): make the report an attribute of the class
                throw; // continue throwing the stick
            }

            std::vector<double> residual_norms;
            perfTimer.reset();
            perfTimer.start();
            // the step is not considered converged until at least minIter iterations is done
            {
                auto convrep = getConvergence(timer, iteration,residual_norms);
                report.converged = convrep.converged()  && iteration > nonlinear_solver.minIter();;
                ConvergenceReport::Severity severity = convrep.severityOfWorstFailure();
                convergence_reports_.back().report.push_back(std::move(convrep));

                // Throw if any NaN or too large residual found.
                if (severity == ConvergenceReport::Severity::NotANumber) {
                    EWOMS_THROW(Ewoms::NumericalIssue, "NaN residual found!");
                } else if (severity == ConvergenceReport::Severity::TooLarge) {
                    EWOMS_THROW(Ewoms::NumericalIssue, "Too large residual found!");
                }
            }
            report.update_time += perfTimer.stop();
            residual_norms_history_.push_back(residual_norms);
            if (!report.converged) {
                perfTimer.reset();
                perfTimer.start();
                report.total_newton_iterations = 1;

                // enable single precision for solvers when dt is smaller then 20 days
                //residual_.singlePrecision = (unit::convert::to(dt, unit::day) < 20.) ;

                // Compute the nonlinear update.
                const int nc = UgGridHelpers::numCells(grid_);
                BVector x(nc);

                // apply the Schur compliment of the well model to the reservoir linearized
                // equations
                wellModel().linearize(eebosSimulator().model().linearizer().jacobian(),
                                      eebosSimulator().model().linearizer().residual());

                // Solve the linear system.
                linear_solve_setup_time_ = 0.0;
                try {
                    solveJacobianSystem(x);
                    report.linear_solve_setup_time += linear_solve_setup_time_;
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();
                }
                catch (...) {
                    report.linear_solve_setup_time += linear_solve_setup_time_;
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();

                    failureReport_ += report;
                    throw; // re-throw up
                }

                perfTimer.reset();
                perfTimer.start();

                // handling well state update before oscillation treatment is a decision based
                // on observation to avoid some big performance degeneration under some circumstances.
                // there is no theorectical explanation which way is better for sure.
                wellModel().postSolve(x);

                if (param_.use_update_stabilization_) {
                    // Stabilize the nonlinear update.
                    bool isOscillate = false;
                    bool isStagnate = false;
                    nonlinear_solver.detectOscillations(residual_norms_history_, iteration, isOscillate, isStagnate);
                    if (isOscillate) {
                        current_relaxation_ -= nonlinear_solver.relaxIncrement();
                        current_relaxation_ = std::max(current_relaxation_, nonlinear_solver.relaxMax());
                        if (terminalOutputEnabled()) {
                            std::string msg = "    Oscillating behavior detected: Relaxation set to "
                                    + std::to_string(current_relaxation_);
                            OpmLog::info(msg);
                        }
                    }
                    nonlinear_solver.stabilizeNonlinearUpdate(x, dx_old_, current_relaxation_);
                }

                // Apply the update, with considering model-dependent limitations and
                // chopping of the update.
                updateSolution(x);

                report.update_time += perfTimer.stop();
            }

            return report;
        }

        void printIf(int c, double x, double y, double eps, std::string type) {
            if (std::abs(x-y) > eps) {
                std::cout << type << " " <<c << ": "<<x << " " << y << std::endl;
            }
        }

        /// Called once after each time step.
        /// In this class, this function does nothing.
        /// \param[in] timer                  simulation timer
        void afterStep(const SimulatorTimerInterface& timer EWOMS_UNUSED)
        {
            eebosSimulator_.problem().endTimeStep();
        }

        /// Assemble the residual and Jacobian of the nonlinear system.
        /// \param[in]      reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /// \param[in]      initial_assembly  pass true if this is the first call to assemble() in this timestep
        SimulatorReportSingle assembleReservoir(const SimulatorTimerInterface& /* timer */,
                                                const int iterationIdx)
        {
            // -------- Mass balance equations --------
            eebosSimulator_.model().newtonMethod().setIterationIndex(iterationIdx);
            eebosSimulator_.problem().beginIteration();
            eebosSimulator_.model().linearizer().linearizeDomain();
            eebosSimulator_.problem().endIteration();

            return wellModel().lastReport();
        }

        // compute the "relative" change of the solution between time steps
        double relativeChange() const
        {
            Scalar resultDelta = 0.0;
            Scalar resultDenom = 0.0;

            const auto& elemMapper = eebosSimulator_.model().elementMapper();
            const auto& gridView = eebosSimulator_.gridView();
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;
                if (elem.partitionType() != Dune::InteriorEntity)
                    continue;

                unsigned globalElemIdx = elemMapper.index(elem);
                const auto& priVarsNew = eebosSimulator_.model().solution(/*timeIdx=*/0)[globalElemIdx];

                Scalar pressureNew;
                pressureNew = priVarsNew[Indices::pressureSwitchIdx];

                Scalar saturationsNew[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationNew = 1.0;
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    saturationsNew[FluidSystem::waterPhaseIdx] = priVarsNew[Indices::waterSaturationIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::waterPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && priVarsNew.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
                    saturationsNew[FluidSystem::gasPhaseIdx] = priVarsNew[Indices::compositionSwitchIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::gasPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    saturationsNew[FluidSystem::oilPhaseIdx] = oilSaturationNew;
                }

                const auto& priVarsOld = eebosSimulator_.model().solution(/*timeIdx=*/1)[globalElemIdx];

                Scalar pressureOld;
                pressureOld = priVarsOld[Indices::pressureSwitchIdx];

                Scalar saturationsOld[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationOld = 1.0;

                // NB fix me! adding pressures changes to satutation changes does not make sense
                Scalar tmp = pressureNew - pressureOld;
                resultDelta += tmp*tmp;
                resultDenom += pressureNew*pressureNew;

                if (FluidSystem::numActivePhases() > 1) {
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        saturationsOld[FluidSystem::waterPhaseIdx] = priVarsOld[Indices::waterSaturationIdx];
                        oilSaturationOld -= saturationsOld[FluidSystem::waterPhaseIdx];
                    }

                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                        priVarsOld.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg)
                    {
                        saturationsOld[FluidSystem::gasPhaseIdx] = priVarsOld[Indices::compositionSwitchIdx];
                        oilSaturationOld -= saturationsOld[FluidSystem::gasPhaseIdx];
                    }

                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                        saturationsOld[FluidSystem::oilPhaseIdx] = oilSaturationOld;
                    }
                    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++ phaseIdx) {
                        Scalar tmpSat = saturationsNew[phaseIdx] - saturationsOld[phaseIdx];
                        resultDelta += tmpSat*tmpSat;
                        resultDenom += saturationsNew[phaseIdx]*saturationsNew[phaseIdx];
                        assert(std::isfinite(resultDelta));
                        assert(std::isfinite(resultDenom));
                    }
                }
            }

            resultDelta = gridView.comm().sum(resultDelta);
            resultDenom = gridView.comm().sum(resultDenom);

            if (resultDenom > 0.0)
                return resultDelta/resultDenom;
            return 0.0;
        }

        /// Number of linear iterations used in last call to solveJacobianSystem().
        int linearIterationsLastSolve() const
        {
            return eebosSimulator_.model().newtonMethod().linearSolver().iterations ();
        }

        /// Solve the Jacobian system Jx = r where J is the Jacobian and
        /// r is the residual.
        void solveJacobianSystem(BVector& x)
        {

            auto& eebosJac = eebosSimulator_.model().linearizer().jacobian();
            auto& eebosResid = eebosSimulator_.model().linearizer().residual();

            // set initial guess
            x = 0.0;

            auto& eebosSolver = eebosSimulator_.model().newtonMethod().linearSolver();
            Dune::Timer perfTimer;
            perfTimer.start();
            eebosSolver.prepare(eebosJac, eebosResid);
            linear_solve_setup_time_ = perfTimer.stop();
            eebosSolver.setResidual(eebosResid);
            // actually, the error needs to be calculated after setResidual in order to
            // account for parallelization properly. since the residual of ECFV
            // discretizations does not need to be synchronized across processes to be
            // consistent, this is not relevant for eflow...
            eebosSolver.setMatrix(eebosJac);
            eebosSolver.solve(x);
       }

        /// Apply an update to the primary variables.
        void updateSolution(const BVector& dx)
        {
            auto& eebosNewtonMethod = eebosSimulator_.model().newtonMethod();
            SolutionVector& solution = eebosSimulator_.model().solution(/*timeIdx=*/0);

            eebosNewtonMethod.update_(/*nextSolution=*/solution,
                                     /*curSolution=*/solution,
                                     /*update=*/dx,
                                     /*resid=*/dx); // the update routines of the black
                                                    // oil model do not care about the
                                                    // residual

            // if the solution is updated, the intensive quantities need to be recalculated
            eebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
        }

        /// Return true if output to cout is wanted.
        bool terminalOutputEnabled() const
        {
            return terminal_output_;
        }

        template <class CollectiveCommunication>
        double convergenceReduction(const CollectiveCommunication& comm,
                                    const double pvSumLocal,
                                    std::vector< Scalar >& R_sum,
                                    std::vector< Scalar >& maxCoeff,
                                    std::vector< Scalar >& B_avg)
        {
            // Compute total pore volume (use only owned entries)
            double pvSum = pvSumLocal;

            if( comm.size() > 1 )
            {
                // global reduction
                std::vector< Scalar > sumBuffer;
                std::vector< Scalar > maxBuffer;
                const int numComp = B_avg.size();
                sumBuffer.reserve( 2*numComp + 1 ); // +1 for pvSum
                maxBuffer.reserve( numComp );
                for( int compIdx = 0; compIdx < numComp; ++compIdx )
                {
                    sumBuffer.push_back( B_avg[ compIdx ] );
                    sumBuffer.push_back( R_sum[ compIdx ] );
                    maxBuffer.push_back( maxCoeff[ compIdx ] );
                }

                // Compute total pore volume
                sumBuffer.push_back( pvSum );

                // compute global sum
                comm.sum( sumBuffer.data(), sumBuffer.size() );

                // compute global max
                comm.max( maxBuffer.data(), maxBuffer.size() );

                // restore values to local variables
                for( int compIdx = 0, buffIdx = 0; compIdx < numComp; ++compIdx, ++buffIdx )
                {
                    B_avg[ compIdx ]    = sumBuffer[ buffIdx ];
                    ++buffIdx;

                    R_sum[ compIdx ]       = sumBuffer[ buffIdx ];
                }

                for( int compIdx = 0; compIdx < numComp; ++compIdx )
                {
                    maxCoeff[ compIdx ] = maxBuffer[ compIdx ];
                }

                // restore global pore volume
                pvSum = sumBuffer.back();
            }

            // return global pore volume
            return pvSum;
        }

        // Get reservoir quantities on this process needed for convergence calculations.
        double localConvergenceData(std::vector<Scalar>& R_sum,
                                    std::vector<Scalar>& maxCoeff,
                                    std::vector<Scalar>& B_avg)
        {
            double pvSumLocal = 0.0;
            const auto& eebosModel = eebosSimulator_.model();
            const auto& eebosProblem = eebosSimulator_.problem();

            const auto& eebosResid = eebosSimulator_.model().linearizer().residual();

            ElementContext elemCtx(eebosSimulator_);
            const auto& gridView = eebosSimulator().gridView();
            const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();

            for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
                 elemIt != elemEndIt;
                 ++elemIt)
            {
                const auto& elem = *elemIt;
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                const double pvValue = eebosProblem.referencePorosity(cell_idx, /*timeIdx=*/0) * eebosModel.dofTotalVolume( cell_idx );
                pvSumLocal += pvValue;

                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }

                    const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));

                    B_avg[ compIdx ] += 1.0 / fs.invB(phaseIdx).value();
                    const auto R2 = eebosResid[cell_idx][compIdx];

                    R_sum[ compIdx ] += R2;
                    maxCoeff[ compIdx ] = std::max( maxCoeff[ compIdx ], std::abs( R2 ) / pvValue );
                }

                if ( has_solvent_ ) {
                    B_avg[ contiSolventEqIdx ] += 1.0 / intQuants.solventInverseFormationVolumeFactor().value();
                    const auto R2 = eebosResid[cell_idx][contiSolventEqIdx];
                    R_sum[ contiSolventEqIdx ] += R2;
                    maxCoeff[ contiSolventEqIdx ] = std::max( maxCoeff[ contiSolventEqIdx ], std::abs( R2 ) / pvValue );
                }
                if (has_polymer_ ) {
                    B_avg[ contiPolymerEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R2 = eebosResid[cell_idx][contiPolymerEqIdx];
                    R_sum[ contiPolymerEqIdx ] += R2;
                    maxCoeff[ contiPolymerEqIdx ] = std::max( maxCoeff[ contiPolymerEqIdx ], std::abs( R2 ) / pvValue );
                }
                if (has_foam_ ) {
                    B_avg[ contiFoamEqIdx ] += 1.0 / fs.invB(FluidSystem::gasPhaseIdx).value();
                    const auto R2 = eebosResid[cell_idx][contiFoamEqIdx];
                    R_sum[ contiFoamEqIdx ] += R2;
                    maxCoeff[ contiFoamEqIdx ] = std::max( maxCoeff[ contiFoamEqIdx ], std::abs( R2 ) / pvValue );
                }
                if (has_brine_ ) {
                    B_avg[ contiBrineEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R2 = eebosResid[cell_idx][contiBrineEqIdx];
                    R_sum[ contiBrineEqIdx ] += R2;
                    maxCoeff[ contiBrineEqIdx ] = std::max( maxCoeff[ contiBrineEqIdx ], std::abs( R2 ) / pvValue );
                }

                if (has_polymermw_) {
                    assert(has_polymer_);

                    B_avg[contiPolymerMWEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    // the residual of the polymer molecular equation is scaled down by a 100, since molecular weight
                    // can be much bigger than 1, and this equation shares the same tolerance with other mass balance equations
                    // TODO: there should be a more general way to determine the scaling-down coefficient
                    const auto R2 = eebosResid[cell_idx][contiPolymerMWEqIdx] / 100.;
                    R_sum[contiPolymerMWEqIdx] += R2;
                    maxCoeff[contiPolymerMWEqIdx] = std::max( maxCoeff[contiPolymerMWEqIdx], std::abs( R2 ) / pvValue );
                }

                if (has_energy_ ) {
                    B_avg[ contiEnergyEqIdx ] += 1.0;
                    const auto R2 = eebosResid[cell_idx][contiEnergyEqIdx];
                    R_sum[ contiEnergyEqIdx ] += R2;
                    maxCoeff[ contiEnergyEqIdx ] = std::max( maxCoeff[ contiEnergyEqIdx ], std::abs( R2 ) / pvValue );
                }

            }

            // compute local average in terms of global number of elements
            const int bSize = B_avg.size();
            for ( int i = 0; i<bSize; ++i )
            {
                B_avg[ i ] /= Scalar( global_nc_ );
            }

            return pvSumLocal;
        }

        double computeCnvErrorPv(const std::vector<Scalar>& B_avg, double dt)
        {
            double errorPV{};
            const auto& eebosModel = eebosSimulator_.model();
            const auto& eebosProblem = eebosSimulator_.problem();
            const auto& eebosResid = eebosSimulator_.model().linearizer().residual();
            const auto& gridView = eebosSimulator().gridView();
            ElementContext elemCtx(eebosSimulator_);

            for (const auto& elem: elements(gridView, Dune::Partitions::interiorBorder))
            {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const double pvValue = eebosProblem.referencePorosity(cell_idx, /*timeIdx=*/0) * eebosModel.dofTotalVolume( cell_idx );
                const auto& cellResidual = eebosResid[cell_idx];
                bool cnvViolated = false;

                for (unsigned eqIdx = 0; eqIdx < cellResidual.size(); ++eqIdx)
                {
                    using std::abs;
                    Scalar CNV = cellResidual[eqIdx] * dt * B_avg[eqIdx] / pvValue;
                    cnvViolated = cnvViolated || (abs(CNV) > param_.tolerance_cnv_);
                }

                if (cnvViolated)
                {
                    errorPV += pvValue;
                }
            }

            return grid_.comm().sum(errorPV);
        }

        ConvergenceReport getReservoirConvergence(const double dt,
                                                  const int iteration,
                                                  std::vector<Scalar>& B_avg,
                                                  std::vector<Scalar>& residual_norms)
        {
            typedef std::vector< Scalar > Vector;

            const int numComp = numEq;
            Vector R_sum(numComp, 0.0 );
            Vector maxCoeff(numComp, std::numeric_limits< Scalar >::lowest() );
            const double pvSumLocal = localConvergenceData(R_sum, maxCoeff, B_avg);

            // compute global sum and max of quantities
            const double pvSum = convergenceReduction(grid_.comm(), pvSumLocal,
                                                      R_sum, maxCoeff, B_avg);

            auto cnvErrorPvFraction = computeCnvErrorPv(B_avg, dt);
            cnvErrorPvFraction /= pvSum;

            const double tol_mb  = param_.tolerance_mb_;
            // Default value of relaxed_max_pv_fraction_ is 1 and
            // max_strict_iter_ is 8. Hence only iteration chooses
            // whether to use relaxed or not.
            // To activate only fraction use fraction below 1 and iter 0.
            const bool use_relaxed = cnvErrorPvFraction < param_.relaxed_max_pv_fraction_ && iteration >= param_.max_strict_iter_;
            const double tol_cnv = use_relaxed ? param_.tolerance_cnv_relaxed_ :  param_.tolerance_cnv_;

            // Finish computation
            std::vector<Scalar> CNV(numComp);
            std::vector<Scalar> mass_balance_residual(numComp);
            for ( int compIdx = 0; compIdx < numComp; ++compIdx )
            {
                CNV[compIdx]                    = B_avg[compIdx] * dt * maxCoeff[compIdx];
                mass_balance_residual[compIdx]  = std::abs(B_avg[compIdx]*R_sum[compIdx]) * dt / pvSum;
                residual_norms.push_back(CNV[compIdx]);
            }

            // Setup component names, only the first time the function is run.
            static std::vector<std::string> compNames;
            if (compNames.empty()) {
                compNames.resize(numComp);
                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }
                    const unsigned canonicalCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
                    const unsigned compIdx = Indices::canonicalToActiveComponentIndex(canonicalCompIdx);
                    compNames[compIdx] = FluidSystem::componentName(canonicalCompIdx);
                }
                if (has_solvent_) {
                    compNames[solventSaturationIdx] = "Solvent";
                }
                if (has_polymer_) {
                    compNames[polymerConcentrationIdx] = "Polymer";
                }
                if (has_polymermw_) {
                    assert(has_polymer_);
                    compNames[polymerMoleWeightIdx] = "MolecularWeightP";
                }
                if (has_energy_) {
                    compNames[temperatureIdx] = "Energy";
                }
                if (has_foam_) {
                    compNames[foamConcentrationIdx] = "Foam";
                }
                if (has_brine_) {
                    compNames[saltConcentrationIdx] = "Brine";
                }
            }

            // Create convergence report.
            ConvergenceReport report;
            using CR = ConvergenceReport;
            for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                double res[2] = { mass_balance_residual[compIdx], CNV[compIdx] };
                CR::ReservoirFailure::Type types[2] = { CR::ReservoirFailure::Type::MassBalance,
                                                        CR::ReservoirFailure::Type::Cnv };
                double tol[2] = { tol_mb, tol_cnv };
                for (int ii : {0, 1}) {
                    if (std::isnan(res[ii])) {
                        report.setReservoirFailed({types[ii], CR::Severity::NotANumber, compIdx});
                        if ( terminal_output_ ) {
                            OpmLog::debug("NaN residual for " + compNames[compIdx] + " equation.");
                        }
                    } else if (res[ii] > maxResidualAllowed()) {
                        report.setReservoirFailed({types[ii], CR::Severity::TooLarge, compIdx});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Too large residual for " + compNames[compIdx] + " equation.");
                        }
                    } else if (res[ii] < 0.0) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Negative residual for " + compNames[compIdx] + " equation.");
                        }
                    } else if (res[ii] > tol[ii]) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                    }
                }
            }

            // Output of residuals.
            if ( terminal_output_ )
            {
                // Only rank 0 does print to std::cout
                if (iteration == 0) {
                    std::string msg = "Iter";
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    MB(";
                        msg += compNames[compIdx][0];
                        msg += ")  ";
                    }
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    CNV(";
                        msg += compNames[compIdx][0];
                        msg += ") ";
                    }
                    OpmLog::debug(msg);
                }
                std::ostringstream ss;
                const std::streamsize oprec = ss.precision(3);
                const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
                ss << std::setw(4) << iteration;
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << mass_balance_residual[compIdx];
                }
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << CNV[compIdx];
                }
                ss.precision(oprec);
                ss.flags(oflags);
                OpmLog::debug(ss.str());
            }

            return report;
        }

        /// Compute convergence based on total mass balance (tol_mb) and maximum
        /// residual mass balance (tol_cnv).
        /// \param[in]   timer       simulation timer
        /// \param[in]   iteration   current iteration number
        /// \param[out]  residual_norms   CNV residuals by phase
        ConvergenceReport getConvergence(const SimulatorTimerInterface& timer,
                                         const int iteration,
                                         std::vector<double>& residual_norms)
        {
            // Get convergence reports for reservoir and wells.
            std::vector<Scalar> B_avg(numEq, 0.0);
            auto report = getReservoirConvergence(timer.currentStepLength(), iteration, B_avg, residual_norms);
            report += wellModel().getWellConvergence(B_avg);

            return report;
        }

        /// The number of active fluid phases in the model.
        int numPhases() const
        {
            return phaseUsage_.num_phases;
        }

        /// Wrapper required due to not following generic API
        template<class T>
        std::vector<std::vector<double> >
        computeFluidInPlace(const T&, const std::vector<int>& fipnum) const
        {
            return computeFluidInPlace(fipnum);
        }

        /// Should not be called
        std::vector<std::vector<double> >
        computeFluidInPlace(const std::vector<int>& /*fipnum*/) const
        {
            //assert(true)
            //return an empty vector
            std::vector<std::vector<double> > regionValues(0, std::vector<double>(0,0.0));
            return regionValues;
        }

        const Simulator& eebosSimulator() const
        { return eebosSimulator_; }

        Simulator& eebosSimulator()
        { return eebosSimulator_; }

        /// return the statistics if the nonlinearIteration() method failed
        const SimulatorReportSingle& failureReport() const
        { return failureReport_; }

        struct StepReport
        {
            int report_step;
            int current_step;
            std::vector<ConvergenceReport> report;
        };

        const std::vector<StepReport>& stepReports() const
        {
            return convergence_reports_;
        }

    protected:
        // ---------  Data members  ---------

        Simulator& eebosSimulator_;
        const Grid&            grid_;
        const PhaseUsage phaseUsage_;
        const bool has_disgas_;
        const bool has_vapoil_;
        const bool has_solvent_;
        const bool has_polymer_;
        const bool has_polymermw_;
        const bool has_energy_;
        const bool has_foam_;
        const bool has_brine_;

        ModelParameters                 param_;
        SimulatorReportSingle failureReport_;

        // Well Model
        BlackoilWellModel<TypeTag>& well_model_;

        /// \brief Whether we print something to std::cout
        bool terminal_output_;
        /// \brief The number of cells of the global grid.
        long int global_nc_;

        std::vector<std::vector<double>> residual_norms_history_;
        double current_relaxation_;
        BVector dx_old_;

        std::vector<StepReport> convergence_reports_;
    public:
        /// return the StandardWells object
        BlackoilWellModel<TypeTag>&
        wellModel() { return well_model_; }

        const BlackoilWellModel<TypeTag>&
        wellModel() const { return well_model_; }

        void beginReportStep()
        {
            eebosSimulator_.problem().beginEpisode();
        }

        void endReportStep()
        {
            eebosSimulator_.problem().endEpisode();
        }

    private:

        double dpMaxRel() const { return param_.dp_max_rel_; }
        double dsMax() const { return param_.ds_max_; }
        double drMaxRel() const { return param_.dr_max_rel_; }
        double maxResidualAllowed() const { return param_.max_residual_allowed_; }
        double linear_solve_setup_time_;
    public:
        std::vector<bool> wasSwitched_;
    };
} // namespace Ewoms

#endif // EWOMS_BLACKOILMODELBASE_IMPL_HH
