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

#ifndef EWOMS_SIMULATORFULLYIMPLICITBLACKOILEEBOS_HH
#define EWOMS_SIMULATORFULLYIMPLICITBLACKOILEEBOS_HH

#include <ewoms/eclsimulators/eflow/nonlinearsolver.hh>
#include <ewoms/eclsimulators/eflow/blackoilmodel.hh>
#include <ewoms/eclsimulators/eflow/blackoilmodelparameters.hh>
#include <ewoms/eclsimulators/wells/wellstatefullyimplicitblackoil.hh>
#include <ewoms/eclsimulators/aquifers/blackoilaquifermodel.hh>
#include <ewoms/eclsimulators/timestepping/adaptivetimestepper.hh>
#include <ewoms/eclgrids/utility/stopwatch.hh>

#include <ewoms/eclio/exceptions.hh>
#include <ewoms/eclio/errormacros.hh>

BEGIN_PROPERTIES

NEW_PROP_TAG(EnableTerminalOutput);
NEW_PROP_TAG(EnableAdaptiveTimeStepping);
NEW_PROP_TAG(EnableTuning);

SET_BOOL_PROP(EclEFlowProblem, EnableTerminalOutput, true);
SET_BOOL_PROP(EclEFlowProblem, EnableAdaptiveTimeStepping, true);
SET_BOOL_PROP(EclEFlowProblem, EnableTuning, false);

END_PROPERTIES

namespace Ewoms {

/// a simulator for the blackoil model
template<class TypeTag>
class SimulatorFullyImplicitBlackoil
{
public:
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using Grid = GET_PROP_TYPE(TypeTag, Grid);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using BlackoilIndices = GET_PROP_TYPE(TypeTag, Indices);
    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
    using SolutionVector = GET_PROP_TYPE(TypeTag, SolutionVector);
    using MaterialLawParams = GET_PROP_TYPE(TypeTag, MaterialLawParams);

    typedef AdaptiveTimeStepper<TypeTag> TimeStepper;
    typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

    typedef WellStateFullyImplicitBlackoil WellState;
    typedef BlackoilModel<TypeTag> Model;
    typedef NonlinearSolver<TypeTag, Model> Solver;
    typedef typename Model::ModelParameters ModelParameters;
    typedef typename Solver::SolverParameters SolverParameters;
    typedef BlackoilWellModel<TypeTag> WellModel;
    typedef BlackoilAquiferModel<TypeTag> AquiferModel;

    /// Initialise from parameters and objects to observe.
    /// \param[in] param       parameters, this class accepts the following:
    ///     parameter (default)            effect
    ///     -----------------------------------------------------------
    ///     output (true)                  write output to files?
    ///     output_dir ("output")          output directoty
    ///     output_interval (1)            output every nth step
    ///     nl_pressure_residual_tolerance (0.0) pressure solver residual tolerance (in Pascal)
    ///     nl_pressure_change_tolerance (1.0)   pressure solver change tolerance (in Pascal)
    ///     nl_pressure_maxiter (10)       max nonlinear iterations in pressure
    ///     nl_maxiter (30)                max nonlinear iterations in transport
    ///     nl_tolerance (1e-9)            transport solver absolute residual tolerance
    ///     num_transport_substeps (1)     number of transport steps per pressure step
    ///     use_segregation_split (false)  solve for gravity segregation (if false,
    ///                                    segregation is ignored).
    ///
    /// \param[in] props         fluid and rock properties
    /// \param[in] linsolver     linear solver
    /// \param[in] eclipse_state the object which represents an internalized ECL deck
    /// \param[in] output_writer
    /// \param[in] threshold_pressures_by_face   if nonempty, threshold pressures that inhibit flow
    SimulatorFullyImplicitBlackoil(Simulator& eebosSimulator)
        : eebosSimulator_(eebosSimulator)
    {
        phaseUsage_ = phaseUsageFromDeck(eclState());

        // Only rank 0 does print to std::cout
        const auto& comm = grid().comm();
        terminalOutput_ = EWOMS_GET_PARAM(TypeTag, bool, EnableTerminalOutput);
        terminalOutput_ = terminalOutput_ && (comm.rank() == 0);
    }

    static void registerParameters()
    {
        ModelParameters::registerParameters();
        SolverParameters::registerParameters();
        TimeStepper::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTerminalOutput,
                             "Print high-level information about the simulation's progress to the terminal");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableAdaptiveTimeStepping,
                             "Use adaptive time stepping between report steps");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTuning,
                             "Honor some aspects of the TUNING keyword.");
    }

    /// Run the simulation.
    /// This will run succesive timesteps until timer.done() is true. It will
    /// modify the reservoir and well states.
    /// \param[in,out] timer       governs the requested reporting timesteps
    /// \param[in,out] state       state of reservoir: pressure, fluxes
    /// \return                    simulation report, with timing data
    SimulatorReport run(SimulatorTimer& timer)
    {
        init(timer);
        // Main simulation loop.
        while (!timer.done()) {
            bool continue_looping = runStep(timer);
            if (!continue_looping) break;
        }
        return finalize();
    }

    void init(SimulatorTimer &timer)
    {
        eebosSimulator_.setEpisodeIndex(-1);

        // Create timers and file for writing timing info.
        solverTimer_ = std::make_unique<Ewoms::time::StopWatch>();
        totalTimer_ = std::make_unique<Ewoms::time::StopWatch>();
        totalTimer_->start();

        // adaptive time stepping
        bool enableAdaptive = EWOMS_GET_PARAM(TypeTag, bool, EnableAdaptiveTimeStepping);
        bool enableTUNING = EWOMS_GET_PARAM(TypeTag, bool, EnableTuning);
        if (enableAdaptive) {
            if (enableTUNING) {
                adaptiveTimeStepping_ = std::make_unique<TimeStepper>(
                    schedule().getTuning(timer.currentStepNum()), terminalOutput_);
            }
            else {
                adaptiveTimeStepping_ = std::make_unique<TimeStepper>(terminalOutput_);
            }

            if (isRestart()) {
                // For restarts the eebosSimulator may have gotten some information
                // about the next timestep size from the OPMEXTRA field
                adaptiveTimeStepping_->setSuggestedNextStep(eebosSimulator_.timeStepSize());
            }
        }
    }

    bool runStep(SimulatorTimer& timer)
    {
        if (schedule().exitStatus().has_value()) {
            if (terminalOutput_) {
                OpmLog::info("Stopping simulation since EXIT was triggered by an action keyword.");
            }
            report_.success.exit_status = schedule().exitStatus().value();
            return false;
        }

        // Report timestep.
        if (terminalOutput_) {
            std::ostringstream ss;
            timer.report(ss);
            OpmLog::debug(ss.str());
        }

        if (terminalOutput_) {
            std::ostringstream stepMsg;
            boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
            stepMsg.imbue(std::locale(std::locale::classic(), facet));
            stepMsg << "\nReport step " << std::setw(2) <<timer.currentStepNum()
                    << "/" << timer.numSteps()
                    << " at day " << (double)unit::convert::to(timer.simulationTimeElapsed(), unit::day)
                    << "/" << (double)unit::convert::to(timer.totalTime(), unit::day)
                    << ", date = " << timer.currentDateTime();
            OpmLog::info(stepMsg.str());
        }

        // write the inital state at the report stage
        if (timer.initialStep()) {
            Dune::Timer perfTimer;
            perfTimer.start();

            eebosSimulator_.setEpisodeIndex(-1);
            eebosSimulator_.setEpisodeLength(0.0);
            eebosSimulator_.setTimeStepSize(0.0);

            wellModel_().beginReportStep(timer.currentStepNum());
            eebosSimulator_.problem().writeOutput();

            report_.success.output_write_time += perfTimer.stop();
        }

        // Run a multiple steps of the solver depending on the time step control.
        solverTimer_->start();

        auto solver = createSolver(wellModel_());

        eebosSimulator_.startNextEpisode(
            eebosSimulator_.startTime()
               + schedule().getTimeMap().getTimePassedUntil(timer.currentStepNum()),
            timer.currentStepLength());
        eebosSimulator_.setEpisodeIndex(timer.currentStepNum());
        solver->model().beginReportStep();
        const auto& events = schedule().getEvents();
        bool enableTUNING = EWOMS_GET_PARAM(TypeTag, bool, EnableTuning);

        // If sub stepping is enabled allow the solver to sub cycle
        // in case the report steps are too large for the solver to converge
        //
        // \Note: The report steps are met in any case
        // \Note: The sub stepping will require a copy of the state variables
        if (adaptiveTimeStepping_) {
            if (enableTUNING) {
                if (events.hasEvent(ScheduleEvents::TUNING_CHANGE,timer.currentStepNum())) {
                    adaptiveTimeStepping_->updateTUNING(schedule().getTuning(timer.currentStepNum()));
                }
            }

            bool event = events.hasEvent(ScheduleEvents::NEW_WELL, timer.currentStepNum()) ||
                events.hasEvent(ScheduleEvents::PRODUCTION_UPDATE, timer.currentStepNum()) ||
                events.hasEvent(ScheduleEvents::INJECTION_UPDATE, timer.currentStepNum()) ||
                events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE, timer.currentStepNum());
            auto stepReport = adaptiveTimeStepping_->step(timer, *solver, event, nullptr);
            report_ += stepReport;
        } else {
            // solve for complete report step
            auto stepReport = solver->step(timer);
            report_ += stepReport;
            if (terminalOutput_) {
                std::ostringstream ss;
                stepReport.reportStep(ss);
                OpmLog::info(ss.str());
            }
        }

        // write simulation state at the report stage
        Dune::Timer perfTimer;
        perfTimer.start();
        const double nextstep = adaptiveTimeStepping_ ? adaptiveTimeStepping_->suggestedNextStep() : -1.0;
        eebosSimulator_.problem().setNextTimeStepSize(nextstep);
        eebosSimulator_.problem().writeOutput();
        report_.success.output_write_time += perfTimer.stop();

        solver->model().endReportStep();

        // take time that was used to solve system for this reportStep
        solverTimer_->stop();

        // update timing.
        report_.success.solver_time += solverTimer_->secsSinceStart();

        // Increment timer, remember well state.
        ++timer;

        if (terminalOutput_) {
            if (!timer.initialStep()) {
                const std::string version = EWOMS_ECLSIMULATORS_VERSION;
                outputTimestampFIP(timer, version);
            }
        }

        if (terminalOutput_) {
            std::string msg =
                "Time step took " + std::to_string(solverTimer_->secsSinceStart()) + " seconds; "
                "total solver time " + std::to_string(report_.success.solver_time) + " seconds.";
            OpmLog::debug(msg);
        }
        return true;
    }

    SimulatorReport finalize()
    {
        // make sure all output is written to disk before run is finished
        {
            Dune::Timer finalOutputTimer;
            finalOutputTimer.start();

            eebosSimulator_.problem().finalizeOutput();
            report_.success.output_write_time += finalOutputTimer.stop();
        }

        // Stop timer and create timing report
        totalTimer_->stop();
        report_.success.total_time = totalTimer_->secsSinceStart();
        report_.success.converged = true;

        return report_;
    }

    const Grid& grid() const
    { return eebosSimulator_.vanguard().grid(); }

protected:

    std::unique_ptr<Solver> createSolver(WellModel& wellModel)
    {
        auto model = std::unique_ptr<Model>(new Model(eebosSimulator_,
                                                      modelParam_,
                                                      wellModel,
                                                      terminalOutput_));

        return std::unique_ptr<Solver>(new Solver(solverParam_, std::move(model)));
    }

    void outputTimestampFIP(const SimulatorTimer& timer, const std::string version)
    {
        std::ostringstream ss;
        boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d %b %Y");
        ss.imbue(std::locale(std::locale::classic(), facet));
        ss << "\n                              **************************************************************************\n"
        << "  Balance  at" << std::setw(10) << (double)unit::convert::to(timer.simulationTimeElapsed(), unit::day) << "  Days"
        << " *" << std::setw(30) << eclState().getTitle() << "                                          *\n"
        << "  Report " << std::setw(4) << timer.reportStepNum() << "    " << timer.currentDateTime()
        << "  *                                             EFlow  version " << std::setw(11) << version << "  *\n"
        << "                              **************************************************************************\n";
        OpmLog::note(ss.str());
    }

    const EclipseState& eclState() const
    { return eebosSimulator_.vanguard().eclState(); }

    const Schedule& schedule() const
    { return eebosSimulator_.vanguard().schedule(); }

    bool isRestart() const
    {
        const auto& initconfig = eclState().getInitConfig();
        return initconfig.restartRequested();
    }

    WellModel& wellModel_()
    { return eebosSimulator_.problem().wellModel(); }

    const WellModel& wellModel_() const
    { return eebosSimulator_.problem().wellModel(); }

    // Data.
    Simulator& eebosSimulator_;
    std::unique_ptr<WellConnectionAuxiliaryModule<TypeTag>> wellAuxMod_;

    ModelParameters modelParam_;
    SolverParameters solverParam_;

    // Observed objects.
    PhaseUsage phaseUsage_;
    // Misc. data
    bool terminalOutput_;

    SimulatorReport report_;
    std::unique_ptr<Ewoms::time::StopWatch> solverTimer_;
    std::unique_ptr<Ewoms::time::StopWatch> totalTimer_;
    std::unique_ptr<TimeStepper> adaptiveTimeStepping_;
};

} // namespace Ewoms

#endif // EWOMS_SIMULATOR_FULLY_IMPLICIT_BLACKOIL_EEBOS_HH
