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
#ifndef EWOMS_EFLOW_MAIN_EEBOS_HH
#define EWOMS_EFLOW_MAIN_EEBOS_HH

#include <sys/utsname.h>

#include <ewoms/eclsimulators/eflow/blackoilmodel.hh>
#include <ewoms/eclsimulators/eflow/simulatorfullyimplicitblackoil.hh>
#include <ewoms/eclsimulators/utils/parallelfilemerger.hh>
#include <ewoms/eclsimulators/linalg/extractparallelgridinformationtoistl.hh>

#include <ewoms/eclsimulators/deprecated/props/satfunc/relpermdiagnostics.hh>

#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/parsecontext.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/ioconfig/ioconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/initconfig/initconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/checkdeck.hh>
#include <ewoms/common/string.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

BEGIN_PROPERTIES

NEW_PROP_TAG(EnableDryRun);
NEW_PROP_TAG(OutputInterval);
NEW_PROP_TAG(UseAmg);
NEW_PROP_TAG(EnableLoggingFalloutWarning);

// TODO: enumeration parameters. we use strings for now.
SET_STRING_PROP(EclEFlowProblem, EnableDryRun, "auto");
// Do not merge parallel output files or warn about them
SET_BOOL_PROP(EclEFlowProblem, EnableLoggingFalloutWarning, false);
SET_INT_PROP(EclEFlowProblem, OutputInterval, 1);

END_PROPERTIES

namespace Ewoms
{
    // The EFlowMain class is the eebos based black-oil simulator.
    template <class TypeTag>
    class EFlowMain
    {
    public:
        using MaterialLawManager = typename GET_PROP(TypeTag, MaterialLaw)::EclMaterialLawManager;
        using EebosSimulator = GET_PROP_TYPE(TypeTag, Simulator);
        using Grid = GET_PROP_TYPE(TypeTag, Grid);
        using GridView = GET_PROP_TYPE(TypeTag, GridView);
        using Problem = GET_PROP_TYPE(TypeTag, Problem);
        using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
        using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);

        typedef Ewoms::SimulatorFullyImplicitBlackoil<TypeTag> Simulator;

        EFlowMain(int argc, char **argv, bool output_cout, bool output_files )
            : argc_{argc}, argv_{argv},
              output_cout_{output_cout}, output_files_{output_files}
        {

        }

        // Read the command line parameters. Throws an exception if something goes wrong.
        static int setupParameters_(int argc, char** argv)
        {
            using ParamsMeta = GET_PROP(TypeTag, ParameterMetaData);
            if (!ParamsMeta::registrationOpen()) {
                // We have already successfully run setupParameters_().
                // For the dynamically chosen runs (as from the main flow
                // executable) we must run this function again with the
                // real typetag to be used, as the first time was with the
                // "EFlowEarlyBird" typetag. However, for the static ones (such
                // as 'flow_onephase_energy') it has already been run with the
                // correct typetag.
                return EXIT_SUCCESS;
            }
            // register the eflow specific parameters
            EWOMS_REGISTER_PARAM(TypeTag, std::string, EnableDryRun,
                                 "Specify if the simulation ought to be actually run, or just pretended to be");
            EWOMS_REGISTER_PARAM(TypeTag, int, OutputInterval,
                                 "Specify the number of report steps between two consecutive writes of restart data");
            EWOMS_REGISTER_PARAM(TypeTag, bool, EnableLoggingFalloutWarning,
                                 "Developer option to see whether logging was on non-root processors. In that case it will be appended to the *.DBG or *.PRT files");

            Simulator::registerParameters();

            // register the parameters inherited from eebos
            Ewoms::registerAllParameters_<TypeTag>(/*finalizeRegistration=*/false);

            // hide the parameters unused by eflow. TODO: this is a pain to maintain
            EWOMS_HIDE_PARAM(TypeTag, EnableGravity);
            EWOMS_HIDE_PARAM(TypeTag, EnableGridAdaptation);

            // this parameter is actually used in eWoms, but the eflow well model
            // hard-codes the assumption that the intensive quantities cache is enabled,
            // so eflow crashes. Let's hide the parameter for that reason.
            EWOMS_HIDE_PARAM(TypeTag, EnableIntensiveQuantityCache);

            // thermodynamic hints are not implemented/required by the eWoms blackoil
            // model
            EWOMS_HIDE_PARAM(TypeTag, EnableThermodynamicHints);

            // in eflow, only the deck file determines the end time of the simulation
            EWOMS_HIDE_PARAM(TypeTag, EndTime);

            // time stepping is not done by the eWoms code in eflow
            EWOMS_HIDE_PARAM(TypeTag, InitialTimeStepSize);
            EWOMS_HIDE_PARAM(TypeTag, MaxTimeStepDivisions);
            EWOMS_HIDE_PARAM(TypeTag, MaxTimeStepSize);
            EWOMS_HIDE_PARAM(TypeTag, MinTimeStepSize);
            EWOMS_HIDE_PARAM(TypeTag, PredeterminedTimeStepsFile);

            EWOMS_HIDE_PARAM(TypeTag, EclMaxTimeStepSizeAfterWellEvent);
            EWOMS_HIDE_PARAM(TypeTag, EclRestartShrinkFactor);
            EWOMS_HIDE_PARAM(TypeTag, EclEnableTuning);

            // eflow also does not use the eWoms Newton method
            EWOMS_HIDE_PARAM(TypeTag, NewtonMaxError);
            EWOMS_HIDE_PARAM(TypeTag, NewtonMaxIterations);
            EWOMS_HIDE_PARAM(TypeTag, NewtonTolerance);
            EWOMS_HIDE_PARAM(TypeTag, NewtonTargetIterations);
            EWOMS_HIDE_PARAM(TypeTag, NewtonVerbose);
            EWOMS_HIDE_PARAM(TypeTag, NewtonWriteConvergence);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonSumTolerance);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonSumToleranceExponent);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonStrictIterations);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonRelaxedVolumeFraction);
            EWOMS_HIDE_PARAM(TypeTag, EclNewtonRelaxedTolerance);

            // the default eWoms checkpoint/restart mechanism does not work with eflow
            EWOMS_HIDE_PARAM(TypeTag, RestartTime);
            EWOMS_HIDE_PARAM(TypeTag, RestartWritingInterval);

            EWOMS_END_PARAM_REGISTRATION(TypeTag);

            int mpiRank = 0;
#if HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif

            // read in the command line parameters
            int status = Ewoms::setupParameters_<TypeTag>(argc, const_cast<const char**>(argv), /*doRegistration=*/false, /*allowUnused=*/true, /*handleHelp=*/(mpiRank==0));
            if (status == 0) {

                // deal with unknown parameters.

                int unknownKeyWords = 0;
                if (mpiRank == 0) {
                    unknownKeyWords = Ewoms::Parameters::printUnused<TypeTag>(std::cerr);
                }
#if HAVE_MPI
                int globalUnknownKeyWords;
                MPI_Allreduce(&unknownKeyWords,  &globalUnknownKeyWords, 1, MPI_INT,  MPI_SUM, MPI_COMM_WORLD);
                unknownKeyWords = globalUnknownKeyWords;
#endif
                if ( unknownKeyWords )
                {
                    if ( mpiRank == 0 )
                    {
                        std::string msg = "Aborting simulation due to unknown "
                            "parameters. Please query \"flow --help\" for "
                            "supported command line parameters.";
                        if (OpmLog::hasBackend("STREAMLOG"))
                        {
                            OpmLog::error(msg);
                        }
                        else {
                            std::cerr << msg << std::endl;
                        }
                    }
                    return EXIT_FAILURE;
                }

                // deal with --print-properties and --print-parameters and unknown parameters.

                bool doExit = false;

                if (EWOMS_GET_PARAM(TypeTag, int, PrintProperties) == 1) {
                    doExit = true;
                    if (mpiRank == 0)
                        EWOMS_PRINT_PROPERTIES(TypeTag);
                }

                if (EWOMS_GET_PARAM(TypeTag, int, PrintParameters) == 1) {
                    doExit = true;
                    if (mpiRank == 0)
                        Ewoms::Parameters::printValues<TypeTag>();
                }

                if (doExit)
                    return -1;
            }

            return status;
        }

        static void printBanner()
        {
            const int lineLen = 71;
            const std::string version = EWOMS_ECLSIMULATORS_VERSION;
            const std::string banner = "This is eflow "+version;
            const int bannerPreLen = (lineLen - 2 - banner.size())/2;
            const int bannerPostLen = bannerPreLen + (lineLen - 2 - banner.size())%2;
            std::cout << "***********************************************************************\n";
            std::cout << "*                                                                     *\n";
            std::cout << "*" << std::string(bannerPreLen, ' ') << banner << std::string(bannerPostLen, ' ') << "*\n";
            std::cout << "*                                                                     *\n";
            std::cout << "* EFlow is a simulator for fully implicit three-phase black-oil flow, *\n";
            std::cout << "*             including solvent and polymer capabilities.             *\n";
            std::cout << "*            For more information, visit https://ewoms.org            *\n";
            std::cout << "*                                                                     *\n";
            std::cout << "*  EFLOW IS NOT ACTIVELY MAINTAINED WITHIN eWoms. ONLY USE IT IF YOU  *\n";
            std::cout << "*       KNOW EXACTLY WHAT YOU ARE DOING AND WHAT IT IS GOOD FOR!      *\n";
            std::cout << "***********************************************************************\n\n";

            int mpiSize = 1;

#if HAVE_MPI
            MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
#endif

            std::cout << "Using "<< mpiSize << " MPI processes\n\n";
        }

        /// This is the main function of EFlow.  It runs a complete simulation with the
        /// given grid and simulator classes, based on the user-specified command-line
        /// input.
        int execute()
        {
            return execute_(&EFlowMain::runSimulator, /*cleanup=*/true);
        }

        int executeInitStep()
        {
            return execute_(&EFlowMain::runSimulatorInit, /*cleanup=*/false);
        }

        // Print an ASCII-art header to the PRT and DEBUG files.
        // \return Whether unkown keywords were seen during parsing.
        static void printPRTHeader(bool output_cout)
        {
          if (output_cout) {
              const std::string version = EWOMS_ECLSIMULATORS_VERSION;
              const double megabyte = 1024 * 1024;
              unsigned num_cpu = std::thread::hardware_concurrency();
              struct utsname arch;
              const char* user = getlogin();
              time_t now = std::time(0);
              struct tm  tstruct;
              char      tmstr[80];
              tstruct = *localtime(&now);
              strftime(tmstr, sizeof(tmstr), "%d-%m-%Y at %X", &tstruct);
              const double mem_size = getTotalSystemMemory() / megabyte;
              std::ostringstream ss;
              ss << "\n\n\n";
              ss << " ########   ########  #          ######   #           #\n";
              ss << " #          #         #         #      #   #         # \n";
              ss << " ######     #####     #         #      #    #   #   #  \n";
              ss << " #          #         #         #      #     # # # #   \n";
              ss << " ########   #         #######    ######       #   #    \n\n";
              ss << "EFlow is a simulator for fully implicit three-phase black-oil flow,";
              ss << "For more information, visit https://ewoms.org\n\n";
              ss << "\n";
              ss << "EFLOW IS NOT ACTIVELY MAINTAINED WITHIN eWoms. ONLY USE IT IF YOU KNOW EXACTLY\n";
              ss << "WHAT YOU ARE DOING AND WHAT IT IS GOOD FOR!\n";
              ss << "\n";
              ss << "EFlow Version     =  " + version + "\n";
              if (uname(&arch) == 0) {
                 ss << "Machine name     =  " << arch.nodename << " (Number of logical cores: " << num_cpu;
                 ss << ", Memory size: " << std::fixed << std::setprecision (2) << mem_size << " MB) \n";
                 ss << "Operating system =  " << arch.sysname << " " << arch.machine << " (Kernel: " << arch.release;
                 ss << ", " << arch.version << " )\n";
                 }
              if (user) {
                 ss << "User             =  " << user << std::endl;
                 }
              ss << "Simulation started on " << tmstr << " hrs\n";

              ss << "Parameters used by EFlow:\n";
              Ewoms::Parameters::printValues<TypeTag>(ss);

              OpmLog::note(ss.str());
          }
        }
    private:
        // called by execute() or executeInitStep()
        int execute_(int (EFlowMain::* runOrInitFunc)(), bool cleanup)
        {
            try {
                // deal with some administrative boilerplate

                int status = setupParameters_(this->argc_, this->argv_);
                if (status)
                    return status;

                setupParallelism();
                setupEebosSimulator();
                runDiagnostics();
                createSimulator();

                // if run, do the actual work, else just initialize
                int exitCode = (this->*runOrInitFunc)();
                if (cleanup) {
                    executeCleanup_();
                }
                return exitCode;
            }
            catch (const std::exception& e) {
                std::ostringstream message;
                message  << "Program threw an exception: " << e.what();

                if (this->output_cout_) {
                    // in some cases exceptions are thrown before the logging system is set
                    // up.
                    if (OpmLog::hasBackend("STREAMLOG")) {
                        OpmLog::error(message.str());
                    }
                    else {
                        std::cout << message.str() << "\n";
                    }
                }

                return EXIT_FAILURE;
            }
        }

        void executeCleanup_() {
            // clean up
            mergeParallelLogFiles();
        }

    protected:
        void setupParallelism()
        {
            // determine the rank of the current process and the number of processes
            // involved in the simulation. MPI must have already been initialized
            // here. (yes, the name of this method is misleading.)
#if HAVE_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
#else
            mpi_rank_ = 0;
            mpi_size_ = 1;
#endif
        }

        void mergeParallelLogFiles()
        {
            // force closing of all log files.
            OpmLog::removeAllBackends();

            if (mpi_rank_ != 0 || mpi_size_ < 2 || !this->output_files_) {
                return;
            }

            namespace fs = Ewoms::filesystem;
            const std::string& output_dir = eclState().getIOConfig().getOutputDir();
            fs::path output_path(output_dir);
            fs::path deck_filename(EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName));
            std::string basename;
            // Strip extension "." and ".DATA"
            std::string extension = uppercase(deck_filename.extension().string());
            if ( extension == ".DATA" || extension == "." )
            {
                basename = uppercase(deck_filename.stem().string());
            }
            else
            {
                basename = uppercase(deck_filename.filename().string());
            }
            std::for_each(fs::directory_iterator(output_path),
                          fs::directory_iterator(),
                          detail::ParallelFileMerger(output_path, basename,
                                                     EWOMS_GET_PARAM(TypeTag, bool, EnableLoggingFalloutWarning)));
        }

        void setupEebosSimulator()
        {
            eebosSimulator_.reset(new EebosSimulator(/*verbose=*/false));
            eebosSimulator_->executionTimer().start();
            eebosSimulator_->model().applyInitialSolution();

            try {
                // Possible to force initialization only behavior (NOSIM).
                const std::string& dryRunString = EWOMS_GET_PARAM(TypeTag, std::string, EnableDryRun);
                if (dryRunString != "" && dryRunString != "auto") {
                    bool yesno;
                    if (dryRunString == "true"
                        || dryRunString == "t"
                        || dryRunString == "1")
                        yesno = true;
                    else if (dryRunString == "false"
                             || dryRunString == "f"
                             || dryRunString == "0")
                        yesno = false;
                    else
                        throw std::invalid_argument("Invalid value for parameter EnableDryRun: '"
                                                    +dryRunString+"'");
                    auto& ioConfig = eclState().getIOConfig();
                    ioConfig.overrideNOSIM(yesno);
                }
            }
            catch (const std::invalid_argument& e) {
                std::cerr << "Failed to create valid EclipseState object" << std::endl;
                std::cerr << "Exception caught: " << e.what() << std::endl;
                throw;
            }
        }

        const Deck& deck() const
        { return eebosSimulator_->vanguard().deck(); }

        Deck& deck()
        { return eebosSimulator_->vanguard().deck(); }

        const EclipseState& eclState() const
        { return eebosSimulator_->vanguard().eclState(); }

        EclipseState& eclState()
        { return eebosSimulator_->vanguard().eclState(); }

        const Schedule& schedule() const
        { return eebosSimulator_->vanguard().schedule(); }

        // Run diagnostics.
        // Writes to:
        //   OpmLog singleton.
        void runDiagnostics()
        {
            if (!this->output_cout_) {
                return;
            }

            // Run relperm diagnostics if we have more than one phase.
            if (FluidSystem::numActivePhases() > 1) {
                RelpermDiagnostics diagnostic;
                if (mpi_size_ > 1) {
#if HAVE_MPI
                    this->grid().switchToGlobalView();
                    static_cast<ParallelEclipseState&>(this->eclState()).switchToGlobalProps();
#endif
                }
                diagnostic.diagnosis(eclState(), this->grid());
                if (mpi_size_ > 1) {
#if HAVE_MPI
                    this->grid().switchToDistributedView();
                    static_cast<ParallelEclipseState&>(this->eclState()).switchToDistributedProps();
#endif
                }
            }
        }

        // Run the simulator.
        int runSimulator()
        {
            return runSimulatorInitOrRun_(&EFlowMain::runSimulatorRunCallback_);
        }

        int runSimulatorInit()
        {
            return runSimulatorInitOrRun_(&EFlowMain::runSimulatorInitCallback_);
        }

    private:
        // Callback that will be called from runSimulatorInitOrRun_().
        int runSimulatorRunCallback_()
        {
            SimulatorReport report = simulator_->run(*simtimer_);
            runSimulatorAfterSim_(report);
            return report.success.exit_status;
        }

        // Callback that will be called from runSimulatorInitOrRun_().
        int runSimulatorInitCallback_()
        {
            simulator_->init(*simtimer_);
            return EXIT_SUCCESS;
        }

        // Output summary after simulation has completed
        void runSimulatorAfterSim_(SimulatorReport &report)
        {
            if (this->output_cout_) {
                std::ostringstream ss;
                ss << "\n\n================    End of simulation     ===============\n\n";
                ss << "Number of MPI processes: " << std::setw(6) << mpi_size_ << "\n";
#if _OPENMP
                int threads = omp_get_max_threads();
#else
                int threads = 1;
#endif
                ss << "Threads per MPI process:  " << std::setw(5) << threads << "\n";
                report.reportFullyImplicit(ss);
                OpmLog::info(ss.str());
                const std::string dir = eclState().getIOConfig().getOutputDir();
                namespace fs = Ewoms::filesystem;
                fs::path output_dir(dir);
                {
                    std::string filename = eclState().getIOConfig().getBaseName() + ".INFOSTEP";
                    fs::path fullpath = output_dir / filename;
                    std::ofstream os(fullpath.string());
                    report.fullReports(os);
                }
            }
        }

        // Run the simulator.
        int runSimulatorInitOrRun_(int (EFlowMain::* initOrRunFunc)())
        {

            const auto& schedule = this->schedule();
            const auto& timeMap = schedule.getTimeMap();
            auto& ioConfig = eclState().getIOConfig();
            simtimer_ = std::make_unique<SimulatorTimer>();

            // initialize variables
            const auto& initConfig = eclState().getInitConfig();
            simtimer_->init(timeMap, (size_t)initConfig.getRestartStep());

            if (this->output_cout_) {
                std::ostringstream oss;

                // This allows a user to catch typos and misunderstandings in the
                // use of simulator parameters.
                if (Ewoms::Parameters::printUnused<TypeTag>(oss)) {
                    std::cout << "-----------------   Unrecognized parameters:   -----------------\n";
                    std::cout << oss.str();
                    std::cout << "----------------------------------------------------------------" << std::endl;
                }
            }

            if (!ioConfig.initOnly()) {
                if (this->output_cout_) {
                    std::string msg;
                    msg = "\n\n================ Starting main simulation loop ===============\n";
                    OpmLog::info(msg);
                }

                return (this->*initOrRunFunc)();
            }
            else {
                if (this->output_cout_) {
                    std::cout << "\n\n================ Simulation turned off ===============\n" << std::flush;
                }
                return EXIT_SUCCESS;
            }
        }

    protected:

        /// This is the main function of EFlow.
        // Create simulator instance.
        // Writes to:
        //   simulator_
        void createSimulator()
        {
            // Create the simulator instance.
            simulator_.reset(new Simulator(*eebosSimulator_));
        }

        static unsigned long long getTotalSystemMemory()
        {
            long pages = sysconf(_SC_PHYS_PAGES);
            long page_size = sysconf(_SC_PAGE_SIZE);
            return pages * page_size;
        }

        Grid& grid()
        { return eebosSimulator_->vanguard().grid(); }

    private:
        std::unique_ptr<EebosSimulator> eebosSimulator_;
        int  mpi_rank_ = 0;
        int  mpi_size_ = 1;
        Ewoms::any parallel_information_;
        std::unique_ptr<Simulator> simulator_;
        std::unique_ptr<SimulatorTimer> simtimer_;
        int argc_;
        char **argv_;
        bool output_cout_;
        bool output_files_;
    };
} // namespace Ewoms

#endif // EWOMS_EFLOW_MAIN_EEBOS_HH
