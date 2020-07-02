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
#ifndef EWOMS_MAIN_HH
#define EWOMS_MAIN_HH

#include <eflow/eflow_blackoil.hh>

# ifndef FLOW_BLACKOIL_ONLY
#  include <eflow/eflow_gasoil.hh>
#  include <eflow/eflow_oilwater.hh>
#  include <eflow/eflow_solvent.hh>
#  include <eflow/eflow_polymer.hh>
#  include <eflow/eflow_foam.hh>
#  include <eflow/eflow_brine.hh>
#  include <eflow/eflow_oilwater_brine.hh>
#  include <eflow/eflow_energy.hh>
#  include <eflow/eflow_oilwater_polymer.hh>
#  include <eflow/eflow_oilwater_polymer_injectivity.hh>
# endif

#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/opmlog/eclipseprtlog.hh>
#include <ewoms/eclio/opmlog/logutil.hh>

#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/checkdeck.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/arraydimchecker.hh>

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <ewoms/eclsimulators/eflow/eflowmain.hh>
#include <ewoms/eclsimulators/eflow/missingfeatures.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

#if HAVE_MPI
#include <ewoms/eclsimulators/utils/paralleleclipsestate.hh>
#include <ewoms/eclsimulators/utils/parallelserialization.hh>
#endif

#include <string>
#include <type_traits>

BEGIN_PROPERTIES

// this is a dummy type tag that is used to setup the parameters before the actual
// simulator.
NEW_TYPE_TAG(EFlowEarlyBird, INHERITS_FROM(EclEFlowProblem));

END_PROPERTIES

namespace Ewoms {
  template <class TypeTag>
  void eflowSetDeck(Deck *deck, EclipseState& eclState, Schedule& schedule, SummaryConfig& summaryConfig)
  {
    using Vanguard = GET_PROP_TYPE(TypeTag, Vanguard);
    Vanguard::setExternalDeck(deck);
    Vanguard::setExternalEclState(&eclState);
    Vanguard::setExternalSchedule(&schedule);
    Vanguard::setExternalSummaryConfig(&summaryConfig);
  }

// ----------------- Main program -----------------
  template <class TypeTag>
  int eflowMain(int argc, char** argv, bool outputCout, bool outputFiles)
  {
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    Ewoms::resetLocale();

# if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
# else
    Dune::MPIHelper::instance(argc, argv);
# endif
    Ewoms::EFlowMain<TypeTag> mainfunc;
    return mainfunc.execute(argc, argv, outputCout, outputFiles);
  }
}

namespace Ewoms
{
    // ----------------- NIH main class for eflow -----------------
    // sane persons would extend FlowMain instead of adding yet another class (FlowMain
    // already is a case of NIH: everything except eFlow uses the eWoms::start()
    // function). Unfortunately the OPM people seemed to have decided to shoot themselfs
    // into their feet...
    //
    //   For now, we will either be instantiated from main() in flow.cpp,
    //   or from a Python pybind11 module..
    // NOTE (March 2020): When used from a pybind11 module, we do not neccessarily
    //   want to run the whole simulation by calling run(), it is also
    //   useful to just run one report step at a time. According to these different
    //   usage scenarios, we refactored the original run() in flow.cpp into this class.
    class EFlowNihMain
    {
    private:
        using EFlowMainType = Ewoms::EFlowMain<TTAG(EclEFlowProblem)>;
        enum class FileOutputMode {
            //! \brief No output to files.
            OUTPUT_NONE = 0,
            //! \brief Output only to log files, no eclipse output.
            OUTPUT_LOG_ONLY = 1,
            //! \brief Output to all files.
            OUTPUT_ALL = 3
        };
    public:
        EFlowNihMain(int argc, char** argv) : argc_(argc), argv_(argv)  {  }

        EFlowNihMain(const std::string &filename)
        {
            deckFilename_.assign(filename);
            flowProgName_.assign("flow");
            argc_ = 2;
            saveArgs_[0] = const_cast<char *>(flowProgName_.c_str());
            saveArgs_[1] = const_cast<char *>(deckFilename_.c_str());
            argv_ = saveArgs_;
        }

        EFlowNihMain(int argc,
             char** argv,
             std::shared_ptr<Ewoms::Deck> deck,
             std::shared_ptr<Ewoms::EclipseState> eclipseState,
             std::shared_ptr<Ewoms::Schedule> schedule,
             std::shared_ptr<Ewoms::SummaryConfig> summaryConfig)
            : argc_(argc)
            , argv_(argv)
            , deck_(deck)
            , eclipseState_(eclipseState)
            , schedule_(schedule)
            , summaryConfig_(summaryConfig)
        {
        }

        int runDynamic()
        {
            int exitCode = EXIT_SUCCESS;
            if (initialize_(exitCode)) {
                return dispatchDynamic_();
            } else {
                return exitCode;
            }
        }

        template <class TypeTag>
        int runStatic()
        {
            int exitCode = EXIT_SUCCESS;
            if (initialize_(exitCode)) {
                return dispatchStatic_<TypeTag>();
            } else {
                return exitCode;
            }
        }

        // To be called from the Python interface code. Only do the
        // initialization and then return a pointer to the EFlowEebosMain
        // object that can later be accessed directly from the Python interface
        // to e.g. advance the simulator one report step
        std::unique_ptr<EFlowMainType> initEFlowEebosBlackoil(int& exitCode)
        {
            exitCode = EXIT_SUCCESS;
            if (initialize_(exitCode)) {
                // TODO: check that this deck really represents a blackoil
                // case. E.g. check that number of phases == 3
                Ewoms::eflowBlackoilSetDeck(
                    setupTime_,
                    deck_.get(),
                    *eclipseState_,
                    *schedule_,
                    *summaryConfig_);
                return Ewoms::eflowBlackoilMainInit(
                    argc_, argv_, outputCout_, outputFiles_);
            } else {
                //NOTE: exitCode was set by initialize_() above;
                return std::unique_ptr<EFlowMainType>(); // nullptr
            }
        }

    private:
        int dispatchDynamic_()
        {
            const auto& phases = eclipseState_->runspec().phases();
            // run the actual simulator
            //
            // TODO: make sure that no illegal combinations like thermal and twophase are
            //       requested.

            if ( false ) {}
#ifndef FLOW_BLACKOIL_ONLY
            // Twophase cases
            else if( phases.size() == 2 ) {
                // oil-gas
                if (phases.active( Ewoms::Phase::GAS )) {
                    Ewoms::eflowGasOilSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Ewoms::eflowGasOilMain(argc_, argv_, outputCout_, outputFiles_);
                }
                // oil-water
                else if ( phases.active( Ewoms::Phase::WATER ) ) {
                    Ewoms::eflowOilWaterSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Ewoms::eflowOilWaterMain(argc_, argv_, outputCout_, outputFiles_);
                }
                else {
                    if (outputCout_)
                        std::cerr << "No suitable configuration found, valid are Twophase (oilwater and oilgas), polymer, solvent, or blackoil" << std::endl;
                    return EXIT_FAILURE;
                }
            }
            // Polymer case
            else if ( phases.active( Ewoms::Phase::POLYMER ) ) {
                if ( !phases.active( Ewoms::Phase::WATER) ) {
                    if (outputCout_)
                        std::cerr << "No valid configuration is found for polymer simulation, valid options include "
                                  << "oilwater + polymer and blackoil + polymer" << std::endl;
                    return EXIT_FAILURE;
                }

                // Need to track the polymer molecular weight
                // for the injectivity study
                if ( phases.active( Ewoms::Phase::POLYMW ) ) {
                    // only oil water two phase for now
                    assert( phases.size() == 4);
                    return Ewoms::eflowOilWaterPolymerInjectivityMain(argc_, argv_, outputCout_, outputFiles_);
                }

                if ( phases.size() == 3 ) { // oil water polymer case
                    Ewoms::eflowOilWaterPolymerSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Ewoms::eflowOilWaterPolymerMain(argc_, argv_, outputCout_, outputFiles_);
                } else {
                    Ewoms::eflowPolymerSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Ewoms::eflowPolymerMain(argc_, argv_, outputCout_, outputFiles_);
                }
            }
            // Foam case
            else if ( phases.active( Ewoms::Phase::FOAM ) ) {
                Ewoms::eflowFoamSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                return Ewoms::eflowFoamMain(argc_, argv_, outputCout_, outputFiles_);
            }
            // Brine case
            else if ( phases.active( Ewoms::Phase::BRINE ) ) {
                if ( !phases.active( Ewoms::Phase::WATER) ) {
                    if (outputCout_)
                        std::cerr << "No valid configuration is found for brine simulation, valid options include "
                                  << "oilwater + brine and blackoil + brine" << std::endl;
                    return EXIT_FAILURE;
                }
                if ( phases.size() == 3 ) { // oil water brine case
                    Ewoms::eflowOilWaterBrineSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Ewoms::eflowOilWaterBrineMain(argc_, argv_, outputCout_, outputFiles_);
                } else {
                    Ewoms::eflowBrineSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                    return Ewoms::eflowBrineMain(argc_, argv_, outputCout_, outputFiles_);
                }
            }
            // Solvent case
            else if ( phases.active( Ewoms::Phase::SOLVENT ) ) {
                Ewoms::eflowSolventSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                return Ewoms::eflowSolventMain(argc_, argv_, outputCout_, outputFiles_);
            }
            // Energy case
            else if (eclipseState_->getSimulationConfig().isThermal()) {
                Ewoms::eflowEnergySetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                return Ewoms::eflowEnergyMain(argc_, argv_, outputCout_, outputFiles_);
            }
#endif // FLOW_BLACKOIL_ONLY
            // Blackoil case
            else if( phases.size() == 3 ) {
                Ewoms::eflowBlackoilSetDeck(setupTime_, deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
                return Ewoms::eflowBlackoilMain(argc_, argv_, outputCout_, outputFiles_);
            }
            else {
                if (outputCout_)
                    std::cerr << "No suitable configuration found, valid are Twophase, polymer, foam, brine, solvent, energy, blackoil." << std::endl;
                return EXIT_FAILURE;
            }
        }

        template <class TypeTag>
        int dispatchStatic_()
        {
            Ewoms::eflowSetDeck<TypeTag>(deck_.get(), *eclipseState_, *schedule_, *summaryConfig_);
            return Ewoms::eflowMain<TypeTag>(argc_, argv_, outputCout_, outputFiles_);
        }

        bool initialize_(int& exitCode)
        {
            Dune::Timer externalSetupTimer;
            externalSetupTimer.start();

            handleVersionCmdLine_(argc_, argv_);
            // MPI setup.
#if HAVE_DUNE_FEM
            Dune::Fem::MPIManager::initialize(argc_, argv_);
            int mpiRank = Dune::Fem::MPIManager::rank();
#else
            // the design of the plain dune MPIHelper class is quite flawed: there is no way to
            // get the instance without having the argc and argv parameters available and it is
            // not possible to determine the MPI rank and size without an instance. (IOW: the
            // rank() and size() methods are supposed to be static.)
            const auto& mpiHelper = Dune::MPIHelper::instance(argc_, argv_);
            int mpiRank = mpiHelper.rank();
#endif

            // we always want to use the default locale, and thus spare us the trouble
            // with incorrect locale settings.
            Ewoms::resetLocale();

            // this is a work-around for a catch 22: we do not know what code path to use without
            // parsing the deck, but we don't know the deck without having access to the
            // parameters and this requires to know the type tag to be used. To solve this, we
            // use a type tag just for parsing the parameters before we instantiate the actual
            // simulator object. (Which parses the parameters again, but since this is done in an
            // identical manner it does not matter.)
            typedef TTAG(EFlowEarlyBird) PreTypeTag;
            typedef GET_PROP_TYPE(PreTypeTag, Problem) PreProblem;

            PreProblem::setBriefDescription("EFlow, an advanced reservoir simulator for ECL-decks provided by the eWoms project.");
            int status = Ewoms::EFlowMain<PreTypeTag>::setupParameters_(argc_, argv_);
            if (status != 0) {
                // if setupParameters_ returns a value smaller than 0, there was no error, but
                // the program should abort. This is the case e.g. for the --help and the
                // --print-properties parameters.
#if HAVE_MPI
                if (status < 0)
                    MPI_Finalize(); // graceful stop for --help or --print-properties command line.
                else
                    MPI_Abort(MPI_COMM_WORLD, status);
#endif
                exitCode = (status > 0) ? status : EXIT_SUCCESS;
                return false; //  Whether to run the simulator
            }

            FileOutputMode outputMode = FileOutputMode::OUTPUT_NONE;
            outputCout_ = false;
            if (mpiRank == 0)
                outputCout_ = EWOMS_GET_PARAM(PreTypeTag, bool, EnableTerminalOutput);

            std::string deckFilename;
            std::string outputDir;
            if ( eclipseState_ ) {
                deckFilename = eclipseState_->getIOConfig().fullBasePath();
                outputDir = eclipseState_->getIOConfig().getOutputDir();
            }
            else {
                deckFilename = EWOMS_GET_PARAM(PreTypeTag, std::string, EclDeckFileName);
            }

            typedef GET_PROP_TYPE(PreTypeTag, Vanguard) PreVanguard;
            try {
                deckFilename = PreVanguard::canonicalDeckPath(deckFilename).string();
            }
            catch (const std::exception& e) {
                if ( mpiRank == 0 ) {
                    std::cerr << "Exception received: " << e.what() << ". Try '--help' for a usage description.\n";
                }
#if HAVE_MPI
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
                exitCode = EXIT_FAILURE;
                return false;
            }
            if (outputCout_) {
                Ewoms::EFlowMain<PreTypeTag>::printBanner();
            }
            // Create Deck and EclipseState.
            try {
                if (outputCout_) {
                    std::cout << "Reading deck file '" << deckFilename << "'\n";
                    std::cout.flush();

                    Ewoms::Parser parser;
                    Ewoms::ParseContext parseContext({{Ewoms::ParseContext::PARSE_RANDOM_SLASH, Ewoms::InputError::IGNORE},
                                                    {Ewoms::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Ewoms::InputError::WARN},
                                                    {Ewoms::ParseContext::SUMMARY_UNKNOWN_WELL, Ewoms::InputError::WARN},
                                                    {Ewoms::ParseContext::SUMMARY_UNKNOWN_GROUP, Ewoms::InputError::WARN}});
                    Ewoms::ErrorGuard errorGuard;
                    if (outputDir.empty())
                        outputDir = EWOMS_GET_PARAM(PreTypeTag, std::string, OutputDir);
                    outputMode = setupLogging_(mpiRank,
                                      deckFilename,
                                      outputDir,
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                                      outputCout_, "STDOUT_LOGGER");

                    if (EWOMS_GET_PARAM(PreTypeTag, bool, EclStrictParsing))
                        parseContext.update( Ewoms::InputError::DELAYED_EXIT1);

                    Ewoms::EFlowMain<PreTypeTag>::printPRTHeader(outputCout_);

                    std::string failureMessage;

                    if (mpiRank == 0) {
                        try
                        {
                            if (!deck_)
                                deck_.reset( new Ewoms::Deck( parser.parseFile(deckFilename , parseContext, errorGuard)));
                            Ewoms::MissingFeatures::checkKeywords(*deck_, parseContext, errorGuard);
                            if ( outputCout_ )
                                Ewoms::checkDeck(*deck_, parser, parseContext, errorGuard);

                            if (!eclipseState_) {
#if HAVE_MPI
                                eclipseState_.reset(new Ewoms::ParallelEclipseState(*deck_));
#else
                                eclipseState_.reset(new Ewoms::EclipseState(*deck_));
#endif
                            }
                            /*
                              For the time being initializing wells and groups from the
                              restart file is not possible, but work is underways and it is
                              included here as a switch.
                            */
                            const bool init_from_restart_file = !EWOMS_GET_PARAM(PreTypeTag, bool, SchedRestart);
                            const auto& init_config = eclipseState_->getInitConfig();
                            if (init_config.restartRequested() && init_from_restart_file) {
                                int report_step = init_config.getRestartStep();
                                const auto& rst_filename = eclipseState_->getIOConfig().getRestartFileName( init_config.getRestartRootName(), report_step, false );
                                Ewoms::EclIO::ERst rst_file(rst_filename);
                                const auto& rst_state = Ewoms::RestartIO::RstState::load(rst_file, report_step);
                                if (!schedule_)
                                    schedule_.reset(new Ewoms::Schedule(*deck_, *eclipseState_, parseContext, errorGuard, &rst_state) );
                            }
                            else {
                                if (!schedule_)
                                    schedule_.reset(new Ewoms::Schedule(*deck_, *eclipseState_, parseContext, errorGuard));
                            }
                            setupMessageLimiter_(schedule_->getMessageLimits(), "STDOUT_LOGGER");
                            if (!summaryConfig_)
                                summaryConfig_.reset( new Ewoms::SummaryConfig(*deck_, *schedule_, eclipseState_->getTableManager(), parseContext, errorGuard));
                        }
                        catch(const std::exception& e)
                        {
                            failureMessage = e.what();
                        }
                    }
#if HAVE_MPI
                    else {
                        if (!summaryConfig_)
                            summaryConfig_.reset(new Ewoms::SummaryConfig);
                        if (!schedule_)
                            schedule_.reset(new Ewoms::Schedule());
                        if (!eclipseState_)
                            eclipseState_.reset(new Ewoms::ParallelEclipseState);
                    }

                    auto comm = Dune::MPIHelper::getCollectiveCommunication();
                    if (false) // hack for smaller delta to OPM
                    {
                        if (errorGuard) {
                            errorGuard.dump();
                            errorGuard.clear();
                        }
                        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
                    }

                    Ewoms::eclStateBroadcast(*eclipseState_, *schedule_, *summaryConfig_);
#endif

                    Ewoms::checkConsistentArrayDimensions(*eclipseState_, *schedule_, parseContext, errorGuard);

                    if (errorGuard) {
                        errorGuard.dump();
                        errorGuard.clear();

                        throw std::runtime_error("Unrecoverable errors were encountered while loading input.");
                    }
                }
                setupTime_ = externalSetupTimer.elapsed();
                outputFiles_ = (outputMode != FileOutputMode::OUTPUT_NONE);
            }
            catch (const std::invalid_argument& e)
            {
                if (outputCout_) {
                    std::cerr << "Failed to create valid EclipseState object." << std::endl;
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                }
#if HAVE_MPI
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
                exitCode = EXIT_FAILURE;
                return false;
            }

            exitCode = EXIT_SUCCESS;
            return true;
        }

        Ewoms::filesystem::path simulationCaseName_( const std::string& casename ) {
            namespace fs = Ewoms::filesystem;

            const auto exists = []( const fs::path& f ) -> bool {
                if( !fs::exists( f ) ) return false;

                if( fs::is_regular_file( f ) ) return true;

                return fs::is_symlink( f )
                && fs::is_regular_file( fs::read_symlink( f ) );
            };

            auto simcase = fs::path( casename );

            if( exists( simcase ) ) {
                return simcase;
            }

            for( const auto& ext : { std::string("data"), std::string("DATA") } ) {
                if( exists( simcase.replace_extension( ext ) ) ) {
                    return simcase;
                }
            }

            throw std::invalid_argument( "Cannot find input case " + casename );
        }

        // This function is an extreme special case, if the program has been invoked
        // *exactly* as:
        //
        //    eflow   --version
        //
        // the call is intercepted by this function which will print "eflow $version"
        // on stdout and exit(0).
        void handleVersionCmdLine_(int argc, char** argv) {
            for ( int i = 1; i < argc; ++i )
            {
                if (std::strcmp(argv[i], "--version") == 0) {
                    std::cout << "eflow " << EWOMS_ECLSIMULATORS_VERSION << std::endl;
                    std::exit(EXIT_SUCCESS);
                }
            }
        }

        void ensureOutputDirExists_(const std::string& cmdline_output_dir)
        {
            if (!Ewoms::filesystem::is_directory(cmdline_output_dir)) {
                try {
                    Ewoms::filesystem::create_directories(cmdline_output_dir);
                }
                catch (...) {
                    throw std::runtime_error("Creation of output directory '" + cmdline_output_dir + "' failed\n");
                }
            }
        }

        // Setup the OpmLog backends
        FileOutputMode setupLogging_(int mpi_rank_, const std::string& deck_filename, const std::string& cmdline_output_dir, const std::string& cmdline_output, bool output_cout_, const std::string& stdout_log_id) {

            if (!cmdline_output_dir.empty()) {
                ensureOutputDirExists_(cmdline_output_dir);
            }

            // create logFile
            using Ewoms::filesystem::path;
            path fpath(deck_filename);
            std::string baseName;
            std::ostringstream debugFileStream;
            std::ostringstream logFileStream;

            // Strip extension "." or ".DATA"
            std::string extension = boost::to_upper_copy(fpath.extension().string());
            if (extension == ".DATA" || extension == ".") {
                baseName = boost::to_upper_copy(fpath.stem().string());
            } else {
                baseName = boost::to_upper_copy(fpath.filename().string());
            }

            std::string output_dir = cmdline_output_dir;
            if (output_dir.empty()) {
                output_dir = fpath.has_parent_path()
                    ? absolute(fpath.parent_path()).generic_string()
                    : Ewoms::filesystem::current_path().generic_string();
            }

            logFileStream << output_dir << "/" << baseName;
            debugFileStream << output_dir << "/" << baseName;

            if (mpi_rank_ != 0) {
                // Added rank to log file for non-zero ranks.
                // This prevents message loss.
                debugFileStream << "." << mpi_rank_;
                // If the following file appears then there is a bug.
                logFileStream << "." << mpi_rank_;
            }
            logFileStream << ".PRT";
            debugFileStream << ".DBG";

            FileOutputMode output;
            {
                static std::map<std::string, FileOutputMode> stringToOutputMode =
                    { {"none", FileOutputMode::OUTPUT_NONE },
                      {"false", FileOutputMode::OUTPUT_LOG_ONLY },
                      {"log", FileOutputMode::OUTPUT_LOG_ONLY },
                      {"all" , FileOutputMode::OUTPUT_ALL },
                      {"true" , FileOutputMode::OUTPUT_ALL }};
                auto outputModeIt = stringToOutputMode.find(cmdline_output);
                if (outputModeIt != stringToOutputMode.end()) {
                    output = outputModeIt->second;
                }
                else {
                    output = FileOutputMode::OUTPUT_ALL;
                    std::cerr << "Value " << cmdline_output <<
                        " is not a recognized output mode. Using \"all\" instead."
                              << std::endl;
                }
            }

            if (output > FileOutputMode::OUTPUT_NONE) {
                std::shared_ptr<Ewoms::EclipsePRTLog> prtLog = std::make_shared<Ewoms::EclipsePRTLog>(logFileStream.str(), Ewoms::Log::NoDebugMessageTypes, false, output_cout_);
                Ewoms::OpmLog::addBackend("ECLIPSEPRTLOG", prtLog);
                prtLog->setMessageLimiter(std::make_shared<Ewoms::MessageLimiter>());
                prtLog->setMessageFormatter(std::make_shared<Ewoms::SimpleMessageFormatter>(false));
            }

            if (output >= FileOutputMode::OUTPUT_LOG_ONLY) {
                std::string debugFile = debugFileStream.str();
                std::shared_ptr<Ewoms::StreamLog> debugLog = std::make_shared<Ewoms::EclipsePRTLog>(debugFileStream.str(), Ewoms::Log::DefaultMessageTypes, false, output_cout_);
                Ewoms::OpmLog::addBackend("DEBUGLOG", debugLog);
            }

            if (mpi_rank_ == 0) {
                std::shared_ptr<Ewoms::StreamLog> streamLog = std::make_shared<Ewoms::StreamLog>(std::cout, Ewoms::Log::StdoutMessageTypes);
                Ewoms::OpmLog::addBackend(stdout_log_id, streamLog);
                streamLog->setMessageFormatter(std::make_shared<Ewoms::SimpleMessageFormatter>(true));
            }
            return output;
        }

        void setupMessageLimiter_(const Ewoms::MessageLimits msgLimits,  const std::string& stdout_log_id) {
            std::shared_ptr<Ewoms::StreamLog> stream_log = Ewoms::OpmLog::getBackend<Ewoms::StreamLog>(stdout_log_id);

            const std::map<int64_t, int> limits = {{Ewoms::Log::MessageType::Note,
                                            msgLimits.getCommentPrintLimit(0)},
                                           {Ewoms::Log::MessageType::Info,
                                            msgLimits.getMessagePrintLimit(0)},
                                           {Ewoms::Log::MessageType::Warning,
                                            msgLimits.getWarningPrintLimit(0)},
                                           {Ewoms::Log::MessageType::Error,
                                            msgLimits.getErrorPrintLimit(0)},
                                           {Ewoms::Log::MessageType::Problem,
                                            msgLimits.getProblemPrintLimit(0)},
                                           {Ewoms::Log::MessageType::Bug,
                                            msgLimits.getBugPrintLimit(0)}};
            stream_log->setMessageLimiter(std::make_shared<Ewoms::MessageLimiter>(10, limits));
        }

        int argc_;
        char** argv_;
        bool outputCout_;
        bool outputFiles_;
        double setupTime_;
        std::string deckFilename_;
        std::string flowProgName_;
        char *saveArgs_[2];
        std::shared_ptr<Ewoms::Deck> deck_;
        std::shared_ptr<Ewoms::EclipseState> eclipseState_;
        std::shared_ptr<Ewoms::Schedule> schedule_;
        std::shared_ptr<Ewoms::SummaryConfig> summaryConfig_;
    };

} // namespace Ewoms

#endif // EWOMS_MAIN_HH
