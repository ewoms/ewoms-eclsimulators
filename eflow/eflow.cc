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

#include <eebos/eclmpiserializer.hh>

#include <eflow/eflow_blackoil.hh>

#ifndef FLOW_BLACKOIL_ONLY
#include <eflow/eflow_gasoil.hh>
#include <eflow/eflow_oilwater.hh>
#include <eflow/eflow_solvent.hh>
#include <eflow/eflow_polymer.hh>
#include <eflow/eflow_foam.hh>
#include <eflow/eflow_brine.hh>
#include <eflow/eflow_energy.hh>
#include <eflow/eflow_oilwater_polymer.hh>
#include <eflow/eflow_oilwater_polymer_injectivity.hh>
#endif

#include <ewoms/eclsimulators/eflow/simulatorfullyimplicitblackoil.hh>
#include <ewoms/eclsimulators/eflow/eflowmain.hh>
#include <ewoms/eclsimulators/utils/paralleleclipsestate.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>
#include <ewoms/eclsimulators/eflow/missingfeatures.hh>
#include <ewoms/common/resetlocale.hh>

#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/opmlog/eclipseprtlog.hh>
#include <ewoms/eclio/opmlog/logutil.hh>

#include <ewoms/eclio/io/rst/state.hh>
#include <ewoms/eclio/io/erst.hh>

#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/parsecontext.hh>
#include <ewoms/eclio/parser/errorguard.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/checkdeck.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/summaryconfig/summaryconfig.hh>

#include <ewoms/eclio/parser/eclipsestate/schedule/arraydimchecker.hh>

#include <ewoms/common/filesystem.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

BEGIN_PROPERTIES

// this is a dummy type tag that is used to setup the parameters before the actual
// simulator.
NEW_TYPE_TAG(EFlowEarlyBird, INHERITS_FROM(EclEFlowProblem));

END_PROPERTIES

namespace detail
{
    Ewoms::filesystem::path simulationCaseName( const std::string& casename ) {
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
    void handleVersionCmdLine(int argc, char** argv) {
        for ( int i = 1; i < argc; ++i )
        {
            if (std::strcmp(argv[i], "--version") == 0) {
                std::cout << "eflow " << EWOMS_ECLSIMULATORS_VERSION << std::endl;
                std::exit(EXIT_SUCCESS);
            }
        }
    }

}

enum class FileOutputMode {
    //! \brief No output to files.
    OUTPUT_NONE = 0,
    //! \brief Output only to log files, no eclipse output.
    OUTPUT_LOG_ONLY = 1,
    //! \brief Output to all files.
    OUTPUT_ALL = 3
};

void ensureOutputDirExists(const std::string& cmdline_output_dir)
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
FileOutputMode setupLogging(int mpi_rank_, const std::string& deck_filename, const std::string& cmdline_output_dir, const std::string& cmdline_output, bool output_cout_, const std::string& stdout_log_id) {

    if (!cmdline_output_dir.empty()) {
        ensureOutputDirExists(cmdline_output_dir);
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
        output_dir = absolute(path(baseName).parent_path()).string();
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

void setupMessageLimiter(const Ewoms::MessageLimits msgLimits,  const std::string& stdout_log_id) {
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

// ----------------- Main program -----------------
int main(int argc, char** argv)
{
    Dune::Timer externalSetupTimer;
    externalSetupTimer.start();

    detail::handleVersionCmdLine(argc, argv);
    // MPI setup.
#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
    int mpiRank = Dune::Fem::MPIManager::rank();
#else
    // the design of the plain dune MPIHelper class is quite flawed: there is no way to
    // get the instance without having the argc and argv parameters available and it is
    // not possible to determine the MPI rank and size without an instance. (IOW: the
    // rank() and size() methods are supposed to be static.)
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);
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
    int status = Ewoms::EFlowMain<PreTypeTag>::setupParameters_(argc, argv);
    if (status != 0) {
        // if setupParameters_ returns a value smaller than 0, there was no error, but
        // the program should abort. This is the case e.g. for the --help and the
        // --print-properties parameters.
#if HAVE_MPI
        MPI_Finalize();
#endif
        return (status >= 0)?status:0;
    }

    FileOutputMode outputMode = FileOutputMode::OUTPUT_NONE;
    bool outputCout = false;
    if (mpiRank == 0)
        outputCout = EWOMS_GET_PARAM(PreTypeTag, bool, EnableTerminalOutput);

    std::string deckFilename = EWOMS_GET_PARAM(PreTypeTag, std::string, EclDeckFileName);
    typedef typename GET_PROP_TYPE(PreTypeTag, Vanguard) PreVanguard;
    try {
        deckFilename = PreVanguard::canonicalDeckPath(deckFilename).string();
    }
    catch (const std::exception& e) {
        if ( mpiRank == 0 )
            std::cerr << "Exception received: " << e.what() << ". Try '--help' for a usage description.\n";
#if HAVE_MPI
        MPI_Finalize();
#endif
        return 1;
    }

    if (outputCout) {
        Ewoms::EFlowMain<PreTypeTag>::printBanner();
    }

    // Create Deck and EclipseState.
    try {
        if (outputCout) {
            std::cout << "Reading deck file '" << deckFilename << "'\n";
            std::cout.flush();
        }
        std::shared_ptr<Ewoms::Deck> deck;
        std::shared_ptr<Ewoms::EclipseState> eclipseState;
        std::shared_ptr<Ewoms::Schedule> schedule;
        std::shared_ptr<Ewoms::SummaryConfig> summaryConfig;
        {
            Ewoms::Parser parser;
            Ewoms::ParseContext parseContext({{Ewoms::ParseContext::PARSE_RANDOM_SLASH, Ewoms::InputError::IGNORE},
                                            {Ewoms::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Ewoms::InputError::WARN},
                                            {Ewoms::ParseContext::SUMMARY_UNKNOWN_WELL, Ewoms::InputError::WARN},
                                            {Ewoms::ParseContext::SUMMARY_UNKNOWN_GROUP, Ewoms::InputError::WARN}});
            Ewoms::ErrorGuard errorGuard;
            outputMode = setupLogging(mpiRank,
                                      deckFilename,
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputDir),
                                      EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                                      outputCout, "STDOUT_LOGGER");

            if (EWOMS_GET_PARAM(PreTypeTag, bool, EclStrictParsing))
                parseContext.update( Ewoms::InputError::DELAYED_EXIT1);

            Ewoms::EFlowMain<PreTypeTag>::printPRTHeader(outputCout);

#ifdef HAVE_MPI
            Ewoms::ParallelEclipseState* parState;
#endif
            if (mpiRank == 0) {
                deck.reset( new Ewoms::Deck( parser.parseFile(deckFilename , parseContext, errorGuard)));
                Ewoms::MissingFeatures::checkKeywords(*deck, parseContext, errorGuard);
                if ( outputCout )
                    Ewoms::checkDeck(*deck, parser, parseContext, errorGuard);

#ifdef HAVE_MPI
                parState = new Ewoms::ParallelEclipseState(*deck);
                eclipseState.reset(parState);
#else
                eclipseState.reset(new Ewoms::EclipseState(*deck));
#endif
                /*
                  For the time being initializing wells and groups from the
                  restart file is not possible, but work is underways and it is
                  included here as a switch.
                */
                const bool init_from_restart_file = !EWOMS_GET_PARAM(PreTypeTag, bool, SchedRestart);
                const auto& init_config = eclipseState->getInitConfig();
                if (init_config.restartRequested() && init_from_restart_file) {
                    throw std::logic_error("Sorry - the ability to initialize wells and groups from the restart file is currently not ready");

                    int report_step = init_config.getRestartStep();
                    const auto& rst_filename = eclipseState->getIOConfig().getRestartFileName( init_config.getRestartRootName(), report_step, false );
                    Ewoms::EclIO::ERst rst_file(rst_filename);
                    const auto& rst_state = Ewoms::RestartIO::RstState::load(rst_file, report_step);
                    schedule.reset(new Ewoms::Schedule(*deck, *eclipseState, parseContext, errorGuard, &rst_state) );
                } else
                    schedule.reset(new Ewoms::Schedule(*deck, *eclipseState, parseContext, errorGuard));

                setupMessageLimiter(schedule->getMessageLimits(), "STDOUT_LOGGER");
                summaryConfig.reset( new Ewoms::SummaryConfig(*deck, *schedule, eclipseState->getTableManager(), parseContext, errorGuard));
#ifdef HAVE_MPI
                Ewoms::Mpi::packAndSend(*summaryConfig, Dune::MPIHelper::getCollectiveCommunication());
                Ewoms::Mpi::packAndSend(*schedule, Dune::MPIHelper::getCollectiveCommunication());
#endif
            }
#ifdef HAVE_MPI
            else {
                summaryConfig.reset(new Ewoms::SummaryConfig);
                schedule.reset(new Ewoms::Schedule);
                parState = new Ewoms::ParallelEclipseState;
                Ewoms::Mpi::receiveAndUnpack(*summaryConfig, mpiHelper.getCollectiveCommunication());
                Ewoms::Mpi::receiveAndUnpack(*schedule, mpiHelper.getCollectiveCommunication());
                eclipseState.reset(parState);
            }
            Ewoms::EclMpiSerializer ser(mpiHelper.getCollectiveCommunication());
            ser.broadcast(*parState);
#endif

            Ewoms::checkConsistentArrayDimensions(*eclipseState, *schedule, parseContext, errorGuard);

            if (errorGuard) {
                errorGuard.dump();
                errorGuard.clear();

                throw std::runtime_error("Unrecoverable errors were encountered while loading input.");
            }
        }
        const auto& phases = eclipseState->runspec().phases();
        bool outputFiles = (outputMode != FileOutputMode::OUTPUT_NONE);
        // run the actual simulator
        //
        // TODO: make sure that no illegal combinations like thermal and twophase are
        //       requested.

        if ( false ) {}
#ifndef FLOW_BLACKOIL_ONLY
        // Twophase cases
        else if( phases.size() == 2 ) {
            // oil-gas
            if (phases.active( Ewoms::Phase::GAS ))
            {
                Ewoms::eflowGasOilSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                return Ewoms::eflowGasOilMain(argc, argv, outputCout, outputFiles);
            }
            // oil-water
            else if ( phases.active( Ewoms::Phase::WATER ) )
            {
                Ewoms::eflowOilWaterSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                return Ewoms::eflowOilWaterMain(argc, argv, outputCout, outputFiles);
            }
            else {
                if (outputCout)
                    std::cerr << "No suitable configuration found, valid are Twophase (oilwater and oilgas), polymer, solvent, or blackoil" << std::endl;
                return EXIT_FAILURE;
            }
        }
        // Polymer case
        else if ( phases.active( Ewoms::Phase::POLYMER ) ) {

            if ( !phases.active( Ewoms::Phase::WATER) ) {
                if (outputCout)
                    std::cerr << "No valid configuration is found for polymer simulation, valid options include "
                              << "oilwater + polymer and blackoil + polymer" << std::endl;
                return EXIT_FAILURE;
            }

            // Need to track the polymer molecular weight
            // for the injectivity study
            if ( phases.active( Ewoms::Phase::POLYMW ) ) {
                // only oil water two phase for now
                assert( phases.size() == 4);
                return Ewoms::eflowOilWaterPolymerInjectivityMain(argc, argv, outputCout, outputFiles);
            }

            if ( phases.size() == 3 ) { // oil water polymer case
                Ewoms::eflowOilWaterPolymerSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                return Ewoms::eflowOilWaterPolymerMain(argc, argv, outputCout, outputFiles);
            } else {
                Ewoms::eflowPolymerSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
                return Ewoms::eflowPolymerMain(argc, argv, outputCout, outputFiles);
            }
        }
        // Foam case
        else if ( phases.active( Ewoms::Phase::FOAM ) ) {
            Ewoms::eflowFoamSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
            return Ewoms::eflowFoamMain(argc, argv, outputCout, outputFiles);
        }
        // Brine case
        else if ( phases.active( Ewoms::Phase::BRINE ) ) {
            Ewoms::eflowBrineSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
            return Ewoms::eflowBrineMain(argc, argv, outputCout, outputFiles);
        }
        // Solvent case
        else if ( phases.active( Ewoms::Phase::SOLVENT ) ) {
            Ewoms::eflowSolventSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
            return Ewoms::eflowSolventMain(argc, argv, outputCout, outputFiles);
        }
        // Energy case
        else if (eclipseState->getSimulationConfig().isThermal()) {
            Ewoms::eflowEnergySetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
            return Ewoms::eflowEnergyMain(argc, argv, outputCout, outputFiles);
        }
#endif // FLOW_BLACKOIL_ONLY
        // Blackoil case
        else if( phases.size() == 3 ) {
            Ewoms::eflowBlackoilSetDeck(externalSetupTimer.elapsed(), deck.get(), *eclipseState, *schedule, *summaryConfig);
            return Ewoms::eflowBlackoilMain(argc, argv, outputCout, outputFiles);
        }
        else
        {
            if (outputCout)
                std::cerr << "No suitable configuration found, valid are Twophase, polymer, foam, brine, solvent, energy, blackoil." << std::endl;
            return EXIT_FAILURE;
        }
    }
    catch (const std::invalid_argument& e)
    {
        if (outputCout) {
            std::cerr << "Failed to create valid EclipseState object." << std::endl;
            std::cerr << "Exception caught: " << e.what() << std::endl;
        }
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
