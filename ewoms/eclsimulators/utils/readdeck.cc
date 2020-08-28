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

#if HAVE_MPI
#include "mpi.h"
#endif
#include "readdeck.hh"

#include <ewoms/common/string.hh>

#include <ewoms/eclio/io/ecliodata.hh>

#include <ewoms/eclio/output/restartio.hh>
#include <ewoms/eclio/io/rst/state.hh>

#include <ewoms/eclio/parser/eclipsestate/checkdeck.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/arraydimchecker.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/errorguard.hh>

#include <ewoms/eclsimulators/eflow/missingfeatures.hh>
#include <ewoms/eclsimulators/utils/paralleleclipsestate.hh>
#include <ewoms/eclsimulators/utils/parallelserialization.hh>

namespace Ewoms
{

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
FileOutputMode setupLogging(int mpi_rank_, const std::string& deck_filename, const std::string& cmdline_output_dir, const std::string& cmdline_output, bool output_cout_, const std::string& stdout_log_id) {

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
    std::string extension = uppercase(fpath.extension().string());
    if (extension == ".DATA" || extension == ".") {
        baseName = uppercase(fpath.stem().string());
    } else {
        baseName = uppercase(fpath.filename().string());
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

void readDeck(int rank, std::string& deckFilename, std::unique_ptr<Ewoms::Deck>& deck, std::unique_ptr<Ewoms::EclipseState>& eclipseState,
              std::unique_ptr<Ewoms::Schedule>& schedule, std::unique_ptr<Ewoms::SummaryConfig>& summaryConfig,
              std::unique_ptr<ErrorGuard> errorGuard, std::unique_ptr<ParseContext> parseContext,
              bool initFromRestart, bool checkDeck)
{
    if (!errorGuard)
    {
        errorGuard = std::make_unique<ErrorGuard>();
    }

#if HAVE_MPI
    int parseSuccess = 0;
#endif
    std::string failureMessage;

    if (rank==0) {
        try
        {
            if ( (!deck || !schedule || !summaryConfig ) && !parseContext)
            {
                EWOMS_THROW(std::logic_error, "We need a parse context if deck, schedule, or summaryConfig are not initialized");
            }

            if (!deck)
            {
                Ewoms::Parser parser;
                deck = std::make_unique<Ewoms::Deck>( parser.parseFile(deckFilename , *parseContext, *errorGuard));
                Ewoms::MissingFeatures::checkKeywords(*deck, *parseContext, *errorGuard);
                if ( checkDeck )
                    Ewoms::checkDeck(*deck, parser, *parseContext, *errorGuard);
            }

            if (!eclipseState) {
#if HAVE_MPI
                eclipseState = std::make_unique<Ewoms::ParallelEclipseState>(*deck);
#else
                eclipseState = std::make_unique<Ewoms::EclipseState>(*deck);
#endif
            }
            /*
              For the time being initializing wells and groups from the
              restart file is not possible, but work is underways and it is
              included here as a switch.
            */
            const auto& init_config = eclipseState->getInitConfig();
            if (init_config.restartRequested() && initFromRestart) {
                int report_step = init_config.getRestartStep();
                const auto& rst_filename = eclipseState->getIOConfig().getRestartFileName( init_config.getRestartRootName(), report_step, false );
                Ewoms::EclIO::ERst rst_file(rst_filename);
                const auto& rst_state = Ewoms::RestartIO::RstState::load(rst_file, report_step);
                if (!schedule)
                    schedule = std::make_unique<Ewoms::Schedule>(*deck, *eclipseState, *parseContext, *errorGuard, &rst_state);
            }
            else {
                if (!schedule)
                    schedule = std::make_unique<Ewoms::Schedule>(*deck, *eclipseState, *parseContext, *errorGuard);
            }
            if (Ewoms::OpmLog::hasBackend("STDOUT_LOGGER")) // loggers might not be set up!
            {
                setupMessageLimiter(schedule->getMessageLimits(), "STDOUT_LOGGER");
            }
            if (!summaryConfig)
                summaryConfig = std::make_unique<Ewoms::SummaryConfig>(*deck, *schedule, eclipseState->getTableManager(), *parseContext, *errorGuard);
#if HAVE_MPI
            parseSuccess = 1;
#endif
        }
        catch(const std::exception& e)
        {
            failureMessage = e.what();
        }
    }
#if HAVE_MPI
    else {
        if (!summaryConfig)
            summaryConfig = std::make_unique<Ewoms::SummaryConfig>();
        if (!schedule)
            schedule = std::make_unique<Ewoms::Schedule>();
        if (!eclipseState)
            eclipseState = std::make_unique<Ewoms::ParallelEclipseState>();
    }

    auto comm = Dune::MPIHelper::getCollectiveCommunication();
    parseSuccess = comm.max(parseSuccess);
    if (!parseSuccess)
    {
        if (*errorGuard) {
            errorGuard->dump();
            errorGuard->clear();
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    Ewoms::eclStateBroadcast(*eclipseState, *schedule, *summaryConfig);
#endif

    Ewoms::checkConsistentArrayDimensions(*eclipseState, *schedule, *parseContext, *errorGuard);

    if (*errorGuard) {
        errorGuard->dump();
        errorGuard->clear();

        throw std::runtime_error("Unrecoverable errors were encountered while loading input.");
    }
}
} // end namespace Ewoms
