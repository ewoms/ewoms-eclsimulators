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
 * \copydoc Ewoms::EclBaseVanguard
 */
#ifndef EWOMS_ECL_BASE_VANGUARD_HH
#define EWOMS_ECL_BASE_VANGUARD_HH

#include <ewoms/eclsimulators/utils/paralleleclipsestate.hh>
#include <ewoms/eclsimulators/utils/parallelserialization.hh>

#include <ewoms/numerics/io/basevanguard.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <ewoms/eclgrids/cpgrid.hh>
#include <ewoms/eclgrids/cpgrid/gridhelpers.hh>
#include <ewoms/eclsimulators/deprecated/props/satfunc/relpermdiagnostics.hh>

#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/parsecontext.hh>
#include <ewoms/eclio/parser/errorguard.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/eclipsestate/checkdeck.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/eclipsegrid.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/summaryconfig/summaryconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/summarystate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/state.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqstate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqconfig.hh>

#include <ewoms/common/string.hh>
#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/opmlog/eclipseprtlog.hh>
#include <ewoms/eclio/opmlog/logutil.hh>

#include <ewoms/common/filesystem.hh>

#include <dune/grid/common/mcmgmapper.hh>

#if HAVE_MPI
#include <mpi.h>
#endif // HAVE_MPI

#include <array>
#include <chrono>
#include <unordered_set>
#include <vector>

namespace Ewoms {
template <class TypeTag>
class EclBaseVanguard;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclBaseVanguard);

// declare the properties required by the for the ecl simulator vanguard
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(EquilGrid);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(EclDeckFileName);
NEW_PROP_TAG(OutputDir);
NEW_PROP_TAG(OutputMode);
NEW_PROP_TAG(EnableEwomsRstFile);
NEW_PROP_TAG(EclStrictParsing);
NEW_PROP_TAG(SchedRestart);
NEW_PROP_TAG(EclOutputInterval);
NEW_PROP_TAG(IgnoreKeywords);
NEW_PROP_TAG(EnableExperiments);
NEW_PROP_TAG(EdgeWeightsMethod);
NEW_PROP_TAG(OwnerCellsFirst);
NEW_PROP_TAG(ElementMapper);

NEW_PROP_TAG(SerialPartitioning);

SET_STRING_PROP(EclBaseVanguard, IgnoreKeywords, "");
SET_STRING_PROP(EclBaseVanguard, EclDeckFileName, "");
SET_INT_PROP(EclBaseVanguard, EclOutputInterval, -1); // use the deck-provided value
SET_BOOL_PROP(EclBaseVanguard, EnableEwomsRstFile, false);
SET_BOOL_PROP(EclBaseVanguard, EclStrictParsing, false);
SET_BOOL_PROP(EclBaseVanguard, SchedRestart, true);
SET_INT_PROP(EclBaseVanguard, EdgeWeightsMethod, 1);
SET_BOOL_PROP(EclBaseVanguard, OwnerCellsFirst, true);
SET_BOOL_PROP(EclBaseVanguard, SerialPartitioning, false);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 */
template <class TypeTag>
class EclBaseVanguard : public BaseVanguard<TypeTag>
{
    using ParentType = BaseVanguard<TypeTag>;
    using Implementation = GET_PROP_TYPE(TypeTag, Vanguard);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using ElementMapper = GET_PROP_TYPE(TypeTag, ElementMapper);

    enum { enableExperiments = GET_PROP_VALUE(TypeTag, EnableExperiments) };

public:
    using Grid = GET_PROP_TYPE(TypeTag, Grid);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);

protected:
    static const int dimension = Grid::dimension;
    using Element = typename GridView::template Codim<0>::Entity;

    static void ensureOutputDirExists_(const std::string& outputDir)
    {
        if (!Ewoms::filesystem::is_directory(outputDir)) {
            try {
                Ewoms::filesystem::create_directories(outputDir);
            }
            catch (const std::exception& e) {
                throw std::runtime_error("Creation of output directory '" + outputDir + "' failed: " + e.what() + "\n");
            }
            catch (...) {
                throw std::runtime_error("Creation of output directory '" + outputDir + "' failed for unknown reasons\n");
            }
        }
    }

public:
    /*!
     * \brief This method is called before the parameters are processed.
     */
    static void preDawn()
    {}

    /*!
     * \brief This method is called after the parameters are processed, but before the
     *        vanguard and simulator objects are instantiated.
     *
     * For ECL simulators, the logging infrastructure is set up here.
     */
    static void dawn()
    {
        int myRank = 0;
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

        // make sure that the directory specified via the --output-dir="..." parameter
        // exists
        const std::string& outputDir = EWOMS_GET_PARAM(TypeTag, std::string, OutputDir);
        if (!outputDir.empty())
            ensureOutputDirExists_(outputDir);

        // find the base name of the case, i.e., the file name without the extension
        namespace fs = Ewoms::filesystem;
        const std::string& deckFileNameParam = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
        fs::path deckFileNamePath(deckFileNameParam);

        // Strip extension "." or ".DATA"
        std::string extension = Ewoms::uppercase(deckFileNamePath.extension().string());
        std::string baseName;
        if (extension == ".DATA" || extension == ".")
            baseName = Ewoms::uppercase(deckFileNamePath.stem().string());
        else
            baseName = Ewoms::uppercase(deckFileNamePath.filename().string());

        // Stitch together the names of the debug and PRT log files
        std::string debugLogFileName = outputDir + "/" + baseName;
        std::string prtLogFileName = outputDir + "/" + baseName;
        if (myRank != 0) {
            // Added rank to log file for non-zero ranks.
            // This prevents message loss.
            debugLogFileName += "." + std::to_string(myRank);
            // If the following file appears then there is a bug.
            prtLogFileName += "." + std::to_string(myRank);
        }
        prtLogFileName += ".PRT";
        debugLogFileName += ".DBG";

        // parse the --output-mode="..." command line parameter
        enum class FileOutputMode {
            // No output to files.
            OutputNone = 0,
            // Output only to log files, no eclipse output.
            OutputLogOnly = 1,
            // Output to all files.
            OutputAll = 3
        };
        const std::string& outputModeParam = EWOMS_GET_PARAM(TypeTag, std::string, OutputMode);
        FileOutputMode outputMode = FileOutputMode::OutputAll;
        if (outputModeParam == "none" || outputModeParam == "false")
            outputMode = FileOutputMode::OutputNone;
        else if (outputModeParam == "log")
            outputMode = FileOutputMode::OutputLogOnly;
        else if (outputModeParam == "all" || outputModeParam == "true")
            outputMode = FileOutputMode::OutputAll;
        else
            std::cerr << "Value '" << outputModeParam << "' is not a recognized output mode. Enabling all output.\n";

        // add the actual logging backends
        if (outputMode > FileOutputMode::OutputNone) {
            std::shared_ptr<Ewoms::EclipsePRTLog> prtLog =
                std::make_shared<Ewoms::EclipsePRTLog>(prtLogFileName,
                                                     /*mask=*/Ewoms::Log::NoDebugMessageTypes,
                                                     /*append=*/false,
                                                     /*printSummary=*/true);
            prtLog->setMessageLimiter(std::make_shared<Ewoms::MessageLimiter>());
            prtLog->setMessageFormatter(std::make_shared<Ewoms::SimpleMessageFormatter>(/*useColorCoding=*/false));
            Ewoms::OpmLog::addBackend("ECLIPSEPRTLOG", prtLog);
        }

        if (outputMode >= FileOutputMode::OutputLogOnly) {
            std::shared_ptr<Ewoms::StreamLog> debugLog =
                std::make_shared<Ewoms::EclipsePRTLog>(debugLogFileName,
                                                     /*mask=*/Ewoms::Log::DefaultMessageTypes,
                                                     /*append=*/false,
                                                     /*printSummary=*/true);
            Ewoms::OpmLog::addBackend("DEBUGLOG", debugLog);
        }

        // on the master rank, print something to stdout
        if (myRank == 0) {
            auto stdoutLog = std::make_shared<Ewoms::StreamLog>(std::cout, /*mask=*/Ewoms::Log::StdoutMessageTypes);
            stdoutLog->setMessageFormatter(std::make_shared<Ewoms::SimpleMessageFormatter>(/*useColorCoding=*/true));
            Ewoms::OpmLog::addBackend("STDOUT_LOGGER", stdoutLog);
        }
    }

    /*!
     * \brief Register the common run-time parameters for all ECL simulator vanguards.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, std::string, EclDeckFileName,
                             "The name of the file which contains the ECL deck to be simulated");
        EWOMS_REGISTER_PARAM(TypeTag, int, EclOutputInterval,
                             "The number of report steps that ought to be skipped between two writes of ECL results");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableEwomsRstFile,
                             "Include eWoms-specific keywords in the ECL restart file to enable more accurate restarts of eWoms-based simulators from these files");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, IgnoreKeywords,
                             "List of Eclipse keywords which should be ignored. As a ':' separated string.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclStrictParsing,
                             "Use strict mode for parsing - all errors are collected before the applicaton exists.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, SchedRestart,
                             "When restarting: should we try to initialize wells and groups from historical SCHEDULE section.");
        EWOMS_REGISTER_PARAM(TypeTag, int, EdgeWeightsMethod,
                             "Choose edge-weighing strategy: 0=uniform, 1=trans, 2=log(trans).");
        EWOMS_REGISTER_PARAM(TypeTag, bool, OwnerCellsFirst,
                             "Order cells owned by rank before ghost/overlap cells.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, SerialPartitioning,
                             "Perform partitioning for parallel runs on a single process.");
    }

    /*!
     * \brief Returns the canonical path to a deck file.
     *
     * The input can either be the canonical deck file name or the name of the case
     * (i.e., without the .DATA extension)
     */
    static Ewoms::filesystem::path canonicalDeckPath(const std::string& caseName)
    {
        const auto fileExists = [](const Ewoms::filesystem::path& f) -> bool
            {
                if (!Ewoms::filesystem::exists(f))
                    return false;

                if (Ewoms::filesystem::is_regular_file(f))
                    return true;

                return Ewoms::filesystem::is_symlink(f) && Ewoms::filesystem::is_regular_file(Ewoms::filesystem::read_symlink(f));
            };

        auto simcase = Ewoms::filesystem::path(caseName);
        if (fileExists(simcase))
            return simcase;

        for (const auto& ext : { std::string("data"), std::string("DATA") }) {
            if (fileExists(simcase.replace_extension(ext)))
                return simcase;
        }

        throw std::invalid_argument("Cannot find input case '"+caseName+"'");
    }

    /*!
     * \brief Creates an Ewoms::parseContext object assuming that the parameters are ready.
     */
    static std::unique_ptr<Ewoms::ParseContext> createParseContext()
    {
        typedef std::pair<std::string, Ewoms::InputError::Action> ParseModePair;
        typedef std::vector<ParseModePair> ParseModePairs;
        ParseModePairs tmp;
        tmp.emplace_back(Ewoms::ParseContext::PARSE_RANDOM_SLASH, Ewoms::InputError::IGNORE);
        tmp.emplace_back(Ewoms::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Ewoms::InputError::WARN);
        tmp.emplace_back(Ewoms::ParseContext::SUMMARY_UNKNOWN_WELL, Ewoms::InputError::WARN);
        tmp.emplace_back(Ewoms::ParseContext::SUMMARY_UNKNOWN_GROUP, Ewoms::InputError::WARN);
        tmp.emplace_back(Ewoms::ParseContext::PARSE_EXTRA_RECORDS, Ewoms::InputError::WARN);
        tmp.emplace_back(Ewoms::ParseContext::PARSE_WGNAME_SPACE, Ewoms::InputError::WARN);

        auto parseContext = std::make_unique<Ewoms::ParseContext>(tmp);

        const std::string ignoredKeywords = EWOMS_GET_PARAM(TypeTag, std::string, IgnoreKeywords);
        if (ignoredKeywords.size() > 0) {
            size_t pos;
            size_t offset = 0;
            while (true) {
                pos = ignoredKeywords.find(':', offset);
                if (pos == std::string::npos) {
                    parseContext->ignoreKeyword(ignoredKeywords.substr(offset));
                    break;
                }
                parseContext->ignoreKeyword(ignoredKeywords.substr(offset, pos - offset));
                offset = pos + 1;
            }
        }

        if (EWOMS_GET_PARAM(TypeTag, bool, EclStrictParsing))
            parseContext->update(Ewoms::InputError::DELAYED_EXIT1);

        return parseContext;
    }

    /*!
     * \brief Set the wall time which was spend externally to set up the external data structures
     *
     * i.e., the objects specified via the other setExternal*() methods.
     */
    static void setExternalSetupTime(Scalar t)
    { externalSetupTime_ = t; }

    /*!
     * \brief Returns the wall time required to set up the simulator before it was born.
     */
    static Scalar externalSetupTime()
    { return externalSetupTime_; }

    /*!
     * \brief Set the Ewoms::ParseContext object which ought to be used for parsing the deck and creating the Ewoms::EclipseState object.
     */
    static void setExternalParseContext(Ewoms::ParseContext* parseContext)
    { externalParseContext_ = parseContext; }

    /*!
     * \brief Set the Ewoms::ErrorGuard object which ought to be used for parsing the deck and creating the Ewoms::EclipseState object.
     */
    static void setExternalErrorGuard(Ewoms::ErrorGuard* errorGuard)
    { externalErrorGuard_ = errorGuard; }

    /*!
     * \brief Set the Ewoms::Deck object which ought to be used when the simulator vanguard
     *        is instantiated.
     *
     * This is basically an optimization: In cases where the ECL input deck must be
     * examined to decide which simulator ought to be used, this avoids having to parse
     * the input twice. When this method is used, the caller is responsible for lifetime
     * management of these two objects, i.e., they are not allowed to be deleted as long
     * as the simulator vanguard object is alive.
     */
    static void setExternalDeck(Ewoms::Deck* deck)
    { externalDeck_ = deck; }

    /*!
     * \brief Set the Ewoms::EclipseState object which ought to be used when the simulator
     *        vanguard is instantiated.
     */
    static void setExternalEclState(Ewoms::EclipseState* eclState)
    { externalEclState_ = eclState; }

    /*!
     * \brief Create the grid for problem data files which use the ECL file format.
     *
     * This is the file format used by the commercial ECLIPSE simulator. Usually it uses
     * a cornerpoint description of the grid.
     */
    EclBaseVanguard(Simulator& simulator)
        : ParentType(simulator)
    {
        int myRank = 0;
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

        std::string fileName = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
        edgeWeightsMethod_   = Dune::EdgeWeightMethod(EWOMS_GET_PARAM(TypeTag, int, EdgeWeightsMethod));
        ownersFirst_ = EWOMS_GET_PARAM(TypeTag, bool, OwnerCellsFirst);
        serialPartitioning_ = EWOMS_GET_PARAM(TypeTag, bool, SerialPartitioning);

        // Make proper case name.
        if (fileName == "")
            throw std::runtime_error("No input deck file has been specified as a command line argument,"
                                     " or via '--ecl-deck-file-name=CASE.DATA'");

        fileName = canonicalDeckPath(fileName).string();

        // compute the base name of the input file name
        const char directorySeparator = '/';
        long int i;
        for (i = fileName.size(); i >= 0; -- i)
            if (fileName[i] == directorySeparator)
                break;
        std::string baseName = fileName.substr(i + 1, fileName.size());

        // remove the extension from the input file
        for (i = baseName.size(); i >= 0; -- i)
            if (baseName[i] == '.')
                break;
        std::string rawCaseName;
        if (i < 0)
            rawCaseName = baseName;
        else
            rawCaseName = baseName.substr(0, i);

        // transform the result to ALL_UPPERCASE
        caseName_ = rawCaseName;
        std::transform(caseName_.begin(), caseName_.end(), caseName_.begin(), ::toupper);

        // create the parser objects for the deck or use their externally specified
        // versions (if desired)
        if (!externalParseContext_) {
            internalParseContext_ = createParseContext();
            parseContext_ = internalParseContext_.get();
        }
        else
            parseContext_ = externalParseContext_;

        if (!externalParseContext_) {
            internalErrorGuard_.reset(new Ewoms::ErrorGuard);
            errorGuard_ = internalErrorGuard_.get();
        }
        else
            errorGuard_ = externalErrorGuard_;

        if (myRank == 0)
            setupEclObjectsOnMaster_(fileName);
        else {
            internalEclSummaryConfig_.reset(new Ewoms::SummaryConfig);
            internalEclSchedule_.reset(new Ewoms::Schedule);
            internalEclState_.reset(new Ewoms::ParallelEclipseState);

            eclSummaryConfig_ = internalEclSummaryConfig_.get();
            eclSchedule_ = internalEclSchedule_.get();
            eclState_ = internalEclState_.get();
        }
#if HAVE_MPI
        Ewoms::eclStateBroadcast(*eclState_, *eclSchedule_, *eclSummaryConfig_);
#endif

        summaryState_.reset(new Ewoms::SummaryState(std::chrono::system_clock::from_time_t(this->eclSchedule_->getStartTime())));
        actionState_.reset(new Ewoms::Action::State());

        // Possibly override IOConfig setting for how often RESTART files should get
        // written to disk (every N report steps)
        int outputInterval = EWOMS_GET_PARAM(TypeTag, int, EclOutputInterval);
        if (outputInterval >= 0)
            schedule().restart().overrideRestartWriteInterval(outputInterval);
    }

    void setupEclObjectsOnMaster_(const std::string& fileName)
    {
        if (!externalDeck_) {
            std::cout << "Reading the deck file '" << fileName << "'" << std::endl;

            Ewoms::Parser parser;
            internalDeck_.reset(new Ewoms::Deck(parser.parseFile(fileName, *parseContext_, *errorGuard_)));
            deck_ = internalDeck_.get();

            if (enableExperiments)
                Ewoms::checkDeck(*deck_, parser,  *parseContext_, *errorGuard_);
        }
        else
            deck_ = externalDeck_;

        if (!externalEclState_) {
#if HAVE_MPI
            internalEclState_.reset(new Ewoms::ParallelEclipseState(*deck_));
#else
            internalEclState_.reset(new Ewoms::EclipseState(*deck_));
#endif

            eclState_ = internalEclState_.get();
        }
        else {
            assert(externalEclState_);

            deck_ = externalDeck_;
            eclState_ = externalEclState_;
        }

        if (!externalEclSchedule_) {
            // create the schedule object. Note that if eclState is supposed to represent
            // the internalized version of the deck, this constitutes a layering
            // violation.

            const bool initFromRestartFile = !EWOMS_GET_PARAM(TypeTag, bool, SchedRestart);
            const auto& initConfig = eclState_->getInitConfig();
            if (initConfig.restartRequested() && initFromRestartFile) {
                int reportStepNum = initConfig.getRestartStep();
                const auto& rstFileName = eclState_->getIOConfig().getRestartFileName(initConfig.getRestartRootName(), reportStepNum, false);
                Ewoms::EclIO::ERst rstFile(rstFileName);
                const auto& rstState = Ewoms::RestartIO::RstState::load(rstFile, reportStepNum);
                internalEclSchedule_.reset(new Ewoms::Schedule(*deck_, *eclState_, *parseContext_, *errorGuard_, &rstState));
            }
            else
                internalEclSchedule_.reset(new Ewoms::Schedule(*deck_, *eclState_, *parseContext_, *errorGuard_));

            eclSchedule_ = internalEclSchedule_.get();
        }
        else
            eclSchedule_ = externalEclSchedule_;
        this->summaryState_.reset( new Ewoms::SummaryState( std::chrono::system_clock::from_time_t(this->eclSchedule_->getStartTime() )));
        this->actionState_.reset( new Ewoms::Action::State() );
        this->aquiferConfig_.reset( new Ewoms::AquiferConfig() );

        if (!externalEclSummaryConfig_) {
            // create the schedule object. Note that if eclState is supposed to represent
            // the internalized version of the deck, this constitutes a layering
            // violation.
            internalEclSummaryConfig_.reset(new Ewoms::SummaryConfig(*deck_,
                                                                   *eclSchedule_,
                                                                   eclState_->getTableManager(),
                                                                   *aquiferConfig_,
                                                                   *parseContext_,
                                                                   *errorGuard_));

            eclSummaryConfig_ = internalEclSummaryConfig_.get();
        }
        else
            eclSummaryConfig_ = externalEclSummaryConfig_;

        if (*errorGuard_) {
            errorGuard_->dump();
            errorGuard_->clear();

            throw std::runtime_error("Unrecoverable errors were encountered while loading input.");
        }

        udqState_.reset(new Ewoms::UDQState( this->eclSchedule_->getUDQConfig(0).params().undefinedValue()));

        // Possibly override IOConfig setting for how often RESTART files should get
        // written to disk (every N report step)
        int outputInterval = EWOMS_GET_PARAM(TypeTag, int, EclOutputInterval);
        if (outputInterval >= 0)
            schedule().restart().overrideRestartWriteInterval(outputInterval);

        // Initialize parallelWells with all local wells
        const auto& schedule_wells = schedule().getWellsatEnd();
        parallelWells_.reserve(schedule_wells.size());

        for (const auto& well: schedule_wells)
        {
            parallelWells_.emplace_back(well.name(), true);
        }
        std::sort(parallelWells_.begin(), parallelWells_.end());
    }

    /*!
     * \brief Return a reference to the parsed ECL deck.
     */
    const Ewoms::Deck& deck() const
    { return *deck_; }

    Ewoms::Deck& deck()
    { return *deck_; }

    /*!
     * \brief Return a reference to the internalized ECL deck.
     */
    const Ewoms::EclipseState& eclState() const
    { return *eclState_; }

    Ewoms::EclipseState& eclState()
    { return *eclState_; }

    /*!
     * \brief Return a reference to the object that managages the ECL schedule.
     */
    const Ewoms::Schedule& schedule() const
    { return *eclSchedule_; }

    Ewoms::Schedule& schedule()
    { return *eclSchedule_; }

    /*!
     * \brief Set the schedule object.
     *
     * The lifetime of this object is not managed by the vanguard, i.e., the object must
     * stay valid until after the vanguard gets destroyed.
     */
    static void setExternalSchedule(Ewoms::Schedule* schedule)
    { externalEclSchedule_ = schedule; }

    /*!
     * \brief Return a reference to the object that determines which quantities ought to
     *        be put into the ECL summary output.
     */
    const Ewoms::SummaryConfig& summaryConfig() const
    { return *eclSummaryConfig_; }

    /*!
     * \brief Set the summary configuration object.
     *
     * The lifetime of this object is not managed by the vanguard, i.e., the object must
     * stay valid until after the vanguard gets destroyed.
     */
    static void setExternalSummaryConfig(Ewoms::SummaryConfig* summaryConfig)
    { externalEclSummaryConfig_ = summaryConfig; }

    /*!
    * \brief Returns the summary state
    *
    * The summary state is a small container object for
    * computed, ready to use summary values. The values will typically be used by
    * the UDQ, WTEST and ACTIONX calculations.
    */
    Ewoms::SummaryState& summaryState()
    { return *summaryState_; }

    const Ewoms::SummaryState& summaryState() const
    { return *summaryState_; }

    /*!
     * \brief Returns the action state
     *
     * The ActionState keeps track of how many times the various actions have run.
     */
    Ewoms::Action::State& actionState()
    { return *actionState_; }

    const Ewoms::Action::State& actionState() const
    { return *actionState_; }

    /*!
     * \brief Returns the udq state
     *
     * The UDQState keeps track of the result of user defined quantity (UDQ) evaluations that can
     * be specified in the ECL format.
     */
    Ewoms::UDQState& udqState()
    { return *udqState_; }

    const Ewoms::UDQState& udqState() const
    { return *udqState_; }

    /*!
     * \brief Parameter deciding the edge-weight strategy of the load balancer.
     */
    Dune::EdgeWeightMethod edgeWeightsMethod() const
    { return edgeWeightsMethod_; }

    /*!
     * \brief Parameter that decide if cells owned by rank are ordered before ghost cells.
     */
    bool ownersFirst() const
    { return ownersFirst_; }

    /*!
     * \brief Parameter that decides if partitioning for parallel runs
     *        should be performed on a single process only.
     */
    bool serialPartitioning() const
    { return serialPartitioning_; }

    /*!
     * \brief Returns the name of the case.
     *
     * i.e., the all-uppercase version of the file name from which the
     * deck is loaded with the ".DATA" suffix removed.
     */
    const std::string& caseName() const
    { return caseName_; }

    /*!
     * \brief Returns the number of logically Cartesian cells in each direction
     */
    const std::array<int, dimension>& cartesianDimensions() const
    { return asImp_().cartesianIndexMapper().cartesianDimensions(); }

    /*!
     * \brief Returns the overall number of cells of the logically Cartesian grid
     */
    int cartesianSize() const
    { return asImp_().cartesianIndexMapper().cartesianSize(); }

    /*!
     * \brief Returns the overall number of cells of the logically EquilCartesian grid
     */
    int equilCartesianSize() const
    { return asImp_().equilCartesianIndexMapper().cartesianSize(); }

    /*!
     * \brief Returns the Cartesian cell id for identifaction with ECL data
     */
    unsigned cartesianIndex(unsigned compressedCellIdx) const
    { return asImp_().cartesianIndexMapper().cartesianIndex(compressedCellIdx); }

    /*!
     * \brief Return the index of the cells in the logical Cartesian grid
     */
    unsigned cartesianIndex(const std::array<int,dimension>& coords) const
    {
        unsigned cartIndex = coords[0];
        int factor = cartesianDimensions()[0];
        for (unsigned i = 1; i < dimension; ++i) {
            cartIndex += coords[i]*factor;
            factor *= cartesianDimensions()[i];
        }

        return cartIndex;
    }

    /*!
     * \brief Return compressed index from cartesian index
     *
     */
    int compressedIndex(int cartesianCellIdx) const
    {
        int index = cartesianToCompressed_[cartesianCellIdx];
        return index;
    }

    /*!
     * \brief Extract Cartesian index triplet (i,j,k) of an active cell.
     *
     * \param [in] cellIdx Active cell index.
     * \param [out] ijk Cartesian index triplet
     */
    void cartesianCoordinate(unsigned cellIdx, std::array<int,3>& ijk) const
    { return asImp_().cartesianIndexMapper().cartesianCoordinate(cellIdx, ijk); }

    /*!
     * \brief Returns the Cartesian cell id given an element index for the grid used for equilibration
     */
    unsigned equilCartesianIndex(unsigned compressedEquilCellIdx) const
    { return asImp_().equilCartesianIndexMapper().cartesianIndex(compressedEquilCellIdx); }

    /*!
     * \brief Extract Cartesian index triplet (i,j,k) of an active cell of the grid used for EQUIL.
     *
     * \param [in] cellIdx Active cell index.
     * \param [out] ijk Cartesian index triplet
     */
    void equilCartesianCoordinate(unsigned cellIdx, std::array<int,3>& ijk) const
    { return asImp_().equilCartesianIndexMapper().cartesianCoordinate(cellIdx, ijk); }

    /*!
     * \brief Returns vector with name and whether the has local perforated cells
     *        for all wells.
     *
     * Will only have usable values for CpGrid.
     */
    const std::vector<std::pair<std::string,bool>>& parallelWells() const
    { return parallelWells_; }

    /*!
     * \brief Get the cell centroids for a distributed grid.
     *
     * Currently this only non-empty for a loadbalanced CpGrid.
     */
    const std::vector<double>& cellCentroids() const
    {
        return centroids_;
    }

    /*!
     * \brief Returns the depth of an degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depth of an element and is
     * thus slightly different from the depth of an element's centroid.
     */
    Scalar cellCenterDepth(unsigned globalSpaceIdx) const
    {
        return cellCenterDepth_[globalSpaceIdx];
    }

    /*!
     * \brief Get the number of cells in the global leaf grid view.
     * \warn This is a collective operation that needs to be called
     * on all ranks.
     */
    std::size_t globalNumCells() const
    {
        const auto& grid = asImp_().grid();
        if (grid.comm().size() == 1)
        {
            return grid.leafGridView().size(0);
        }
        const auto& gridView = grid.leafGridView();
        constexpr int codim = 0;
        constexpr auto Part = Dune::Interior_Partition;
        auto local_cells = std::distance(gridView.template begin<codim, Part>(),
                                         gridView.template end<codim, Part>());
        return grid.comm().sum(local_cells);
    }

protected:
    void callImplementationInit()
    {
        asImp_().createGrids_();
        asImp_().filterConnections_();
        asImp_().updateOutputDir_();
        asImp_().finalizeInit_();

        int myRank = 0;
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
#endif

        if (enableExperiments && myRank == 0) {
            Ewoms::RelpermDiagnostics relpermDiagnostics;
            relpermDiagnostics.diagnosis(*eclState_, asImp_().grid());
        }
    }
    void updateCartesianToCompressedMapping_()
    {
        size_t num_cells = asImp_().grid().leafGridView().size(0);
        cartesianToCompressed_.resize(cartesianSize(), -1);
        for (unsigned i = 0; i < num_cells; ++i) {
            unsigned cartesianCellIdx = cartesianIndex(i);
            cartesianToCompressed_[cartesianCellIdx] = i;
        }
    }

    void updateCellDepths_()
    {
        int numCells = this->gridView().size(/*codim=*/0);
        cellCenterDepth_.resize(numCells);

#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());
#else
        ElementMapper elemMapper(this->gridView());
#endif
        auto elemIt = this->gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = this->gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& element = *elemIt;
            const unsigned int elemIdx = elemMapper.index(element);
            cellCenterDepth_[elemIdx] = cellCenterDepth(element);
        }
    }

private:
    void updateOutputDir_()
    {
        // update the location for output
        std::string outputDir = EWOMS_GET_PARAM(TypeTag, std::string, OutputDir);
        auto& ioConfig = eclState_->getIOConfig();
        if (outputDir == "")
            // If no output directory parameter is specified, use the output directory
            // which Ewoms::IOConfig thinks that should be used. Normally this is the
            // directory in which the input files are located.
            outputDir = ioConfig.getOutputDir();

        // ensure that the output directory exists and that it is a directory
        if (!Ewoms::filesystem::is_directory(outputDir)) {
            try {
                Ewoms::filesystem::create_directories(outputDir);
            }
            catch (...) {
                 throw std::runtime_error("Creation of output directory '"+outputDir+"' failed\n");
            }
        }

        // specify the directory output. This is not a very nice mechanism because
        // the eclState is supposed to be immutable here, IMO.
        ioConfig.setOutputDir(outputDir);

        ioConfig.setEclCompatibleRST(!EWOMS_GET_PARAM(TypeTag, bool, EnableEwomsRstFile));
    }

    Scalar cellCenterDepth(const Element& element) const
    {
        typedef typename Element::Geometry Geometry;
        static constexpr int zCoord = Element::dimension - 1;
        Scalar zz = 0.0;

        const Geometry geometry = element.geometry();
        const int corners = geometry.corners();
        for (int i=0; i < corners; ++i)
            zz += geometry.corner(i)[zCoord];

        return zz/Scalar(corners);
    }

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    const Implementation& asImp_() const
    { return *static_cast<const Implementation*>(this); }

    std::string caseName_;

    static Scalar externalSetupTime_;

    static Ewoms::ParseContext* externalParseContext_;
    static Ewoms::ErrorGuard* externalErrorGuard_;
    static Ewoms::Deck* externalDeck_;
    static Ewoms::EclipseState* externalEclState_;
    static Ewoms::Schedule* externalEclSchedule_;
    static Ewoms::SummaryConfig* externalEclSummaryConfig_;

    std::unique_ptr<Ewoms::ParseContext> internalParseContext_;
    std::unique_ptr<Ewoms::ErrorGuard> internalErrorGuard_;
    std::unique_ptr<Ewoms::Deck> internalDeck_;
    std::unique_ptr<Ewoms::EclipseState> internalEclState_;
    std::unique_ptr<Ewoms::Schedule> internalEclSchedule_;
    std::unique_ptr<Ewoms::AquiferConfig> aquiferConfig_;
    std::unique_ptr<Ewoms::SummaryConfig> internalEclSummaryConfig_;
    std::unique_ptr<Ewoms::SummaryState> summaryState_;
    std::unique_ptr<Ewoms::Action::State> actionState_;
    std::unique_ptr<Ewoms::UDQState> udqState_;

    // these attributes point  either to the internal  or to the external version of the
    // parser objects.
    Ewoms::ParseContext* parseContext_;
    Ewoms::ErrorGuard* errorGuard_;
    Ewoms::Deck* deck_;
    Ewoms::EclipseState* eclState_;
    Ewoms::Schedule* eclSchedule_;
    Ewoms::SummaryConfig* eclSummaryConfig_;

    Dune::EdgeWeightMethod edgeWeightsMethod_;
    bool ownersFirst_;
    bool serialPartitioning_;

protected:
    /*! \brief The cell centroids after loadbalance was called.
     * Empty otherwise. Used by EclTransmissibilty.
     */
    std::vector<double> centroids_;

    /*! \brief Mapping between cartesian and compressed cells.
     *  It is initialized the first time it is called
     */
    std::vector<int> cartesianToCompressed_;

    /*! \brief Cell center depths computed
     *  from averaging cell corner depths
     */
    std::vector<Scalar> cellCenterDepth_;

    /*! \brief information about wells in parallel
     *
     * For each well in the model there is an entry with its name
     * and a boolean indicating whether it perforates local cells.
     */
    std::vector<std::pair<std::string,bool>> parallelWells_;

};

template <class TypeTag>
typename EclBaseVanguard<TypeTag>::Scalar EclBaseVanguard<TypeTag>::externalSetupTime_ = 0.0;

template <class TypeTag>
Ewoms::ParseContext* EclBaseVanguard<TypeTag>::externalParseContext_ = nullptr;

template <class TypeTag>
Ewoms::ErrorGuard* EclBaseVanguard<TypeTag>::externalErrorGuard_ = nullptr;

template <class TypeTag>
Ewoms::Deck* EclBaseVanguard<TypeTag>::externalDeck_ = nullptr;

template <class TypeTag>
Ewoms::EclipseState* EclBaseVanguard<TypeTag>::externalEclState_;

template <class TypeTag>
Ewoms::Schedule* EclBaseVanguard<TypeTag>::externalEclSchedule_ = nullptr;

template <class TypeTag>
Ewoms::SummaryConfig* EclBaseVanguard<TypeTag>::externalEclSummaryConfig_ = nullptr;

} // namespace Ewoms

#endif
