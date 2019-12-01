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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Ewoms::EclBaseVanguard
 */
#ifndef EWOMS_ECL_BASE_VANGUARD_HH
#define EWOMS_ECL_BASE_VANGUARD_HH

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

#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/opmlog/eclipseprtlog.hh>
#include <ewoms/eclio/opmlog/logutil.hh>

#include <boost/filesystem.hpp>

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
NEW_PROP_TAG(EclOutputInterval);
NEW_PROP_TAG(IgnoreKeywords);
NEW_PROP_TAG(EnableExperiments);
NEW_PROP_TAG(EdgeWeightsMethod);

SET_STRING_PROP(EclBaseVanguard, IgnoreKeywords, "");
SET_STRING_PROP(EclBaseVanguard, EclDeckFileName, "");
SET_INT_PROP(EclBaseVanguard, EclOutputInterval, -1); // use the deck-provided value
SET_BOOL_PROP(EclBaseVanguard, EnableEwomsRstFile, false);
SET_BOOL_PROP(EclBaseVanguard, EclStrictParsing, false);
SET_INT_PROP(EclBaseVanguard, EdgeWeightsMethod, 1);

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
    typedef BaseVanguard<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Vanguard) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

    enum { enableExperiments = GET_PROP_VALUE(TypeTag, EnableExperiments) };

public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

protected:
    static const int dimension = Grid::dimension;

    static void ensureOutputDirExists_(const std::string& outputDir)
    {
        if (!boost::filesystem::is_directory(outputDir)) {
            try {
                boost::filesystem::create_directories(outputDir);
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
        namespace fs = boost::filesystem;
        const std::string& deckFileNameParam = EWOMS_GET_PARAM(TypeTag, std::string, EclDeckFileName);
        fs::path deckFileNamePath(deckFileNameParam);

        // Strip extension "." or ".DATA"
        std::string extension = boost::to_upper_copy(deckFileNamePath.extension().string());
        std::string baseName;
        if (extension == ".DATA" || extension == ".")
            baseName = boost::to_upper_copy(deckFileNamePath.stem().string());
        else
            baseName = boost::to_upper_copy(deckFileNamePath.filename().string());

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
        EWOMS_REGISTER_PARAM(TypeTag, int, EdgeWeightsMethod,
                             "Choose edge-weighing strategy: 0=uniform, 1=trans, 2=log(trans).");
    }

    /*!
     * \brief Returns the canonical path to a deck file.
     *
     * The input can either be the canonical deck file name or the name of the case
     * (i.e., without the .DATA extension)
     */
    static boost::filesystem::path canonicalDeckPath(const std::string& caseName)
    {
        const auto fileExists = [](const boost::filesystem::path& f) -> bool
            {
                if (!boost::filesystem::exists(f))
                    return false;

                if (boost::filesystem::is_regular_file(f))
                    return true;

                return boost::filesystem::is_symlink(f) && boost::filesystem::is_regular_file(boost::filesystem::read_symlink(f));
            };

        auto simcase = boost::filesystem::path(caseName);
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

        std::unique_ptr<Ewoms::ParseContext> parseContext(new Ewoms::ParseContext(tmp));

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

        // Skip processing of filename if external deck already exists.
        if (!externalDeck_)
        {
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
        }

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

        if (!externalDeck_) {
            if (myRank == 0)
                std::cout << "Reading the deck file '" << fileName << "'" << std::endl;

            Ewoms::Parser parser;
            internalDeck_.reset(new Ewoms::Deck(parser.parseFile(fileName, *parseContext_, *errorGuard_)));
            deck_ = internalDeck_.get();

            if (enableExperiments && myRank == 0)
                Ewoms::checkDeck(*deck_, parser,  *parseContext_, *errorGuard_);
        }
        else {
            assert(externalDeck_);
            deck_ = externalDeck_;
        }

        if (!externalEclState_) {
            internalEclState_.reset(new Ewoms::EclipseState(*deck_, *parseContext_, *errorGuard_));
            eclState_ = internalEclState_.get();
        }
        else {
            assert(externalDeck_);
            assert(externalEclState_);

            deck_ = externalDeck_;
            eclState_ = externalEclState_;
        }

        if (!externalEclSchedule_) {
            // create the schedule object. Note that if eclState is supposed to represent
            // the internalized version of the deck, this constitutes a layering
            // violation.
            internalEclSchedule_.reset(new Ewoms::Schedule(*deck_, *eclState_, *parseContext_, *errorGuard_));
            eclSchedule_ = internalEclSchedule_.get();
        }
        else
            eclSchedule_ = externalEclSchedule_;
        this->summaryState_.reset( new Ewoms::SummaryState( std::chrono::system_clock::from_time_t(this->eclSchedule_->getStartTime() )));

        if (!externalEclSummaryConfig_) {
            // create the schedule object. Note that if eclState is supposed to represent
            // the internalized version of the deck, this constitutes a layering
            // violation.
            internalEclSummaryConfig_.reset(new Ewoms::SummaryConfig(*deck_,
                                                                   *eclSchedule_,
                                                                   eclState_->getTableManager(),
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

        // Possibly override IOConfig setting for how often RESTART files should get
        // written to disk (every N report step)
        int outputInterval = EWOMS_GET_PARAM(TypeTag, int, EclOutputInterval);
        if (outputInterval >= 0)
            eclState_->getRestartConfig().overrideRestartWriteInterval(outputInterval);
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
     * \brief Parameter deciding the edge-weight strategy of the load balancer.
     */
    Dune::EdgeWeightMethod edgeWeightsMethod() const
    { return edgeWeightsMethod_; }
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
     * \brief Return the names of the wells which do not penetrate any cells on the local
     *        process.
     *
     * This is a kludge around the fact that for distributed grids, not all wells are
     * seen by all proccesses.
     */
    std::unordered_set<std::string> defunctWellNames() const
    { return std::unordered_set<std::string>(); }

protected:
    void callImplementationInit()
    {
        asImp_().createGrids_();
        asImp_().filterConnections_();
        asImp_().updateOutputDir_();
        asImp_().finalizeInit_();

        if (enableExperiments) {
            Ewoms::RelpermDiagnostics relpermDiagnostics;
            relpermDiagnostics.diagnosis(*eclState_, *deck_, asImp_().grid());
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
        if (!boost::filesystem::is_directory(outputDir)) {
            try {
                boost::filesystem::create_directories(outputDir);
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
    std::unique_ptr<Ewoms::SummaryConfig> internalEclSummaryConfig_;
    std::unique_ptr<Ewoms::SummaryState> summaryState_;

    // these attributes point  either to the internal  or to the external version of the
    // parser objects.
    Ewoms::ParseContext* parseContext_;
    Ewoms::ErrorGuard* errorGuard_;
    Ewoms::Deck* deck_;
    Ewoms::EclipseState* eclState_;
    Ewoms::Schedule* eclSchedule_;
    Ewoms::SummaryConfig* eclSummaryConfig_;

    Dune::EdgeWeightMethod edgeWeightsMethod_;
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