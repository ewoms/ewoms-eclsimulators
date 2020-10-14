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

#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/checkdeck.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/arraydimchecker.hh>

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <ewoms/eclsimulators/eflow/eflowmain.hh>
#include <ewoms/eclsimulators/utils/readdeck.hh>

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
  void eflowSetDeck(Deck& deck, EclipseState& eclState, Schedule& schedule, SummaryConfig& summaryConfig)
  {
    using Vanguard = GET_PROP_TYPE(TypeTag, Vanguard);
    Vanguard::setExternalDeck(&deck);
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
    Ewoms::EFlowMain<TypeTag> mainfunc(argc, argv, outputCout, outputFiles);
    return mainfunc.execute();
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
             std::unique_ptr<Ewoms::Deck> deck,
             std::unique_ptr<Ewoms::EclipseState> eclipseState,
             std::unique_ptr<Ewoms::Schedule> schedule,
             std::unique_ptr<Ewoms::SummaryConfig> summaryConfig)
            : argc_(argc)
            , argv_(argv)
            , deck_(std::move(deck))
            , eclipseState_(std::move(eclipseState))
            , schedule_(std::move(schedule))
            , summaryConfig_(std::move(summaryConfig))
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
                    *deck_,
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
                    Ewoms::eflowGasOilSetDeck(setupTime_, *deck_, *eclipseState_,
                                               *schedule_, *summaryConfig_);
                    return Ewoms::eflowGasOilMain(argc_, argv_, outputCout_, outputFiles_);
                }
                // oil-water
                else if ( phases.active( Ewoms::Phase::WATER ) ) {
                    Ewoms::eflowOilWaterSetDeck(setupTime_, *deck_, *eclipseState_, *schedule_, *summaryConfig_);
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
                    Ewoms::eflowOilWaterPolymerSetDeck(setupTime_, *deck_,
                                                        *eclipseState_,
                                                        *schedule_,
                                                        *summaryConfig_);
                    return Ewoms::eflowOilWaterPolymerMain(argc_, argv_, outputCout_, outputFiles_);
                } else {
                    Ewoms::eflowPolymerSetDeck(setupTime_, *deck_,
                                                *eclipseState_,
                                                *schedule_,
                                                *summaryConfig_);
                    return Ewoms::eflowPolymerMain(argc_, argv_, outputCout_, outputFiles_);
                }
            }
            // Foam case
            else if ( phases.active( Ewoms::Phase::FOAM ) ) {
                Ewoms::eflowFoamSetDeck(setupTime_, *deck_,
                                         *eclipseState_,
                                         *schedule_,
                                         *summaryConfig_);
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
                    Ewoms::eflowOilWaterBrineSetDeck(setupTime_, *deck_,
                                                      *eclipseState_,
                                                      *schedule_,
                                                      *summaryConfig_);
                    return Ewoms::eflowOilWaterBrineMain(argc_, argv_, outputCout_, outputFiles_);
                } else {
                    Ewoms::eflowBrineSetDeck(setupTime_, *deck_,
                                              *eclipseState_,
                                              *schedule_,
                                              *summaryConfig_);
                    return Ewoms::eflowBrineMain(argc_, argv_, outputCout_, outputFiles_);
                }
            }
            // Solvent case
            else if ( phases.active( Ewoms::Phase::SOLVENT ) ) {
                Ewoms::eflowSolventSetDeck(setupTime_, *deck_,
                                            *eclipseState_,
                                            *schedule_,
                                            *summaryConfig_);
                return Ewoms::eflowSolventMain(argc_, argv_, outputCout_, outputFiles_);
            }
            // Energy case
            else if (eclipseState_->getSimulationConfig().isThermal()) {
                Ewoms::eflowEnergySetDeck(setupTime_, *deck_,
                                           *eclipseState_,
                                           *schedule_,
                                           *summaryConfig_);
                return Ewoms::eflowEnergyMain(argc_, argv_, outputCout_, outputFiles_);
            }
#endif // FLOW_BLACKOIL_ONLY
            // Blackoil case
            else if( phases.size() == 3 ) {
                Ewoms::eflowBlackoilSetDeck(setupTime_, *deck_,
                                             *eclipseState_,
                                             *schedule_,
                                             *summaryConfig_);
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
            Ewoms::eflowSetDeck<TypeTag>(*deck_,
                                         *eclipseState_,
                                         *schedule_,
                                         *summaryConfig_);
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
            using PreTypeTag = TTAG(EFlowEarlyBird);
            using PreProblem = GET_PROP_TYPE(PreTypeTag, Problem);

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

            using PreVanguard = GET_PROP_TYPE(PreTypeTag, Vanguard);
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
                const bool init_from_restart_file = !EWOMS_GET_PARAM(PreTypeTag, bool, SchedRestart);
                if (outputDir.empty())
                    outputDir = EWOMS_GET_PARAM(PreTypeTag, std::string, OutputDir);
                outputMode = setupLogging(mpiRank,
                                          deckFilename,
                                          outputDir,
                                          EWOMS_GET_PARAM(PreTypeTag, std::string, OutputMode),
                                          outputCout_, "STDOUT_LOGGER");
                auto parseContext =
                    std::make_unique<Ewoms::ParseContext>(std::vector<std::pair<std::string , InputError::Action>>
                                                        {{Ewoms::ParseContext::PARSE_RANDOM_SLASH, Ewoms::InputError::IGNORE},
                                                         {Ewoms::ParseContext::PARSE_MISSING_DIMS_KEYWORD, Ewoms::InputError::WARN},
                                                         {Ewoms::ParseContext::SUMMARY_UNKNOWN_WELL, Ewoms::InputError::WARN},
                                                         {Ewoms::ParseContext::SUMMARY_UNKNOWN_GROUP, Ewoms::InputError::WARN}});
                if (EWOMS_GET_PARAM(PreTypeTag, bool, EclStrictParsing))
                    parseContext->update(Ewoms::InputError::DELAYED_EXIT1);

                Ewoms::EFlowMain<PreTypeTag>::printPRTHeader(outputCout_);

                if (outputCout_) {
                    OpmLog::info("Reading deck file '" + deckFilename + "'");
                }

                readDeck(mpiRank, deckFilename, deck_, eclipseState_, schedule_,
                         summaryConfig_, nullptr, std::move(parseContext),
                         init_from_restart_file, outputCout_);

                if (outputCout_) {
                    OpmLog::info("Done reading deck file.");
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

        int argc_;
        char** argv_;
        bool outputCout_;
        bool outputFiles_;
        double setupTime_;
        std::string deckFilename_;
        std::string flowProgName_;
        char *saveArgs_[2];
        std::unique_ptr<Ewoms::Deck> deck_;
        std::unique_ptr<Ewoms::EclipseState> eclipseState_;
        std::unique_ptr<Ewoms::Schedule> schedule_;
        std::unique_ptr<Ewoms::SummaryConfig> summaryConfig_;
    };

} // namespace Ewoms

#endif // EWOMS_MAIN_HH
