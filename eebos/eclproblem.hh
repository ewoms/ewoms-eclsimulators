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
 *
 * \copydoc Ewoms::EclProblem
 */
#ifndef EWOMS_ECL_PROBLEM_HH
#define EWOMS_ECL_PROBLEM_HH

#include "eclcpgridvanguard.hh"

#include "eclwriter.hh"
#include "eclwellmanager.hh"
#include "eclequilinitializer.hh"
#include "ecloutputblackoilmodule.hh"
#include "ecltransmissibility.hh"
#include "eclthresholdpressure.hh"
#include "ecldummygradientcalculator.hh"
#include "eclfluxmodule.hh"
#include "eclbaseaquifermodel.hh"
#include "eclnewtonmethod.hh"
#include "ecltracermodel.hh"
#include "vtkecltracermodule.hh"

#include <ewoms/numerics/utils/pffgridvector.hh>
#include <ewoms/numerics/models/blackoil/blackoilmodel.hh>
#include <ewoms/numerics/discretizations/ecfv/ecfvdiscretization.hh>

#include <ewoms/material/fluidmatrixinteractions/eclmateriallawmanager.hh>
#include <ewoms/material/thermal/eclthermallawmanager.hh>

#include <ewoms/material/fluidstates/compositionalfluidstate.hh>
#include <ewoms/material/fluidsystems/blackoilfluidsystem.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/drygaspvt.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/wetgaspvt.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/liveoilpvt.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/deadoilpvt.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/constantcompressibilityoilpvt.hh>
#include <ewoms/material/fluidsystems/blackoilpvt/constantcompressibilitywaterpvt.hh>
#include <ewoms/common/intervaltabulated2dfunction.hh>
#include <ewoms/common/uniformxtabulated2dfunction.hh>
#include <ewoms/common/tabulated1dfunction.hh>

#include <ewoms/common/valgrind.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/simulationconfig/rockconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/eqldims.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/actioncontext.hh>
#include <ewoms/eclio/parser/eclipsestate/summaryconfig/summaryconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rockwnodtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/overburdtable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/rocktabtable.hh>
#include <ewoms/common/conditionalstorage.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <ewoms/eclio/output/eclipseio.hh>

#include <ewoms/eclio/opmlog/opmlog.hh>

#include <boost/date_time.hpp>

#include <set>
#include <vector>
#include <string>
#include <algorithm>

namespace Ewoms {
template <class TypeTag>
class EclProblem;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclBaseProblem, INHERITS_FROM(EclCpGridVanguard, EclOutputBlackOil, VtkEclTracer));
//NEW_TYPE_TAG(EclBaseProblem, INHERITS_FROM(EclPolyhedralGridVanguard, EclOutputBlackOil, VtkEclTracer));

// The class which deals with ECL wells
NEW_PROP_TAG(EclWellModel);

// Write all solutions for visualization, not just the ones for the
// report steps...
NEW_PROP_TAG(EnableWriteAllSolutions);

// The number of time steps skipped between writing two consequtive restart files
NEW_PROP_TAG(RestartWritingInterval);

// Enable partial compensation of systematic mass losses via the source term of the next time
// step
NEW_PROP_TAG(EclEnableDriftCompensation);

// Enable the additional checks even if compiled in debug mode (i.e., with the NDEBUG
// macro undefined). Next to a slightly better performance, this also eliminates some
// print statements in debug mode.
NEW_PROP_TAG(EnableDebuggingChecks);

// if thermal flux boundaries are enabled an effort is made to preserve the initial
// thermal gradient specified via the TEMPVD keyword
NEW_PROP_TAG(EnableThermalFluxBoundaries);

// Specify whether API tracking should be enabled (replaces PVT regions).
// TODO: This is not yet implemented
NEW_PROP_TAG(EnableApiTracking);

// The class which deals with ECL aquifers
NEW_PROP_TAG(EclAquiferModel);

// In experimental mode, decides if the aquifer model should be enabled or not
NEW_PROP_TAG(EclEnableAquifers);

// Global switch to enable or disable all ECL-style simulation output
NEW_PROP_TAG(EnableEclOutput);

// Switch to determining whether to write ECL-style asynchronously
//
// I.e., in a non-blocking way. Writing asynsynchronously may lead to race conditions and
// other "fun" issues for debugging and only really makes a noticeable difference if the
// mass storage to be written to is slow (e.g., if it is on a network file system), but
// the code for this feature is quite mature, so it is enabled by default.
NEW_PROP_TAG(EnableAsyncEclOutput);

// Determines whether to write ECL output using double or single-precision floating point
// value.
NEW_PROP_TAG(EclOutputDoublePrecision);

// Global switch to enable or disable all ECL-style simulation output
NEW_PROP_TAG(EnableEclOutput);

// Switch to determining whether to write ECL-style asynchronously
//
// I.e., in a non-blocking way. Writing asynsynchronously may lead to race conditions and
// other "fun" issues for debugging and only really makes a noticeable difference if the
// mass storage to be written to is slow (e.g., if it is on a network file system), but
// the code for this feature is quite mature, so it is enabled by default.
NEW_PROP_TAG(EnableAsyncEclOutput);

// Determines whether to write ECL output using double or single-precision floating point
// value.
NEW_PROP_TAG(EclOutputDoublePrecision);

// time stepping parameters
NEW_PROP_TAG(EclMaxTimeStepSizeAfterWellEvent);
NEW_PROP_TAG(EclRestartShrinkFactor);
NEW_PROP_TAG(EclEnableTuning);
NEW_PROP_TAG(OutputMode);

// Set the problem property
SET_TYPE_PROP(EclBaseProblem, Problem, Ewoms::EclProblem<TypeTag>);

// Select the element centered finite volume method as spatial discretization
SET_TAG_PROP(EclBaseProblem, SpatialDiscretizationSplice, EcfvDiscretization);

//! for eebos, use automatic differentiation to linearize the system of PDEs
SET_TAG_PROP(EclBaseProblem, LocalLinearizerSplice, AutoDiffLocalLinearizer);

// Set the material law for fluid fluxes
SET_PROP(EclBaseProblem, MaterialLaw)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);

    typedef Ewoms::ThreePhaseMaterialTraits<Scalar,
                                          /*wettingPhaseIdx=*/FluidSystem::waterPhaseIdx,
                                          /*nonWettingPhaseIdx=*/FluidSystem::oilPhaseIdx,
                                          /*gasPhaseIdx=*/FluidSystem::gasPhaseIdx> Traits;

public:
    typedef Ewoms::EclMaterialLawManager<Traits> EclMaterialLawManager;

    typedef typename EclMaterialLawManager::MaterialLaw type;
};

// Set the material law for energy storage in rock
SET_PROP(EclBaseProblem, SolidEnergyLaw)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    typedef Ewoms::EclThermalLawManager<Scalar, FluidSystem> EclThermalLawManager;

    typedef typename EclThermalLawManager::SolidEnergyLaw type;
};

// Set the material law for thermal conduction
SET_PROP(EclBaseProblem, ThermalConductionLaw)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    typedef Ewoms::EclThermalLawManager<Scalar, FluidSystem> EclThermalLawManager;

    typedef typename EclThermalLawManager::ThermalConductionLaw type;
};

// eebos can use a slightly faster stencil class because it does not need the normals and
// the integration points of intersections
SET_PROP(EclBaseProblem, Stencil)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);

public:
    typedef Ewoms::EcfvStencil<Scalar,
                             GridView,
                             /*needIntegrationPos=*/false,
                             /*needNormal=*/false> type;
};

// by default use the dummy aquifer "model"
SET_TYPE_PROP(EclBaseProblem, EclAquiferModel, Ewoms::EclBaseAquiferModel<TypeTag>);

// use the built-in proof of concept well model by default
SET_TYPE_PROP(EclBaseProblem, EclWellModel, Ewoms::EclWellManager<TypeTag>);

// Enable aquifers by default in experimental mode
SET_BOOL_PROP(EclBaseProblem, EclEnableAquifers, true);

// Enable gravity
SET_BOOL_PROP(EclBaseProblem, EnableGravity, true);

// only write the solutions for the report steps to disk
SET_BOOL_PROP(EclBaseProblem, EnableWriteAllSolutions, false);

// disable API tracking
SET_BOOL_PROP(EclBaseProblem, EnableApiTracking, false);

// The default for the end time of the simulation [s]
//
// By default, stop it after the universe will probably have stopped
// to exist. (the ECL problem will finish the simulation explicitly
// after it simulated the last episode specified in the deck.)
SET_SCALAR_PROP(EclBaseProblem, EndTime, 1e100);

// The default for the initial time step size of the simulation [s].
//
// The chosen value means that the size of the first time step is the
// one of the initial episode (if the length of the initial episode is
// not millions of trillions of years, that is...)
SET_SCALAR_PROP(EclBaseProblem, InitialTimeStepSize, 3600.0*24);

// the default for the allowed volumetric error for oil per second [kg/(m^3 s)]
SET_SCALAR_PROP(EclBaseProblem, NewtonTolerance, 10.0);

// the tolerated amount of "incorrect" amount of oil per second for the complete
// reservoir. this is scaled by the pore volume of the reservoir, i.e., larger reservoirs
// will tolerate larger residuals.
SET_SCALAR_PROP(EclBaseProblem, EclNewtonSumTolerance, 1.0);

// set the exponent for the volume scaling of the sum tolerance: larger reservoirs can
// tolerate a higher amount of mass lost per time step than smaller ones! since this is
// not linear, we use the cube root of the overall pore volume by default, i.e., the
// value specified by the NewtonSumTolerance parameter is the "incorrect" mass per
// timestep for an reservoir that exhibits 1 m^3 of pore volume. A reservoir with a total
// pore volume of 10^3 m^3 will tolerate 10 times as much.
SET_SCALAR_PROP(EclBaseProblem, EclNewtonSumToleranceExponent, 1.0/3.0);

// set number of Newton iterations where the volumetric residual is considered for
// convergence
SET_INT_PROP(EclBaseProblem, EclNewtonStrictIterations, 8);

// set fraction of the pore volume where the volumetric residual may be violated during
// strict Newton iterations
SET_SCALAR_PROP(EclBaseProblem, EclNewtonRelaxedVolumeFraction, 0.03);

// the maximum volumetric error of a cell in the relaxed region
SET_SCALAR_PROP(EclBaseProblem, EclNewtonRelaxedTolerance, 1e9);

// Ignore the maximum error mass for early termination of the newton method.
SET_SCALAR_PROP(EclBaseProblem, NewtonMaxError, 10e9);

// set the maximum number of Newton iterations to 14 because the likelyhood that a time
// step succeeds at more than 14 Newton iteration is rather small
SET_INT_PROP(EclBaseProblem, NewtonMaxIterations, 14);

// also, reduce the target for the "optimum" number of Newton iterations to 6. Note that
// this is only relevant if the time step is reduced from the report step size for some
// reason. (because eebos first tries to do a report step using a single time step.)
SET_INT_PROP(EclBaseProblem, NewtonTargetIterations, 6);

// Disable the VTK output by default for this problem ...
SET_BOOL_PROP(EclBaseProblem, EnableVtkOutput, false);

// ... but enable the ECL output by default
SET_BOOL_PROP(EclBaseProblem, EnableEclOutput, true);

// If available, write the ECL output in a non-blocking manner
SET_BOOL_PROP(EclBaseProblem, EnableAsyncEclOutput, true);

// By default, use single precision for the ECL formated results
SET_BOOL_PROP(EclBaseProblem, EclOutputDoublePrecision, false);

// The default location for the ECL output files
SET_STRING_PROP(EclBaseProblem, OutputDir, ".");

// the cache for intensive quantities can be used for ECL problems and also yields a
// decent speedup...
SET_BOOL_PROP(EclBaseProblem, EnableIntensiveQuantityCache, true);

// the cache for the storage term can also be used and also yields a decent speedup
SET_BOOL_PROP(EclBaseProblem, EnableStorageCache, true);

// Use the "velocity module" which uses the Eclipse "NEWTRAN" transmissibilities
SET_TYPE_PROP(EclBaseProblem, FluxModule, Ewoms::EclTransFluxModule<TypeTag>);

// Use the dummy gradient calculator in order not to do unnecessary work.
SET_TYPE_PROP(EclBaseProblem, GradientCalculator, Ewoms::EclDummyGradientCalculator<TypeTag>);

// Use a custom Newton-Raphson method class for eebos in order to attain more
// sophisticated update and error computation mechanisms
SET_TYPE_PROP(EclBaseProblem, NewtonMethod, Ewoms::EclNewtonMethod<TypeTag>);

// The frequency of writing restart (*.ers) files. This is the number of time steps
// between writing restart files
SET_INT_PROP(EclBaseProblem, RestartWritingInterval, 0xffffff); // disable

// Drift compensation is an experimental feature, i.e., systematic errors in the
// conservation quantities are only compensated for
// as default if experimental mode is enabled.
SET_BOOL_PROP(EclBaseProblem, EclEnableDriftCompensation, GET_PROP_VALUE(TypeTag, EnableExperiments));

// By default, we enable the debugging checks if we're compiled in debug mode
SET_BOOL_PROP(EclBaseProblem, EnableDebuggingChecks, true);

// store temperature (but do not conserve energy, as long as EnableEnergy is false)
SET_BOOL_PROP(EclBaseProblem, EnableTemperature, true);

// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
SET_BOOL_PROP(EclBaseProblem, EnablePolymer, false);
SET_BOOL_PROP(EclBaseProblem, EnableSolvent, false);
SET_BOOL_PROP(EclBaseProblem, EnableEnergy, false);
SET_BOOL_PROP(EclBaseProblem, EnableFoam, false);

// disable thermal flux boundaries by default
SET_BOOL_PROP(EclBaseProblem, EnableThermalFluxBoundaries, false);

SET_BOOL_PROP(EclBaseProblem, EnableTracerModel, false);

// By default, simulators derived from the EclBaseProblem are production simulators,
// i.e., experimental features must be explicitly enabled at compile time
SET_BOOL_PROP(EclBaseProblem, EnableExperiments, false);

// set defaults for the time stepping parameters
SET_SCALAR_PROP(EclBaseProblem, EclMaxTimeStepSizeAfterWellEvent, 3600*24*365.25);
SET_SCALAR_PROP(EclBaseProblem, EclRestartShrinkFactor, 3);
SET_BOOL_PROP(EclBaseProblem, EclEnableTuning, false);

SET_STRING_PROP(EclBaseProblem, OutputMode, "all");

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief This problem simulates an input file given in the data format used by the
 *        commercial ECLIPSE simulator.
 */
template <class TypeTag>
class EclProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    using ParentType = GET_PROP_TYPE(TypeTag, BaseProblem);
    using Implementation = GET_PROP_TYPE(TypeTag, Problem);

    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);
    using Stencil = GET_PROP_TYPE(TypeTag, Stencil);
    using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
    using GlobalEqVector = GET_PROP_TYPE(TypeTag, GlobalEqVector);
    using EqVector = GET_PROP_TYPE(TypeTag, EqVector);

    // Grid and world dimension
    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    enum { numPhases = FluidSystem::numPhases };
    enum { enableExperiments = GET_PROP_VALUE(TypeTag, EnableExperiments) };
    enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
    enum { enablePolymer = GET_PROP_VALUE(TypeTag, EnablePolymer) };
    enum { enableBrine = GET_PROP_VALUE(TypeTag, EnableBrine) };
    enum { enablePolymerMolarWeight = GET_PROP_VALUE(TypeTag, EnablePolymerMW) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableThermalFluxBoundaries = GET_PROP_VALUE(TypeTag, EnableThermalFluxBoundaries) };
    enum { enableApiTracking = GET_PROP_VALUE(TypeTag, EnableApiTracking) };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { waterCompIdx = FluidSystem::waterCompIdx };

    using PrimaryVariables = GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);
    using BoundaryRateVector = GET_PROP_TYPE(TypeTag, BoundaryRateVector);
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using EclMaterialLawManager = typename GET_PROP(TypeTag, MaterialLaw)::EclMaterialLawManager;
    using EclThermalLawManager = typename GET_PROP(TypeTag, SolidEnergyLaw)::EclThermalLawManager;
    using MaterialLawParams = typename EclMaterialLawManager::MaterialLawParams;
    using SolidEnergyLawParams = typename EclThermalLawManager::SolidEnergyLawParams;
    using ThermalConductionLawParams = typename EclThermalLawManager::ThermalConductionLawParams;
    using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
    using DofMapper = GET_PROP_TYPE(TypeTag, DofMapper);
    using Evaluation = GET_PROP_TYPE(TypeTag, Evaluation);
    using Indices = GET_PROP_TYPE(TypeTag, Indices);
    using IntensiveQuantities = GET_PROP_TYPE(TypeTag, IntensiveQuantities);
    using EclWellModel = GET_PROP_TYPE(TypeTag, EclWellModel);
    using EclAquiferModel = GET_PROP_TYPE(TypeTag, EclAquiferModel);

    typedef BlackOilSolventModule<TypeTag> SolventModule;
    typedef BlackOilPolymerModule<TypeTag> PolymerModule;
    typedef BlackOilFoamModule<TypeTag> FoamModule;
    typedef BlackOilBrineModule<TypeTag> BrineModule;

    typedef typename EclEquilInitializer<TypeTag>::ScalarFluidState InitialFluidState;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef EclWriter<TypeTag> EclWriterType;
    typedef EclTracerModel<TypeTag> TracerModel;
    typedef Ewoms::UniformXTabulated2DFunction<Scalar> TabulatedTwoDFunction;
    typedef Ewoms::Tabulated1DFunction<Scalar> TabulatedFunction;

    struct RockParams {
        Scalar referencePressure;
        Scalar compressibility;
    };

public:
    /*!
     * \copydoc FvBaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();
        EclWriterType::registerParameters();
        VtkEclTracerModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableWriteAllSolutions,
                             "Write all solutions to disk instead of only the ones for the "
                             "report steps");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableEclOutput,
                             "Write binary output which is compatible with the commercial "
                             "Eclipse simulator");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclOutputDoublePrecision,
                             "Tell the output writer to use double precision. Useful for 'perfect' restarts");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, RestartWritingInterval,
                             "The frequencies of which time steps are serialized to disk");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableTracerModel,
                             "Transport tracers found in the deck.");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclEnableDriftCompensation,
                             "Enable partial compensation of systematic mass losses via the source term of the next time step");
        if (enableExperiments)
            EWOMS_REGISTER_PARAM(TypeTag, bool, EclEnableAquifers,
                                 "Enable analytic and numeric aquifer models");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclMaxTimeStepSizeAfterWellEvent,
                             "Maximum time step size after an well event");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclRestartShrinkFactor,
                             "Factor by which the time step is reduced after convergence failure");
        EWOMS_REGISTER_PARAM(TypeTag, bool, EclEnableTuning,
                             "Honor some aspects of the TUNING keyword from the ECL deck.");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, OutputMode,
                             "Specify which messages are going to be printed. Valid values are: none, log, all (default)");

    }

    /*!
     * \copydoc FvBaseProblem::prepareOutputDir
     */
    std::string prepareOutputDir() const
    { return this->simulator().vanguard().eclState().getIOConfig().getOutputDir(); }

    /*!
     * \copydoc FvBaseProblem::handlePositionalParameter
     */
    static int handlePositionalParameter(std::set<std::string>& seenParams,
                                         std::string& errorMsg,
                                         int argc EWOMS_UNUSED,
                                         const char** argv,
                                         int paramIdx,
                                         int posParamIdx EWOMS_UNUSED)
    {
        using ParamsMeta = GET_PROP(TypeTag, ParameterMetaData);
        Dune::ParameterTree& tree = ParamsMeta::tree();

        std::string param = argv[paramIdx];
        size_t i = param.find('=');
        if (i != std::string::npos) {
            std::string oldParamName = param.substr(0, i);
            std::string oldParamValue = param.substr(i+1);
            std::string newParamName = "--" + oldParamName;
            for (size_t j = 0; j < newParamName.size(); ++j)
                if (newParamName[j] == '_')
                    newParamName[j] = '-';
            errorMsg =
                "The old syntax to specify parameters on the command line is no longer supported: "
                "Try replacing '"+oldParamName+"="+oldParamValue+"' with "+
                "'"+newParamName+"="+oldParamValue+"'!";
            return 0;
        }

        if (seenParams.count("EclDeckFileName") > 0) {
            errorMsg =
                "Parameter 'EclDeckFileName' specified multiple times"
                " as a command line parameter";
            return 0;
        }

        tree["EclDeckFileName"] = argv[paramIdx];
        seenParams.insert("EclDeckFileName");
        return 1;
    }

    /*!
     * \copydoc FvBaseProblem::helpPreamble
     */
    static std::string helpPreamble(int argc EWOMS_UNUSED,
                                    const char **argv)
    {
        std::string desc = Implementation::briefDescription();
        if (!desc.empty())
            desc = desc + "\n";

        return
            "Usage: "+std::string(argv[0]) + " [OPTIONS] [ECL_DECK_FILENAME]\n"
            + desc;
    }

    /*!
     * \copydoc FvBaseProblem::briefDescription
     */
    static std::string briefDescription()
    {
        if (briefDescription_.empty()) {
            // return the default description. this is different for releases and
            // development versions.
            std::string result =
                "The Ewoms Ecl Black-Oil Simulator (eebos); an advanced "
                "ECL-compatible simulation program for hydrocarbon reservoirs "
                "developed by the eWoms project (https://ewoms.org).";

            bool isDevelVersion =
                std::string(EWOMS_COMMON_VERSION).find("pre") != std::string::npos ||
                std::string(EWOMS_MATERIAL_VERSION).find("pre") != std::string::npos ||
                std::string(EWOMS_NUMERICS_VERSION).find("pre") != std::string::npos ||
                std::string(EWOMS_ECLIO_VERSION).find("pre") != std::string::npos ||
                std::string(EWOMS_ECLGRIDS_VERSION).find("pre") != std::string::npos ||
                std::string(EWOMS_ECLSIMULATORS_VERSION).find("pre") != std::string::npos;

            if (isDevelVersion)
                result +=
                    "\n\n"
                    "YOU ARE USING A DEVELOPMENT VERSION OF THE SIMULATOR. Do not use it for production runs!";

            return result;
        }
        else
            return briefDescription_;
    }

    /*!
     * \brief Specifies the string returned by briefDescription()
     *
     * This string appears in the usage message.
     */
    static void setBriefDescription(const std::string& msg)
    { briefDescription_ = msg; }

    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    EclProblem(Simulator& simulator)
        : ParentType(simulator)
        , transmissibilities_(simulator.vanguard())
        , thresholdPressures_(simulator)
        , wellModel_(simulator)
        , aquiferModel_(simulator)
        , pffDofData_(simulator.gridView(), this->elementMapper())
        , tracerModel_(simulator)
    {
        this->model().addOutputModule(new VtkEclTracerModule<TypeTag>(simulator));
        // Tell the black-oil extensions to initialize their internal data structures
        const auto& vanguard = simulator.vanguard();
        SolventModule::initFromEclState(vanguard.eclState(), vanguard.schedule());
        PolymerModule::initFromEclState(vanguard.eclState());
        FoamModule::initFromEclState(vanguard.eclState());
        BrineModule::initFromEclState(vanguard.eclState());

        // create the ECL writer
        eclWriter_.reset(new EclWriterType(simulator));

        enableDriftCompensation_ = EWOMS_GET_PARAM(TypeTag, bool, EclEnableDriftCompensation);

        enableEclOutput_ = EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput);

        if (enableExperiments)
            enableAquifers_ = EWOMS_GET_PARAM(TypeTag, bool, EclEnableAquifers);
        else
            enableAquifers_ = true;

        enableTuning_ = EWOMS_GET_PARAM(TypeTag, bool, EclEnableTuning);
        initialTimeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, InitialTimeStepSize);
        minTimeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, MinTimeStepSize);
        maxTimeStepSize_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxTimeStepSize);
        maxTimeStepAfterWellEvent_ = EWOMS_GET_PARAM(TypeTag, Scalar, EclMaxTimeStepSizeAfterWellEvent);
        restartShrinkFactor_ = EWOMS_GET_PARAM(TypeTag, Scalar, EclRestartShrinkFactor);
        maxFails_ = EWOMS_GET_PARAM(TypeTag, unsigned, MaxTimeStepDivisions);
    }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        auto& simulator = this->simulator();
        const auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();
        const auto& timeMap = schedule.getTimeMap();

        // Set the start time of the simulation
        simulator.setStartTime(timeMap.getStartTime(/*reportStepIdx=*/0));
        simulator.setEndTime(timeMap.getTotalTime());

        // We want the episode index to be the same as the report step index to make
        // things simpler, so we have to set the episode index to -1 because it is
        // incremented by endEpisode(). The size of the initial time step and
        // length of the initial episode is set to zero for the same reason.
        simulator.setEpisodeIndex(-1);
        simulator.setEpisodeLength(0.0);

        // the "NOGRAV" keyword from Frontsim or setting the EnableGravity to false
        // disables gravity, else the standard value of the gravity constant at sea level
        // on earth is used
        this->gravity_ = 0.0;
        if (EWOMS_GET_PARAM(TypeTag, bool, EnableGravity))
            this->gravity_[dim - 1] = 9.80665;
        if (!eclState.getInitConfig().hasGravity())
            this->gravity_[dim - 1] = 0.0;

        if (enableTuning_) {
            // if support for the TUNING keyword is enabled, we get the initial time
            // steping parameters from it instead of from command line parameters
            const auto& tuning = schedule.getTuning(0);
            initialTimeStepSize_ = tuning.TSINIT;
            maxTimeStepAfterWellEvent_ = tuning.TMAXWC;
            maxTimeStepSize_ = tuning.TSMAXZ;
            restartShrinkFactor_ = 1./tuning.TSFCNV;
            minTimeStepSize_ = tuning.TSMINZ;
        }

        initFluidSystem_();

        // deal with DRSDT
        unsigned ntpvt = eclState.runspec().tabdims().getNumPVTTables();
        size_t numDof = this->model().numGridDof();
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            maxDRs_.resize(ntpvt, 1e30);
            dRsDtOnlyFreeGas_.resize(ntpvt, false);
            lastRs_.resize(numDof, 0.0);
            maxDRv_.resize(ntpvt, 1e30);
            lastRv_.resize(numDof, 0.0);
            maxOilSaturation_.resize(numDof, 0.0);
        }

        updateElementDepths_();
        readRockParameters_();
        readMaterialParameters_();
        readThermalParameters_();
        transmissibilities_.finishInit();

        const auto& initconfig = eclState.getInitConfig();
        if (initconfig.restartRequested())
            readEclRestartSolution_();
        else
            readInitialCondition_();

        updatePffDofData_();

        if (GET_PROP_VALUE(TypeTag, EnablePolymer)) {
            const auto& vanguard = this->simulator().vanguard();
            const auto& gridView = vanguard.gridView();
            int numElements = gridView.size(/*codim=*/0);
            maxPolymerAdsorption_.resize(numElements, 0.0);
        }

        tracerModel_.init();

        readBoundaryConditions_();

        if (enableDriftCompensation_) {
            drift_.resize(numDof);
            drift_ = 0.0;
        }

        if (enableExperiments && this->gridView().comm().rank() == 0)
            checkDeckCompatibility_();

        // write the static output files (EGRID, INIT, SMSPEC, etc.)
        if (enableEclOutput_)
            eclWriter_->writeInit();

        simulator.vanguard().releaseGlobalTransmissibilities();

        // after finishing the initialization and writing the initial solution, we move
        // to the first "real" episode/report step
        // for restart the episode index and start is already set
        if (!initconfig.restartRequested()) {
            simulator.startNextEpisode(timeMap.getTimeStepLength(0));
            simulator.setEpisodeIndex(0);
        }
    }

    void prefetch(const Element& elem) const
    { pffDofData_.prefetch(elem); }

    /*!
     * \brief This method restores the complete state of the problem and its sub-objects
     *        from disk.
     *
     * The serialization format used by this method is ad-hoc. It is the inverse of the
     * serialize() method.
     *
     * \tparam Restarter The deserializer type
     *
     * \param res The deserializer object
     */
    template <class Restarter>
    void deserialize(Restarter& res)
    {
        // reload the current episode/report step from the deck
        beginEpisode();

        // deserialize the wells
        wellModel_.deserialize(res);

        if (enableAquifers_)
            // deserialize the aquifer
            aquiferModel_.deserialize(res);
    }

    /*!
     * \brief This method writes the complete state of the problem and its subobjects to
     *        disk.
     *
     * The file format used here is ad-hoc.
     */
    template <class Restarter>
    void serialize(Restarter& res)
    {
        wellModel_.serialize(res);

        if (enableAquifers_)
            aquiferModel_.serialize(res);
    }

    /*!
     * \brief Called by the simulator before an episode begins.
     */
    void beginEpisode()
    {
        // Proceed to the next report step
        auto& simulator = this->simulator();
        int episodeIdx = simulator.episodeIndex();
        auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();
        const auto& events = schedule.getEvents();
        const auto& timeMap = schedule.getTimeMap();

        if (episodeIdx >= 0 && events.hasEvent(Ewoms::ScheduleEvents::GEO_MODIFIER, episodeIdx)) {
            // bring the contents of the keywords to the current state of the SCHEDULE
            // section.
            //
            // TODO (?): make grid topology changes possible (depending on what exactly
            // has changed, the grid may need be re-created which has some serious
            // implications on e.g., the solution of the simulation.)
            const auto& miniDeck = schedule.getModifierDeck(episodeIdx);
            eclState.applyModifierDeck(miniDeck);

            // re-compute all quantities which may possibly be affected.
            transmissibilities_.update(true);
            referencePorosity_[1] = referencePorosity_[0];
            updateReferencePorosity_();
            updatePffDofData_();
        }

        if (enableExperiments && this->gridView().comm().rank() == 0 && episodeIdx >= 0) {
            // print some useful information in experimental mode. (the production
            // simulator does this externally.)
            std::ostringstream ss;
            boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
            boost::posix_time::ptime curDateTime =
                boost::posix_time::from_time_t(timeMap.getStartTime(episodeIdx));
            ss.imbue(std::locale(std::locale::classic(), facet));
            ss << "Report step " << episodeIdx + 1
                      << "/" << timeMap.numTimesteps()
                      << " at day " << timeMap.getTimePassedUntil(episodeIdx)/(24*3600)
                      << "/" << timeMap.getTotalTime()/(24*3600)
                      << ", date = " << curDateTime.date()
                      << "\n ";
            OpmLog::info(ss.str());
        }

        // react to TUNING changes
        bool tuningEvent = false;
        if (episodeIdx > 0 && enableTuning_ && events.hasEvent(Ewoms::ScheduleEvents::TUNING_CHANGE, episodeIdx))
        {
            const auto& tuning = schedule.getTuning(episodeIdx);
            initialTimeStepSize_ = tuning.TSINIT;
            maxTimeStepAfterWellEvent_ = tuning.TMAXWC;
            maxTimeStepSize_ = tuning.TSMAXZ;
            restartShrinkFactor_ = 1./tuning.TSFCNV;
            minTimeStepSize_ = tuning.TSMINZ;
            tuningEvent = true;
        }

        // set up the wells for the next episode.
        wellModel_.beginEpisode();

        // set up the aquifers for the next episode.
        if (enableAquifers_)
            // set up the aquifers for the next episode.
            aquiferModel_.beginEpisode();

        // set the size of the initial time step of the episode
        Scalar dt = limitNextTimeStepSize_(simulator.episodeLength());
        if (episodeIdx == 0 || tuningEvent)
            // allow the size of the initial time step to be set via an external parameter
            // if TUNING is enabled, also limit the time step size after a tuning event to TSINIT
            dt = std::min(dt, initialTimeStepSize_);
        simulator.setTimeStepSize(dt);
    }

    /*!
     * \brief Called by the simulator before each time integration.
     */
    void beginTimeStep()
    {
        const auto& simulator = this->simulator();
        int episodeIdx = simulator.episodeIndex();

        if (enableExperiments && this->gridView().comm().rank() == 0 && episodeIdx >= 0) {
            std::ostringstream ss;
            boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
            boost::posix_time::ptime date = boost::posix_time::from_time_t((this->simulator().startTime())) + boost::posix_time::milliseconds(static_cast<long long>(this->simulator().time() / Ewoms::prefix::milli));
            ss.imbue(std::locale(std::locale::classic(), facet));
            ss <<"Eebos: Begin time step " << this->simulator().timeStepIndex() << ", stepsize "
               << unit::convert::to(this->simulator().timeStepSize(), unit::day) << " days,"
               << " at day " << (double)unit::convert::to(this->simulator().time(), unit::day)
               << "/" << (double)unit::convert::to(this->simulator().endTime(), unit::day)
               << ", date = " << date;
            OpmLog::info(ss.str());
        }

        // update explicit quantities between timesteps.
        const auto& oilVaporizationControl = simulator.vanguard().schedule().getOilVaporizationProperties(episodeIdx);
        if (drsdtActive_())
            // DRSDT is enabled
            for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRs_.size(); ++pvtRegionIdx)
                maxDRs_[pvtRegionIdx] = oilVaporizationControl.getMaxDRSDT(pvtRegionIdx)*simulator.timeStepSize();

        if (drvdtActive_())
            // DRVDT is enabled
            for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRv_.size(); ++pvtRegionIdx)
                maxDRv_[pvtRegionIdx] = oilVaporizationControl.getMaxDRVDT(pvtRegionIdx)*this->simulator().timeStepSize();

        // update maximum water saturation and minimum pressure
        // used when ROCKCOMP is activated
        const bool invalidateFromMaxWaterSat = updateMaxWaterSaturation_();
        const bool invalidateFromMinPressure = updateMinPressure_();

        // update hysteresis and max oil saturation used in vappars
        const bool invalidateFromHyst = updateHysteresis_();
        const bool invalidateFromMaxOilSat = updateMaxOilSaturation_();

        // the derivatives may have change
        bool invalidateIntensiveQuantities = invalidateFromMaxWaterSat || invalidateFromMinPressure || invalidateFromHyst || invalidateFromMaxOilSat;
        if (invalidateIntensiveQuantities)
            this->model().updateIntensiveQuantitiesCache(/*timeIdx=*/0);

        if (GET_PROP_VALUE(TypeTag, EnablePolymer))
            updateMaxPolymerAdsorption_();

        wellModel_.beginTimeStep();
        if (enableAquifers_)
            aquiferModel_.beginTimeStep();
        tracerModel_.beginTimeStep();

    }

    /*!
     * \brief Return if the storage term of the first iteration is identical to the storage
     *        term for the solution of the previous time step.
     *
     * For quite technical reasons, the storage term cannot be recycled if either DRSDT
     * or DRVDT are active in eebos. Nor if the porosity is changes between timesteps
     * using a pore volume multiplier (i.e., poreVolumeMultiplier() != 1.0)
     */
    bool recycleFirstIterationStorage() const
    { return !drsdtActive_() && !drvdtActive_() && rockCompPoroMultWc_.empty() && rockCompPoroMult_.empty();  }

    /*!
     * \brief Called by the simulator before each Newton-Raphson iteration.
     */
    void beginIteration()
    {
        wellModel_.beginIteration();
        if (enableAquifers_)
            aquiferModel_.beginIteration();
    }

    /*!
     * \brief Called by the simulator after each Newton-Raphson iteration.
     */
    void endIteration()
    {
        wellModel_.endIteration();
        if (enableAquifers_)
            aquiferModel_.endIteration();
    }

    /*!
     * \brief Called by the simulator after each time integration.
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        if (GET_PROP_VALUE(TypeTag, EnableDebuggingChecks)) {
            // in debug mode, we don't care about performance, so we check if the model does
            // the right thing (i.e., the mass change inside the whole reservoir must be
            // equivalent to the fluxes over the grid's boundaries plus the source rates
            // specified by the problem)
            int rank = this->simulator().gridView().comm().rank();
            if (rank == 0)
                std::cout << "checking conservativeness of solution\n";
            this->model().checkConservativeness(/*tolerance=*/-1, /*verbose=*/true);
            if (rank == 0)
                std::cout << "solution is sufficiently conservative\n";
        }
#endif // NDEBUG

        auto& simulator = this->simulator();
        wellModel_.endTimeStep();
        if (enableAquifers_)
            aquiferModel_.endTimeStep();
        tracerModel_.endTimeStep();

        // deal with DRSDT and DRVDT
        updateCompositionChangeLimits_();

        if (enableDriftCompensation_) {
            const auto& residual = this->model().linearizer().residual();
            for (unsigned globalDofIdx = 0; globalDofIdx < residual.size(); globalDofIdx ++) {
                drift_[globalDofIdx] = residual[globalDofIdx];
                drift_[globalDofIdx] *= simulator.timeStepSize();
                if (GET_PROP_VALUE(TypeTag, UseVolumetricResidual))
                    drift_[globalDofIdx] *= this->model().dofTotalVolume(globalDofIdx);
            }
        }

        simulator.model().updateWaterPriority();

        bool isSubStep = !EWOMS_GET_PARAM(TypeTag, bool, EnableWriteAllSolutions) && !this->simulator().episodeWillBeOver();
        eclWriter_->evalSummaryState(isSubStep);

        applyActions();
    }

    /*!
     * \brief Called by the simulator after the end of an episode.
     */
    void endEpisode()
    {
        auto& simulator = this->simulator();
        auto& schedule = simulator.vanguard().schedule();
        const auto& timeMap = schedule.getTimeMap();

        int episodeIdx = simulator.episodeIndex();

        // check if we're finished ...
        if (episodeIdx + 1 >= static_cast<int>(timeMap.numTimesteps())) {
            simulator.setFinished(true);
            return;
        }

        // .. if we're not yet done, start the next episode (report step)
        simulator.startNextEpisode(timeMap.getTimeStepLength(episodeIdx + 1));
    }

    /*!
     * \brief Always returns true. The ecl output writer takes care of the rest
     */
    bool shouldWriteOutput() const
    { return true; }

    /*!
     * \brief Returns true if an eWoms restart file should be written to disk.
     *
     * The EclProblem does not write any restart files using the ad-hoc format, only ones
     * using the ECL format.
     */
    bool shouldWriteRestartFile() const
    { return false; }

    /*!
     * \brief Write the requested quantities of the current solution into the output
     *        files.
     */
    void writeOutput(bool verbose = true)
    {
        // use the generic code to prepare the output fields and to
        // write the desired VTK files.
        ParentType::writeOutput(verbose);

        bool isSubStep = !EWOMS_GET_PARAM(TypeTag, bool, EnableWriteAllSolutions) && !this->simulator().episodeWillBeOver();
        if (enableEclOutput_)
            eclWriter_->writeOutput(isSubStep);
    }

    void finalizeOutput() {
        // this will write all pending output to disk
        // to avoid corruption of output files
        eclWriter_.reset();
    }

    void applyActions()
    {
        Simulator& simulator = this->simulator();
        int reportStep = simulator.episodeIndex();
        auto& schedule = simulator.vanguard().schedule();
        const auto& summaryState = simulator.vanguard().summaryState();
        auto& actionState = simulator.vanguard().actionState();
        const auto& actions = schedule.actions(reportStep);
        if (actions.empty())
            return;

        auto scheduleTime = schedule.simTime(reportStep);
        Ewoms::Action::Context context( summaryState, schedule.getWListManager(reportStep) );
        auto now = Ewoms::TimeStampUTC( schedule.getStartTime() ) + std::chrono::duration<double>(scheduleTime);
        std::string ts;
        {
            std::ostringstream os;
            os << std::setw(4) <<                      std::to_string(now.year())  << '/'
               << std::setw(2) << std::setfill('0') << std::to_string(now.month()) << '/'
               << std::setw(2) << std::setfill('0') << std::to_string(now.day()) << "  report:" << std::to_string(reportStep);

            ts = os.str();
        }

        for (const auto& action : actions.pending(actionState, scheduleTime)) {
            const auto& actionResult = action->eval(context);
            if (actionResult) {
                std::string wellsString;
                const auto& matchingWells = actionResult.wells();
                if (matchingWells.size() > 0) {
                    for (std::size_t iw = 0; iw < matchingWells.size() - 1; iw++)
                        wellsString += matchingWells[iw] + ", ";
                    wellsString += matchingWells.back();
                }
                std::string msg = "The action: " + action->name() + " evaluated to true at " + ts + " wells: " + wellsString;
                Ewoms::OpmLog::info(msg);
                schedule.applyAction(reportStep, *action, actionResult);
                actionState.add_run(*action, scheduleTime);
            }
            else {
                std::string msg = "The action: " + action->name() + " evaluated to false at " + ts;
                Ewoms::OpmLog::info(msg);
            }
        }
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context,
                                           unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return transmissibilities_.permeability(globalSpaceIdx);
    }

    /*!
     * \brief This method returns the intrinsic permeability tensor
     *        given a global element index.
     *
     * Its main (only?) usage is the ECL transmissibility calculation code...
     */
    const DimMatrix& intrinsicPermeability(unsigned globalElemIdx) const
    { return transmissibilities_.permeability(globalElemIdx); }

    /*!
     * \copydoc EclTransmissiblity::transmissibility
     */
    template <class Context>
    Scalar transmissibility(const Context& context,
                            unsigned EWOMS_OPTIM_UNUSED fromDofLocalIdx,
                            unsigned toDofLocalIdx) const
    {
        assert(fromDofLocalIdx == 0);
        return pffDofData_.get(context.element(), toDofLocalIdx).transmissibility;
    }

    /*!
     * \copydoc EclTransmissiblity::transmissibilityBoundary
     */
    template <class Context>
    Scalar transmissibilityBoundary(const Context& elemCtx,
                                    unsigned boundaryFaceIdx) const
    {
        unsigned elemIdx = elemCtx.globalSpaceIndex(/*dofIdx=*/0, /*timeIdx=*/0);
        return transmissibilities_.transmissibilityBoundary(elemIdx, boundaryFaceIdx);
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibility(const Context& context,
                                       unsigned faceIdx,
                                       unsigned timeIdx) const
    {
        const auto& face = context.stencil(timeIdx).interiorFace(faceIdx);
        unsigned toDofLocalIdx = face.exteriorIndex();
        return *pffDofData_.get(context.element(), toDofLocalIdx).thermalHalfTrans;
    }

    /*!
     * \copydoc EclTransmissiblity::thermalHalfTransmissibility
     */
    template <class Context>
    Scalar thermalHalfTransmissibilityBoundary(const Context& elemCtx,
                                               unsigned boundaryFaceIdx) const
    {
        unsigned elemIdx = elemCtx.globalSpaceIndex(/*dofIdx=*/0, /*timeIdx=*/0);
        return transmissibilities_.thermalHalfTransBoundary(elemIdx, boundaryFaceIdx);
    }

    /*!
     * \brief Return a reference to the object that handles the "raw" transmissibilities.
     */
    const EclTransmissibility<TypeTag>& eclTransmissibilities() const
    { return transmissibilities_; }

    /*!
     * \copydoc BlackOilBaseProblem::thresholdPressure
     */
    Scalar thresholdPressure(unsigned elem1Idx, unsigned elem2Idx) const
    { return thresholdPressures_.thresholdPressure(elem1Idx, elem2Idx); }

    const EclThresholdPressure<TypeTag>& thresholdPressure() const
    { return thresholdPressures_; }

    EclThresholdPressure<TypeTag>& thresholdPressure()
    { return thresholdPressures_; }

    const EclTracerModel<TypeTag>& tracerModel() const
    { return tracerModel_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     *
     * For the EclProblem, this method is identical to referencePorosity(). The intensive
     * quantities object may apply various multipliers (e.g. ones which model rock
     * compressibility and water induced rock compaction) to it which depend on the
     * current physical conditions.
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return referencePorosity_[timeIdx][globalSpaceIdx];
    }

    /*!
     * \brief Returns the porosity of an element
     *
     * The reference porosity of an element is the porosity of the medium before modified
     * by the current solution. Note that this method is *not* part of the generic eWoms
     * problem API because it would bake the assumption that only the elements are the
     * degrees of freedom into the interface.
     */
    Scalar referencePorosity(unsigned elementIdx, unsigned timeIdx) const
    { return referencePorosity_[timeIdx][elementIdx]; }

    /*!
     * \brief Returns the depth of an degree of freedom [m]
     *
     * For ECL problems this is defined as the average of the depth of an element and is
     * thus slightly different from the depth of an element's centroid.
     */
    template <class Context>
    Scalar dofCenterDepth(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return elementCenterDepth_[globalSpaceIdx];
    }

    /*!
     * \copydoc BlackoilProblem::rockCompressibility
     */
    template <class Context>
    Scalar rockCompressibility(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        if (rockParams_.empty())
            return 0.0;

        unsigned tableIdx = 0;
        if (!rockTableIdx_.empty()) {
            unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            tableIdx = rockTableIdx_[globalSpaceIdx];
        }

        return rockParams_[tableIdx].compressibility;
    }

    /*!
     * \copydoc BlackoilProblem::rockReferencePressure
     */
    template <class Context>
    Scalar rockReferencePressure(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        if (rockParams_.empty())
            return 1e5;

        unsigned tableIdx = 0;
        if (!rockTableIdx_.empty()) {
            unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
            tableIdx = rockTableIdx_[globalSpaceIdx];
        }

        return rockParams_[tableIdx].referencePressure;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return materialLawParams(globalSpaceIdx);
    }

    const MaterialLawParams& materialLawParams(unsigned globalDofIdx) const
    { return materialLawManager_->materialLawParams(globalDofIdx); }

    /*!
     * \brief Return the parameters for the energy storage law of the rock
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context EWOMS_UNUSED,
                         unsigned spaceIdx EWOMS_UNUSED,
                         unsigned timeIdx EWOMS_UNUSED) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return thermalLawManager_->solidEnergyLawParams(globalSpaceIdx);
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionLawParams(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalSpaceIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return thermalLawManager_->thermalConductionLawParams(globalSpaceIdx);
    }

    /*!
     * \brief Returns the ECL material law manager
     *
     * Note that this method is *not* part of the generic eWoms problem API because it
     * would force all problens use the ECL material laws.
     */
    std::shared_ptr<const EclMaterialLawManager> materialLawManager() const
    { return materialLawManager_; }

    /*!
     * \copydoc materialLawManager()
     */
    std::shared_ptr<EclMaterialLawManager> materialLawManager()
    { return materialLawManager_; }

    /*!
     * \brief Returns the initial solvent saturation for a given a cell index
     */
    Scalar solventSaturation(unsigned elemIdx) const
    {
        if (solventSaturation_.empty())
            return 0;

        return solventSaturation_[elemIdx];
    }

    /*!
     * \brief Returns the initial polymer concentration for a given a cell index
     */
    Scalar  polymerConcentration(unsigned elemIdx) const
    {
        if (polymerConcentration_.empty())
            return 0;

        return polymerConcentration_[elemIdx];
    }

    /*!
    * \brief Returns the polymer molecule weight for a given cell index
    */
    // TODO: remove this function if not called
    Scalar polymerMolecularWeight(const unsigned elemIdx) const
    {
        if (polymerMoleWeight_.empty())
            return 0.0;

        return polymerMoleWeight_[elemIdx];
    }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned pvtRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return pvtRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the index the relevant PVT region given a cell index
     */
    unsigned pvtRegionIndex(unsigned elemIdx) const
    {
        if (pvtnum_.empty())
            return 0;

        return pvtnum_[elemIdx];
    }

    const std::vector<int>& pvtRegionArray() const
    { return pvtnum_; }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned satnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return satnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the index the relevant saturation function region given a cell index
     */
    unsigned satnumRegionIndex(unsigned elemIdx) const
    {
        if (satnum_.empty())
            return 0;

        return satnum_[elemIdx];
    }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned miscnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return miscnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the index the relevant MISC region given a cell index
     */
    unsigned miscnumRegionIndex(unsigned elemIdx) const
    {
        if (miscnum_.empty())
            return 0;

        return miscnum_[elemIdx];
    }

    /*!
     * \brief Returns the index of the relevant region for thermodynmic properties
     */
    template <class Context>
    unsigned plmixnumRegionIndex(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return plmixnumRegionIndex(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the index the relevant PLMIXNUM (for polymer module) region given a cell index
     */
    unsigned plmixnumRegionIndex(unsigned elemIdx) const
    {
        if (plmixnum_.empty())
            return 0;

        return plmixnum_[elemIdx];
    }

    /*!
     * \brief Returns the max polymer adsorption value
     */
    template <class Context>
    Scalar maxPolymerAdsorption(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    { return maxPolymerAdsorption(context.globalSpaceIndex(spaceIdx, timeIdx)); }

    /*!
     * \brief Returns the max polymer adsorption value
     */
    Scalar maxPolymerAdsorption(unsigned elemIdx) const
    {
        if (maxPolymerAdsorption_.empty())
            return 0;

        return maxPolymerAdsorption_[elemIdx];
    }

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    { return this->simulator().vanguard().caseName(); }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        // use the initial temperature of the DOF if temperature is not a primary
        // variable
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        return initialFluidStates_[globalDofIdx].temperature(/*phaseIdx=*/0);
    }

    /*!
     * \copydoc FvBaseProblem::boundary
     *
     * ECLIPSE uses no-flow conditions for all boundaries. \todo really?
     */
    template <class Context>
    void boundary(BoundaryRateVector& values,
                  const Context& context,
                  unsigned spaceIdx,
                  unsigned timeIdx) const
    {
        if (!enableEnergy || !enableThermalFluxBoundaries)
            values.setNoFlow();
        else {
            // in the energy case we need to specify a non-trivial boundary condition
            // because the geothermal gradient needs to be maintained. for this, we
            // simply assume the initial temperature at the boundary and specify the
            // thermal flow accordingly. in this context, "thermal flow" means energy
            // flow due to a temerature gradient while assuming no-flow for mass
            unsigned interiorDofIdx = context.interiorScvIndex(spaceIdx, timeIdx);
            unsigned globalDofIdx = context.globalSpaceIndex(interiorDofIdx, timeIdx);
            values.setThermalFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
        }

        if (nonTrivialBoundaryConditions()) {
            unsigned indexInInside = context.intersection(spaceIdx).indexInInside();
            unsigned interiorDofIdx = context.interiorScvIndex(spaceIdx, timeIdx);
            unsigned globalDofIdx = context.globalSpaceIndex(interiorDofIdx, timeIdx);
            unsigned pvtRegionIdx = pvtRegionIndex(context, spaceIdx, timeIdx);
            switch (indexInInside) {
            case 0:
                if (freebcXMinus_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcXMinus_[globalDofIdx], pvtRegionIdx);
                break;
            case 1:
                if (freebcX_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcX_[globalDofIdx], pvtRegionIdx);
                break;
            case 2:
                if (freebcYMinus_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcYMinus_[globalDofIdx], pvtRegionIdx);
                break;
            case 3:
                if (freebcY_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcY_[globalDofIdx], pvtRegionIdx);
                break;
            case 4:
                if (freebcZMinus_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcZMinus_[globalDofIdx], pvtRegionIdx);
                break;
            case 5:
                if (freebcZ_[globalDofIdx])
                    values.setFreeFlow(context, spaceIdx, timeIdx, initialFluidStates_[globalDofIdx]);
                else
                    values.setMassRate(massratebcZ_[globalDofIdx], pvtRegionIdx);
                break;
            default:
                throw std::logic_error("invalid face index for boundary condition");

            }
        }
    }

    /*!
     * \copydoc FvBaseProblem::initial
     *
     * The reservoir problem uses a constant boundary condition for
     * the whole domain.
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        values.setPvtRegionIndex(pvtRegionIndex(context, spaceIdx, timeIdx));
        values.assignNaive(initialFluidStates_[globalDofIdx]);

        if (enableSolvent)
            values[Indices::solventSaturationIdx] = solventSaturation_[globalDofIdx];

        if (enablePolymer)
            values[Indices::polymerConcentrationIdx] = polymerConcentration_[globalDofIdx];

        if (enablePolymerMolarWeight)
            values[Indices::polymerMoleWeightIdx]= polymerMoleWeight_[globalDofIdx];

        if (enableBrine)
            values[Indices::saltConcentrationIdx] = initialFluidStates_[globalDofIdx].saltConcentration();

        values.checkDefined();
    }

    /*!
     * \copydoc FvBaseProblem::initialSolutionApplied()
     */
    void initialSolutionApplied()
    {
        // initialize the wells. Note that this needs to be done after initializing the
        // intrinsic permeabilities and the after applying the initial solution because
        // the well model uses these...
        wellModel_.init();

        // let the object for threshold pressures initialize itself. this is done only at
        // this point, because determining the threshold pressures may require to access
        // the initial solution.
        thresholdPressures_.finishInit();

        updateCompositionChangeLimits_();

        if (enableAquifers_)
            aquiferModel_.initialSolutionApplied();
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0 everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context,
                unsigned spaceIdx,
                unsigned timeIdx) const
    {
        rate = 0.0;

        wellModel_.computeTotalRatesForDof(rate, context, spaceIdx, timeIdx);

        // convert the source term from the total mass rate of the
        // cell to the one per unit of volume as used by the model.
        const unsigned globalDofIdx = context.globalSpaceIndex(spaceIdx, timeIdx);
        for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx) {
            rate[eqIdx] /= this->model().dofTotalVolume(globalDofIdx);

            Ewoms::Valgrind::CheckDefined(rate[eqIdx]);
            assert(Ewoms::isfinite(rate[eqIdx]));
        }

        if (enableAquifers_)
            aquiferModel_.addToSource(rate, context, spaceIdx, timeIdx);

        // if requested, compensate systematic mass loss for cells which were "well
        // behaved" in the last time step
        if (enableDriftCompensation_) {
            const auto& intQuants = context.intensiveQuantities(spaceIdx, timeIdx);
            const auto& simulator = this->simulator();
            const auto& model = this->model();

            // we need a higher maxCompensation than the Newton tolerance because the
            // current time step might be shorter than the last one
            Scalar maxCompensation = 10.0*model.newtonMethod().tolerance();

            Scalar poro = intQuants.referencePorosity();
            Scalar dt = simulator.timeStepSize();

            EqVector dofDriftRate = drift_[globalDofIdx];
            dofDriftRate /= dt*context.dofTotalVolume(spaceIdx, timeIdx);

            // compute the weighted total drift rate
            Scalar totalDriftRate = 0.0;
            for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx)
                totalDriftRate +=
                    std::abs(dofDriftRate[eqIdx])*dt*model.eqWeight(globalDofIdx, eqIdx)/poro;

            // make sure that we do not exceed the maximum rate of drift compensation
            if (totalDriftRate > maxCompensation)
                dofDriftRate *= maxCompensation/totalDriftRate;

            for (unsigned eqIdx = 0; eqIdx < numEq; ++ eqIdx)
                rate[eqIdx] -= dofDriftRate[eqIdx];
        }
    }

    /*!
     * \brief Returns the maximum value of the gas dissolution factor at the current time
     *        for a given degree of freedom.
     */
    Scalar maxGasDissolutionFactor(unsigned timeIdx, unsigned globalDofIdx) const
    {
        int pvtRegionIdx = pvtRegionIndex(globalDofIdx);
        if (!drsdtActive_() || maxDRs_[pvtRegionIdx] < 0.0)
            return std::numeric_limits<Scalar>::max()/2.0;

        // this is a bit hacky because it assumes that a time discretization with only
        // two time indices is used.
        if (timeIdx == 0)
            return lastRs_[globalDofIdx] + maxDRs_[pvtRegionIdx];
        else
            return lastRs_[globalDofIdx];
    }

    /*!
     * \brief Returns the maximum value of the oil vaporization factor at the current
     *        time for a given degree of freedom.
     */
    Scalar maxOilVaporizationFactor(unsigned timeIdx, unsigned globalDofIdx) const
    {
        int pvtRegionIdx = pvtRegionIndex(globalDofIdx);
        if (!drvdtActive_() || maxDRv_[pvtRegionIdx] < 0.0)
            return std::numeric_limits<Scalar>::max()/2.0;

        // this is a bit hacky because it assumes that a time discretization with only
        // two time indices is used.
        if (timeIdx == 0)
            return lastRv_[globalDofIdx] + maxDRv_[pvtRegionIdx];
        else
            return lastRv_[globalDofIdx];
    }

    /*!
     * \brief Returns an element's historic maximum oil phase saturation that was
     *        observed during the simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This is a bit of a hack from the conceptional point of view, but it is required to
     * match the results of the 'eflow' and ECLIPSE 100 simulators.
     */
    Scalar maxOilSaturation(unsigned globalDofIdx) const
    {
        if (!vapparsActive())
            return 0.0;

        return maxOilSaturation_[globalDofIdx];
    }

    /*!
     * \brief Sets an element's maximum oil phase saturation observed during the
     *        simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This a hack on top of the maxOilSaturation() hack but it is currently required to
     * do restart externally, i.e. from the eflow code.
     */
    void setMaxOilSaturation(unsigned globalDofIdx, Scalar value)
    {
        if (!vapparsActive())
            return;

        maxOilSaturation_[globalDofIdx] = value;
    }

    /*!
     * \brief Returns an element's historic maximum water phase saturation that was
     *        observed during the simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This is used for output of the maximum water saturation used as input
     * for water induced rock compation ROCK2D/ROCK2DTR.
     */
    Scalar maxWaterSaturation(unsigned globalDofIdx) const
    {
        if (maxWaterSaturation_.empty())
            return 0.0;

        return maxWaterSaturation_[globalDofIdx];
    }

    /*!
     * \brief Returns an element's historic minimum pressure of the oil phase that was
     *        observed during the simulation.
     *
     * In this context, "historic" means the the time before the current timestep began.
     *
     * This is used for output of the minimum pressure used as input
     * for the irreversible rock compation option.
     */
    Scalar minOilPressure(unsigned globalDofIdx) const
    {
        if (minOilPressure_.empty())
            return 0.0;

        return minOilPressure_[globalDofIdx];
    }

    /*!
     * \brief Returns a reference to the ECL well manager used by the problem.
     *
     * This can be used for inspecting wells outside of the problem.
     */
    const EclWellModel& wellModel() const
    { return wellModel_; }

    EclWellModel& wellModel()
    { return wellModel_; }

    EclAquiferModel& mutableAquiferModel()
    { return aquiferModel_; }

    // temporary solution to facilitate output of initial state from eflow
    const InitialFluidState& initialFluidState(unsigned globalDofIdx) const
    { return initialFluidStates_[globalDofIdx]; }

    const Ewoms::EclipseIO& eclIO() const
    { return eclWriter_->eclIO(); }

    bool vapparsActive() const
    {
        const auto& simulator = this->simulator();
        int epsiodeIdx = std::max(simulator.episodeIndex(), 0);
        const auto& oilVaporizationControl = simulator.vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);
        return (oilVaporizationControl.getType() == Ewoms::OilVaporizationProperties::OilVaporization::VAPPARS);
    }

    bool nonTrivialBoundaryConditions() const
    { return nonTrivialBoundaryConditions_; }

    /*!
     * \brief Propose the size of the next time step to the simulator.
     *
     * This method is only called if the Newton solver does converge, the simulator
     * automatically cuts the time step in half without consultating this method again.
     */
    Scalar nextTimeStepSize() const
    {
        // allow external code to do the timestepping
        if (this->nextTimeStepSize_ > 0.0)
            return this->nextTimeStepSize_;

        const auto& simulator = this->simulator();
        int episodeIdx = simulator.episodeIndex();

        // for the initial episode, we use a fixed time step size
        if (episodeIdx < 0)
            return initialTimeStepSize_;

        // ask the newton method for a suggestion. This suggestion will be based on how
        // well the previous time step converged. After that, apply the runtime time
        // stepping constraints.
        const auto& newtonMethod = this->model().newtonMethod();
        return limitNextTimeStepSize_(newtonMethod.suggestTimeStepSize(simulator.timeStepSize()));
    }

    /*!
     * \brief Returns the minimum allowable size of a time step.
     */
    Scalar minTimeStepSize() const
    { return minTimeStepSize_; }

    /*!
     * \brief Returns the maximum number of subsequent failures for the time integration
     *        before giving up.
     */
    unsigned maxTimeIntegrationFailures() const
    { return maxFails_; }

    /*!
     * \brief Calculate the porosity multiplier due to water induced rock compaction.
     *
     * TODO: The API of this is a bit ad-hoc, it would be better to use context objects.
     */
    template <class LhsEval>
    LhsEval rockCompPoroMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {

        if (rockCompPoroMult_.empty() && rockCompPoroMultWc_.empty())
            return 1.0;

        unsigned tableIdx = 0;
        if (!rockTableIdx_.empty())
            tableIdx = rockTableIdx_[elementIdx];

        const auto& fs = intQuants.fluidState();
        LhsEval effectiveOilPressure = Ewoms::decay<LhsEval>(fs.pressure(oilPhaseIdx));
        if (!minOilPressure_.empty())
            // The pore space change is irreversible
            effectiveOilPressure =
                Ewoms::min(Ewoms::decay<LhsEval>(fs.pressure(oilPhaseIdx)),
                         minOilPressure_[elementIdx]);

        if (!overburdenPressure_.empty())
            effectiveOilPressure -= overburdenPressure_[elementIdx];

        if (!rockCompPoroMult_.empty()) {
            return rockCompPoroMult_[tableIdx].eval(effectiveOilPressure, /*extrapolation=*/true);
        }

        // water compaction
        assert(!rockCompPoroMultWc_.empty());
        LhsEval SwMax = Ewoms::max(Ewoms::decay<LhsEval>(fs.saturation(waterPhaseIdx)), maxWaterSaturation_[elementIdx]);
        LhsEval SwDeltaMax = SwMax - initialFluidStates_[elementIdx].saturation(waterPhaseIdx);

        return rockCompPoroMultWc_[tableIdx].eval(effectiveOilPressure, SwDeltaMax, /*extrapolation=*/true);
    }

    /*!
     * \brief Calculate the transmissibility multiplier due to water induced rock compaction.
     *
     * TODO: The API of this is a bit ad-hoc, it would be better to use context objects.
     */
    template <class LhsEval>
    LhsEval rockCompTransMultiplier(const IntensiveQuantities& intQuants, unsigned elementIdx) const
    {
        if (rockCompTransMult_.empty() && rockCompTransMultWc_.empty())
            return 1.0;

        unsigned tableIdx = 0;
        if (!rockTableIdx_.empty())
            tableIdx = rockTableIdx_[elementIdx];

        const auto& fs = intQuants.fluidState();
        LhsEval effectiveOilPressure = Ewoms::decay<LhsEval>(fs.pressure(oilPhaseIdx));

        if (!minOilPressure_.empty())
            // The pore space change is irreversible
            effectiveOilPressure =
                Ewoms::min(Ewoms::decay<LhsEval>(fs.pressure(oilPhaseIdx)),
                         minOilPressure_[elementIdx]);

        if (overburdenPressure_.size() > 0)
            effectiveOilPressure -= overburdenPressure_[elementIdx];

        if (!rockCompTransMult_.empty())
            return rockCompTransMult_[tableIdx].eval(effectiveOilPressure, /*extrapolation=*/true);

        // water compaction
        assert(!rockCompTransMultWc_.empty());
        LhsEval SwMax = Ewoms::max(Ewoms::decay<LhsEval>(fs.saturation(waterPhaseIdx)), maxWaterSaturation_[elementIdx]);
        LhsEval SwDeltaMax = SwMax - initialFluidStates_[elementIdx].saturation(waterPhaseIdx);

        return rockCompTransMultWc_[tableIdx].eval(effectiveOilPressure, SwDeltaMax, /*extrapolation=*/true);
    }

    /*!
     * \brief Get the pressure of the overburden.
     *
     * This method is mainly for output.
     */
    Scalar overburdenPressure(unsigned elementIdx) const
    {
        if (overburdenPressure_.empty())
            return 0.0;

        return overburdenPressure_[elementIdx];
    }

private:
    void checkDeckCompatibility_() const
    {
        const auto& deck = this->simulator().vanguard().deck();
        const auto& comm = this->simulator().gridView().comm();
        bool beVerbose = comm.rank() == 0;

        if (enableApiTracking)
            throw std::logic_error("API tracking is not yet implemented but requested at compile time.");
        if (!enableApiTracking && deck.hasKeyword("API"))
            throw std::logic_error("The simulator is build with API tracking disabled, but API tracking is requested by the deck.");

        if (enableSolvent && !deck.hasKeyword("SOLVENT"))
            throw std::runtime_error("The simulator requires the solvent option to be enabled, but the deck does not.");
        else if (!enableSolvent && deck.hasKeyword("SOLVENT"))
            throw std::runtime_error("The deck enables the solvent option, but the simulator is compiled without it.");

        if (enablePolymer && !deck.hasKeyword("POLYMER"))
            throw std::runtime_error("The simulator requires the polymer option to be enabled, but the deck does not.");
        else if (!enablePolymer && deck.hasKeyword("POLYMER"))
            throw std::runtime_error("The deck enables the polymer option, but the simulator is compiled without it.");

        if (deck.hasKeyword("TEMP") && deck.hasKeyword("THERMAL"))
            throw std::runtime_error("The deck enables both, the TEMP and the THERMAL options, but they are mutually exclusive.");

        bool deckEnergyEnabled = (deck.hasKeyword("TEMP") || deck.hasKeyword("THERMAL"));
        if (enableEnergy && !deckEnergyEnabled)
            throw std::runtime_error("The simulator requires the TEMP or the THERMAL option to be enabled, but the deck activates neither.");
        else if (!enableEnergy && deckEnergyEnabled)
            throw std::runtime_error("The deck enables the TEMP or the THERMAL option, but the simulator is not compiled to support either.");

        if (deckEnergyEnabled && deck.hasKeyword("TEMP") && beVerbose)
            std::cerr << "WARNING: The deck requests the TEMP option, i.e., treating energy "
                      << "conservation as a post processing step. This is currently unsupported, "
                      << "i.e., energy conservation is always handled fully implicitly." << std::endl;

        int numDeckPhases = FluidSystem::numActivePhases();
        if (numDeckPhases < Indices::numPhases && beVerbose)
            std::cerr << "WARNING: The number of active phases specified by the deck ("
                      << numDeckPhases << ") is smaller than the number of compiled-in phases ("
                      << Indices::numPhases << "). This usually results in a significant "
                      << "performance degradation compared to using a specialized simulator."  << std::endl;
        else if (numDeckPhases < Indices::numPhases)
            throw std::runtime_error("The deck enables "+std::to_string(numDeckPhases)+" phases "
                                     "while this simulator can only handle "+
                                     std::to_string(Indices::numPhases)+".");

        // make sure that the correct phases are active
        if (FluidSystem::phaseIsActive(oilPhaseIdx) && !Indices::oilEnabled)
            throw std::runtime_error("The deck enables oil, but this simulator cannot handle it.");
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && !Indices::gasEnabled)
            throw std::runtime_error("The deck enables gas, but this simulator cannot handle it.");
        if (FluidSystem::phaseIsActive(waterPhaseIdx) && !Indices::waterEnabled)
            throw std::runtime_error("The deck enables water, but this simulator cannot handle it.");
        // the opposite cases should be fine (albeit a bit slower than what's possible)
    }

    bool drsdtActive_() const
    {
        const auto& simulator = this->simulator();
        int epsiodeIdx = std::max(simulator.episodeIndex(), 0);
        const auto& oilVaporizationControl = simulator.vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);
        return (oilVaporizationControl.drsdtActive());
    }

    bool drvdtActive_() const
    {
        const auto& simulator = this->simulator();
        int epsiodeIdx = std::max(simulator.episodeIndex(), 0);
        const auto& oilVaporizationControl = simulator.vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);
        return (oilVaporizationControl.drvdtActive());

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

    void updateElementDepths_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& gridView = vanguard.gridView();
        const auto& elemMapper = this->elementMapper();;

        int numElements = gridView.size(/*codim=*/0);
        elementCenterDepth_.resize(numElements);

        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& element = *elemIt;
            const unsigned int elemIdx = elemMapper.index(element);

            elementCenterDepth_[elemIdx] = cellCenterDepth(element);
        }
    }

    // update the parameters needed for DRSDT and DRVDT
    void updateCompositionChangeLimits_()
    {
        // update the "last Rs" values for all elements, including the ones in the ghost
        // and overlap regions
        const auto& simulator = this->simulator();
        int epsiodeIdx = std::max(simulator.episodeIndex(), 0);
        const auto& oilVaporizationControl = simulator.vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);

        if (oilVaporizationControl.drsdtActive()) {
            ElementContext elemCtx(simulator);
            const auto& vanguard = simulator.vanguard();
            auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
            const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = iq.fluidState();

                typedef typename std::decay<decltype(fs)>::type FluidState;

                int pvtRegionIdx = pvtRegionIndex(compressedDofIdx);
                if (oilVaporizationControl.getOption(pvtRegionIdx) || fs.saturation(gasPhaseIdx) > freeGasMinSaturation_)
                    lastRs_[compressedDofIdx] =
                        Ewoms::BlackOil::template getRs_<FluidSystem,
                                                       FluidState,
                                                       Scalar>(fs, iq.pvtRegionIndex());
                else
                    lastRs_[compressedDofIdx] = std::numeric_limits<Scalar>::infinity();
            }
        }

        // update the "last Rv" values for all elements, including the ones in the ghost
        // and overlap regions
        if (drvdtActive_()) {
            ElementContext elemCtx(simulator);
            const auto& vanguard = simulator.vanguard();
            auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
            const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = iq.fluidState();

                typedef typename std::decay<decltype(fs)>::type FluidState;

                lastRv_[compressedDofIdx] =
                    Ewoms::BlackOil::template getRv_<FluidSystem,
                                                   FluidState,
                                                   Scalar>(fs, iq.pvtRegionIndex());
            }
        }
    }

    bool updateMaxOilSaturation_()
    {
        const auto& simulator = this->simulator();

        // we use VAPPARS
        if (vapparsActive()) {
            ElementContext elemCtx(simulator);
            const auto& vanguard = simulator.vanguard();
            auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
            const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const Element& elem = *elemIt;

                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = iq.fluidState();

                Scalar So = Ewoms::decay<Scalar>(fs.saturation(oilPhaseIdx));

                maxOilSaturation_[compressedDofIdx] = std::max(maxOilSaturation_[compressedDofIdx], So);
            }

            // we need to invalidate the intensive quantities cache here because the
            // derivatives of Rs and Rv will most likely have changed
            return true;
        }

        return false;
    }

    bool updateMaxWaterSaturation_()
    {
        // water compaction is activated in ROCKCOMP
        if (maxWaterSaturation_.size()== 0)
            return false;

        maxWaterSaturation_[/*timeIdx=*/1] = maxWaterSaturation_[/*timeIdx=*/0];
        ElementContext elemCtx(this->simulator());
        const auto& vanguard = this->simulator().vanguard();
        auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq.fluidState();

            Scalar Sw = Ewoms::decay<Scalar>(fs.saturation(waterPhaseIdx));
            maxWaterSaturation_[compressedDofIdx] = std::max(maxWaterSaturation_[compressedDofIdx], Sw);
        }

        return true;
    }

    bool updateMinPressure_()
    {
        // IRREVERS option is used in ROCKCOMP
        if (minOilPressure_.empty())
            return false;

        ElementContext elemCtx(this->simulator());
        const auto& vanguard = this->simulator().vanguard();
        auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& iq = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq.fluidState();

            minOilPressure_[compressedDofIdx] =
                std::min(minOilPressure_[compressedDofIdx],
                         Ewoms::getValue(fs.pressure(oilPhaseIdx)));
        }

        return true;
    }

    void readRockParameters_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& rockConfig = eclState.getSimulationConfig().rock_config();

        // read the rock compressibility parameters
        {
            const auto& comp = rockConfig.comp();
            rockParams_.clear();
            for (const auto& c : comp)
                rockParams_.push_back( { c.pref, c.compressibility } );
        }

        // read the parameters for water-induced rock compaction
        readRockCompactionParameters_();

        unsigned numElem = vanguard.gridView().size(0);
        if (eclState.fieldProps().has_int(rockConfig.rocknum_property())) {
            rockTableIdx_.resize(numElem);
            const auto& num = eclState.fieldProps().get_int(rockConfig.rocknum_property());
            for (size_t elemIdx = 0; elemIdx < numElem; ++ elemIdx) {
                rockTableIdx_[elemIdx] = num[elemIdx] - 1;
            }
        }

        // Store overburden pressure pr element
        const auto& overburdTables = eclState.getTableManager().getOverburdTables();
        if (!overburdTables.empty()) {
            overburdenPressure_.resize(numElem,0.0);
            size_t numRocktabTables = rockConfig.num_rock_tables();

            if (overburdTables.size() != numRocktabTables)
                throw std::runtime_error(std::to_string(numRocktabTables) +" OVERBURD tables is expected, but " + std::to_string(overburdTables.size()) +" is provided");

            std::vector<Ewoms::Tabulated1DFunction<Scalar>> overburdenTables(numRocktabTables);
            for (size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
                const Ewoms::OverburdTable& overburdTable =  overburdTables.template getTable<Ewoms::OverburdTable>(regionIdx);
                overburdenTables[regionIdx].setXYContainers(overburdTable.getDepthColumn(),overburdTable.getOverburdenPressureColumn());
            }

            for (size_t elemIdx = 0; elemIdx < numElem; ++ elemIdx) {
                unsigned tableIdx = 0;
                if (!rockTableIdx_.empty()) {
                    tableIdx = rockTableIdx_[elemIdx];
                }
                overburdenPressure_[elemIdx] = overburdenTables[tableIdx].eval(elementCenterDepth_[elemIdx], /*extrapolation=*/true);
            }
        }

    }

    void readRockCompactionParameters_()
    {
        const auto& vanguard = this->simulator().vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& rockConfig = eclState.getSimulationConfig().rock_config();

        if (!rockConfig.active())
            return; // deck does not enable rock compaction

        unsigned numElem = vanguard.gridView().size(0);
        switch (rockConfig.hysteresis_mode()) {
        case RockConfig::Hysteresis::REVERS:
            break;
        case RockConfig::Hysteresis::IRREVERS:
            // interpolate the porv volume multiplier using the minimum pressure in the cell
            // i.e. don't allow re-inflation.
            minOilPressure_.resize(numElem, 1e99);
            break;
        default:
            throw std::runtime_error("Not support ROCKOMP hysteresis option ");
        }

        size_t numRocktabTables = rockConfig.num_rock_tables();
        bool waterCompaction = rockConfig.water_compaction();

        if (!waterCompaction) {
            const auto& rocktabTables = eclState.getTableManager().getRocktabTables();
            if (rocktabTables.size() != numRocktabTables)
                throw std::runtime_error("ROCKCOMP is activated." + std::to_string(numRocktabTables)
                                         +" ROCKTAB tables is expected, but " + std::to_string(rocktabTables.size()) +" is provided");

            rockCompPoroMult_.resize(numRocktabTables);
            rockCompTransMult_.resize(numRocktabTables);
            for (size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
                const auto& rocktabTable = rocktabTables.template getTable<Ewoms::RocktabTable>(regionIdx);
                const auto& pressureColumn = rocktabTable.getPressureColumn();
                const auto& poroColumn = rocktabTable.getPoreVolumeMultiplierColumn();
                const auto& transColumn = rocktabTable.getTransmissibilityMultiplierColumn();
                rockCompPoroMult_[regionIdx].setXYContainers(pressureColumn, poroColumn);
                rockCompTransMult_[regionIdx].setXYContainers(pressureColumn, transColumn);
            }
        } else {
            const auto& rock2dTables = eclState.getTableManager().getRock2dTables();
            const auto& rock2dtrTables = eclState.getTableManager().getRock2dtrTables();
            const auto& rockwnodTables = eclState.getTableManager().getRockwnodTables();
            maxWaterSaturation_.resize(numElem, 0.0);

            if (rock2dTables.size() != numRocktabTables)
                throw std::runtime_error("Water compation option is selected in ROCKCOMP." + std::to_string(numRocktabTables)
                                         +" ROCK2D tables is expected, but " + std::to_string(rock2dTables.size()) +" is provided");

            if (rockwnodTables.size() != numRocktabTables)
                throw std::runtime_error("Water compation option is selected in ROCKCOMP." + std::to_string(numRocktabTables)
                                         +" ROCKWNOD tables is expected, but " + std::to_string(rockwnodTables.size()) +" is provided");
            //TODO check size match
            rockCompPoroMultWc_.resize(numRocktabTables, TabulatedTwoDFunction(TabulatedTwoDFunction::InterpolationPolicy::Vertical));
            for (size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
                const Ewoms::RockwnodTable& rockwnodTable =  rockwnodTables.template getTable<Ewoms::RockwnodTable>(regionIdx);
                const auto& rock2dTable = rock2dTables[regionIdx];

                if (rockwnodTable.getSaturationColumn().size() != rock2dTable.sizeMultValues())
                    throw std::runtime_error("Number of entries in ROCKWNOD and ROCK2D needs to match.");

                for (size_t xIdx = 0; xIdx < rock2dTable.size(); ++xIdx) {
                    rockCompPoroMultWc_[regionIdx].appendXPos(rock2dTable.getPressureValue(xIdx));
                    for (size_t yIdx = 0; yIdx < rockwnodTable.getSaturationColumn().size(); ++yIdx)
                        rockCompPoroMultWc_[regionIdx].appendSamplePoint(xIdx,
                                                                           rockwnodTable.getSaturationColumn()[yIdx],
                                                                           rock2dTable.getPvmultValue(xIdx, yIdx));
                }
            }

            if (!rock2dtrTables.empty()) {
                rockCompTransMultWc_.resize(numRocktabTables, TabulatedTwoDFunction(TabulatedTwoDFunction::InterpolationPolicy::Vertical));
                for (size_t regionIdx = 0; regionIdx < numRocktabTables; ++regionIdx) {
                    const Ewoms::RockwnodTable& rockwnodTable =  rockwnodTables.template getTable<Ewoms::RockwnodTable>(regionIdx);
                    const auto& rock2dtrTable = rock2dtrTables[regionIdx];

                    if (rockwnodTable.getSaturationColumn().size() != rock2dtrTable.sizeMultValues())
                        throw std::runtime_error("Number of entries in ROCKWNOD and ROCK2DTR needs to match.");

                    for (size_t xIdx = 0; xIdx < rock2dtrTable.size(); ++xIdx) {
                        rockCompTransMultWc_[regionIdx].appendXPos(rock2dtrTable.getPressureValue(xIdx));
                        for (size_t yIdx = 0; yIdx < rockwnodTable.getSaturationColumn().size(); ++yIdx)
                            rockCompTransMultWc_[regionIdx].appendSamplePoint(xIdx,
                                                                                     rockwnodTable.getSaturationColumn()[yIdx],
                                                                                     rock2dtrTable.getTransMultValue(xIdx, yIdx));
                    }
                }
            }
        }
    }

    void readMaterialParameters_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        // the PVT and saturation region numbers
        updatePvtnum_();
        updateSatnum_();

        // the MISC region numbers (solvent model)
        updateMiscnum_();
        // the PLMIX region numbers (polymer model)
        updatePlmixnum_();

        ////////////////////////////////
        // porosity
        updateReferencePorosity_();
        referencePorosity_[1] = referencePorosity_[0];
        ////////////////////////////////

        ////////////////////////////////
        // fluid-matrix interactions (saturation functions; relperm/capillary pressure)
        materialLawManager_ = std::make_shared<EclMaterialLawManager>();
        materialLawManager_->initFromEclState(eclState);
        materialLawManager_->initParamsForElements(eclState, this->model().numGridDof());
        ////////////////////////////////
    }

    void readThermalParameters_()
    {
        if (!enableEnergy)
            return;

        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        // fluid-matrix interactions (saturation functions; relperm/capillary pressure)
        thermalLawManager_ = std::make_shared<EclThermalLawManager>();
        thermalLawManager_->initParamsForElements(eclState, this->model().numGridDof());
    }

    void updateReferencePorosity_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        size_t numDof = this->model().numGridDof();

        referencePorosity_[/*timeIdx=*/0].resize(numDof);

        const auto& fieldProps = eclState.fieldProps();
        const std::vector<double> porvData = fieldProps.porv(true);
        const std::vector<int> actnumData = fieldProps.actnum();
        for (size_t dofIdx = 0; dofIdx < numDof; ++ dofIdx) {
            unsigned cartElemIdx = vanguard.cartesianIndex(dofIdx);
            Scalar poreVolume = porvData[cartElemIdx];

            // we define the porosity as the accumulated pore volume divided by the
            // geometric volume of the element. Note that -- in pathetic cases -- it can
            // be larger than 1.0!
            Scalar dofVolume = simulator.model().dofTotalVolume(dofIdx);
            assert(dofVolume > 0.0);
            referencePorosity_[/*timeIdx=*/0][dofIdx] = poreVolume/dofVolume;
        }
    }

    void initFluidSystem_()
    {
        const auto& simulator = this->simulator();
        const auto& eclState = simulator.vanguard().eclState();
        const auto& schedule = simulator.vanguard().schedule();

        FluidSystem::initFromEclState(eclState, schedule);
   }

    void readInitialCondition_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();

        if (eclState.getInitConfig().hasEquil())
            readEquilInitialCondition_();
        else
            readExplicitInitialCondition_();

        readBlackoilExtentionsInitialConditions_();

        //initialize min/max values
        size_t numElems = this->model().numGridDof();
        for (size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            const auto& fs = initialFluidStates_[elemIdx];
            if(maxWaterSaturation_.size() > 0)
                maxWaterSaturation_[elemIdx] = std::max(maxWaterSaturation_[elemIdx], fs.saturation(waterPhaseIdx));
            if(maxOilSaturation_.size() > 0)
                maxOilSaturation_[elemIdx] = std::max(maxOilSaturation_[elemIdx], fs.saturation(oilPhaseIdx));
            if(minOilPressure_.size() > 0)
                minOilPressure_[elemIdx] = std::min(minOilPressure_[elemIdx], fs.pressure(oilPhaseIdx));
        }

    }

    void readEquilInitialCondition_()
    {
        const auto& simulator = this->simulator();

        // initial condition corresponds to hydrostatic conditions.
        typedef Ewoms::EclEquilInitializer<TypeTag> EquilInitializer;
        EquilInitializer equilInitializer(simulator, *materialLawManager_);

        size_t numElems = this->model().numGridDof();
        initialFluidStates_.resize(numElems);
        for (size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = initialFluidStates_[elemIdx];
            elemFluidState.assign(equilInitializer.initialFluidState(elemIdx));
        }
    }

    void readEclRestartSolution_()
    {
        // Set the start time of the simulation
        auto& simulator = this->simulator();
        const auto& schedule = simulator.vanguard().schedule();
        const auto& eclState = simulator.vanguard().eclState();
        const auto& timeMap = schedule.getTimeMap();
        const auto& initconfig = eclState.getInitConfig();
        int episodeIdx = initconfig.getRestartStep();

        simulator.setStartTime(timeMap.getStartTime(/*timeStepIdx=*/0));
        simulator.setTime(timeMap.getTimePassedUntil(episodeIdx));

        simulator.startNextEpisode(simulator.startTime() + simulator.time(),
                                   timeMap.getTimeStepLength(episodeIdx));
        simulator.setEpisodeIndex(episodeIdx);

        eclWriter_->beginRestart();

        Scalar dt = std::min(eclWriter_->restartTimeStepSize(), simulator.episodeLength());
        simulator.setTimeStepSize(dt);

        size_t numElems = this->model().numGridDof();
        initialFluidStates_.resize(numElems);
        if (enableSolvent)
            solventSaturation_.resize(numElems, 0.0);

        if (enablePolymer)
            polymerConcentration_.resize(numElems, 0.0);

        if (enablePolymerMolarWeight) {
            const std::string msg {"Support of the RESTART for polymer molecular weight "
                                   "is not implemented yet. The polymer weight value will be "
                                   "zero when RESTART begins"};
            Ewoms::OpmLog::warning("NO_POLYMW_RESTART", msg);
            polymerMoleWeight_.resize(numElems, 0.0);
        }

        for (size_t elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            auto& elemFluidState = initialFluidStates_[elemIdx];
            elemFluidState.setPvtRegionIndex(pvtRegionIndex(elemIdx));
            eclWriter_->eclOutputModule().initHysteresisParams(simulator, elemIdx);
            eclWriter_->eclOutputModule().assignToFluidState(elemFluidState, elemIdx);

            // Note: Function processRestartSaturations_() mutates the
            // 'ssol' argument--the value from the restart file--if solvent
            // is enabled.  Then, store the updated solvent saturation into
            // 'solventSaturation_'.  Otherwise, just pass a dummy value to
            // the function and discard the unchanged result.  Do not index
            // into 'solventSaturation_' unless solvent is enabled.
            {
                auto ssol = enableSolvent
                    ? eclWriter_->eclOutputModule().getSolventSaturation(elemIdx)
                    : Scalar(0);

                processRestartSaturations_(elemFluidState, ssol);

                if (enableSolvent)
                    solventSaturation_[elemIdx] = ssol;
            }

            lastRs_[elemIdx] = elemFluidState.Rs();
            lastRv_[elemIdx] = elemFluidState.Rv();

            if (enablePolymer)
                 polymerConcentration_[elemIdx] = eclWriter_->eclOutputModule().getPolymerConcentration(elemIdx);
            // if we need to restart for polymer molecular weight simulation, we need to add related here
        }

        const int epsiodeIdx = simulator.episodeIndex();
        const auto& oilVaporizationControl = simulator.vanguard().schedule().getOilVaporizationProperties(epsiodeIdx);
        if (drsdtActive_())
            // DRSDT is enabled
            for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRs_.size(); ++pvtRegionIdx)
                maxDRs_[pvtRegionIdx] = oilVaporizationControl.getMaxDRSDT(pvtRegionIdx)*simulator.timeStepSize();

        if (drvdtActive_())
            // DRVDT is enabled
            for (size_t pvtRegionIdx = 0; pvtRegionIdx < maxDRv_.size(); ++pvtRegionIdx)
                maxDRv_[pvtRegionIdx] = oilVaporizationControl.getMaxDRVDT(pvtRegionIdx)*simulator.timeStepSize();

        if (tracerModel().numTracers() > 0 && this->gridView().comm().rank() == 0)
            std::cout << "Warning: Restart is not implemented for the tracer model, it will initialize itself "
                      << "with the initial tracer concentration.\n"
                      << std::flush;

        // assign the restart solution to the current solution. note that we still need
        // to compute real initial solution after this because the initial fluid states
        // need to be correct for stuff like boundary conditions.
        auto& sol = this->model().solution(/*timeIdx=*/0);
        const auto& gridView = this->gridView();
        ElementContext elemCtx(simulator);
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            if (elem.partitionType() != Dune::InteriorEntity)
                continue;

            elemCtx.updatePrimaryStencil(elem);
            int elemIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            initial(sol[elemIdx], elemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);
        }

        // make sure that the ghost and overlap entities exhibit the correct
        // solution. alternatively, this could be done in the loop above by also
        // considering non-interior elements. Since the initial() method might not work
        // 100% correctly for such elements, let's play safe and explicitly synchronize
        // using message passing.
        this->model().syncOverlap();

        eclWriter_->endRestart();
    }

    void processRestartSaturations_(InitialFluidState& elemFluidState, Scalar& solventSaturation)
    {
        // each phase needs to be above certain value to be claimed to be existing
        // this is used to recover some RESTART running with the defaulted single-precision format
        const Scalar smallSaturationTolerance = 1.e-6;
        Scalar sumSaturation = 0.0;
        for (size_t phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                if (elemFluidState.saturation(phaseIdx) < smallSaturationTolerance)
                    elemFluidState.setSaturation(phaseIdx, 0.0);

                sumSaturation += elemFluidState.saturation(phaseIdx);
            }

        }
        if (enableSolvent) {
            if (solventSaturation < smallSaturationTolerance)
                solventSaturation = 0.0;

           sumSaturation += solventSaturation;
        }

        assert(sumSaturation > 0.0);

        for (size_t phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (FluidSystem::phaseIsActive(phaseIdx)) {
                const Scalar saturation = elemFluidState.saturation(phaseIdx) / sumSaturation;
                elemFluidState.setSaturation(phaseIdx, saturation);
            }
        }
        if (enableSolvent) {
            solventSaturation = solventSaturation / sumSaturation;
        }
    }

    void readExplicitInitialCondition_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();
        const auto& fieldProps = eclState.fieldProps();
        bool hasSwat = fieldProps.has_double("SWAT");
        bool hasSgas = fieldProps.has_double("SGAS");
        bool hasRs = fieldProps.has_double("RS");
        bool hasRv = fieldProps.has_double("RV");
        bool hasPressure = fieldProps.has_double("PRESSURE");

        // make sure all required quantities are enables
        if (FluidSystem::phaseIsActive(waterPhaseIdx) && !hasSwat)
            throw std::runtime_error("The ECL input file requires the presence of the SWAT keyword if "
                                     "the water phase is active");
        if (FluidSystem::phaseIsActive(gasPhaseIdx) && !hasSgas)
            throw std::runtime_error("The ECL input file requires the presence of the SGAS keyword if "
                                     "the gas phase is active");

        if (!hasPressure)
             throw std::runtime_error("The ECL input file requires the presence of the PRESSURE "
                                      "keyword if the model is initialized explicitly");
        if (FluidSystem::enableDissolvedGas() && !hasRs)
            throw std::runtime_error("The ECL input file requires the RS keyword to be present if"
                                     " dissolved gas is enabled");
        if (FluidSystem::enableVaporizedOil() && !hasRv)
            throw std::runtime_error("The ECL input file requires the RV keyword to be present if"
                                     " vaporized oil is enabled");

        size_t numDof = this->model().numGridDof();

        initialFluidStates_.resize(numDof);

        std::vector<double> waterSaturationData;
        std::vector<double> gasSaturationData;
        std::vector<double> pressureData;
        std::vector<double> rsData;
        std::vector<double> rvData;
        std::vector<double> tempiData;

        if (FluidSystem::phaseIsActive(waterPhaseIdx))
            waterSaturationData = fieldProps.get_double("SWAT");
        else
            waterSaturationData.resize(numDof);

        if (FluidSystem::phaseIsActive(gasPhaseIdx))
            gasSaturationData = fieldProps.get_double("SGAS");
        else
            gasSaturationData.resize(numDof);

        pressureData = fieldProps.get_double("PRESSURE");
        if (FluidSystem::enableDissolvedGas())
            rsData = fieldProps.get_double("RS");

        if (FluidSystem::enableVaporizedOil())
            rvData = fieldProps.get_double("RV");

        // initial reservoir temperature
        tempiData = fieldProps.get_double("TEMPI");

        // calculate the initial fluid states
        for (size_t dofIdx = 0; dofIdx < numDof; ++dofIdx) {
            auto& dofFluidState = initialFluidStates_[dofIdx];

            dofFluidState.setPvtRegionIndex(pvtRegionIndex(dofIdx));

            //////
            // set temperature
            //////
            Scalar temperatureLoc = tempiData[dofIdx];
            if (!std::isfinite(temperatureLoc) || temperatureLoc <= 0)
                temperatureLoc = FluidSystem::surfaceTemperature;
            dofFluidState.setTemperature(temperatureLoc);

            //////
            // set saturations
            //////
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))
                dofFluidState.setSaturation(FluidSystem::waterPhaseIdx,
                                            waterSaturationData[dofIdx]);
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
                dofFluidState.setSaturation(FluidSystem::gasPhaseIdx,
                                            gasSaturationData[dofIdx]);
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx))
                dofFluidState.setSaturation(FluidSystem::oilPhaseIdx,
                                            1.0
                                            - waterSaturationData[dofIdx]
                                            - gasSaturationData[dofIdx]);

            //////
            // set phase pressures
            //////
            Scalar oilPressure = pressureData[dofIdx];

            // this assumes that capillary pressures only depend on the phase saturations
            // and possibly on temperature. (this is always the case for ECL problems.)
            Dune::FieldVector<Scalar, numPhases> pc(0.0);
            const auto& matParams = materialLawParams(dofIdx);
            MaterialLaw::capillaryPressures(pc, matParams, dofFluidState);
            Ewoms::Valgrind::CheckDefined(oilPressure);
            Ewoms::Valgrind::CheckDefined(pc);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                dofFluidState.setPressure(phaseIdx, oilPressure + (pc[phaseIdx] - pc[oilPhaseIdx]));
            }

            if (FluidSystem::enableDissolvedGas())
                dofFluidState.setRs(rsData[dofIdx]);
            else if (Indices::gasEnabled && Indices::oilEnabled)
                dofFluidState.setRs(0.0);

            if (FluidSystem::enableVaporizedOil())
                dofFluidState.setRv(rvData[dofIdx]);
            else if (Indices::gasEnabled && Indices::oilEnabled)
                dofFluidState.setRv(0.0);

            //////
            // set invB_
            //////
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx))
                    continue;

                const auto& b = FluidSystem::inverseFormationVolumeFactor(dofFluidState, phaseIdx, pvtRegionIndex(dofIdx));
                dofFluidState.setInvB(phaseIdx, b);

                const auto& rho = FluidSystem::density(dofFluidState, phaseIdx, pvtRegionIndex(dofIdx));
                dofFluidState.setDensity(phaseIdx, rho);

            }
        }
    }

    void readBlackoilExtentionsInitialConditions_()
    {
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& eclState = vanguard.eclState();
        size_t numDof = this->model().numGridDof();

        if (enableSolvent) {
            if (eclState.fieldProps().has_double("SSOL"))
                solventSaturation_ = eclState.fieldProps().get_double("SSOL");
            else
                solventSaturation_.resize(numDof, 0.0);
        }

        if (enablePolymer) {
            if (eclState.fieldProps().has_double("SPOLY"))
                polymerConcentration_ = eclState.fieldProps().get_double("SPOLY");
            else
                polymerConcentration_.resize(numDof, 0.0);
        }

        if (enablePolymerMolarWeight) {
            if (eclState.fieldProps().has_double("SPOLYMW"))
                polymerMoleWeight_ = eclState.fieldProps().get_double("SPOLYMW");
            else
                polymerMoleWeight_.resize(numDof, 0.0);
        }
    }

    // update the hysteresis parameters of the material laws for the whole grid
    bool updateHysteresis_()
    {
        if (!materialLawManager_->enableHysteresis())
            return false;

        // we need to update the hysteresis data for _all_ elements (i.e., not just the
        // interior ones) to avoid desynchronization of the processes in the parallel case!
        const auto& simulator = this->simulator();
        ElementContext elemCtx(simulator);
        const auto& vanguard = simulator.vanguard();
        auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            materialLawManager_->updateHysteresis(intQuants.fluidState(), compressedDofIdx);
        }
        return true;
    }

    void updateMaxPolymerAdsorption_()
    {
        // we need to update the max polymer adsoption data for all elements
        const auto& simulator = this->simulator();
        ElementContext elemCtx(simulator);
        const auto& vanguard = simulator.vanguard();
        auto elemIt = vanguard.gridView().template begin</*codim=*/0>();
        const auto& elemEndIt = vanguard.gridView().template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

            unsigned compressedDofIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);

            maxPolymerAdsorption_[compressedDofIdx] = std::max(maxPolymerAdsorption_[compressedDofIdx] , Ewoms::scalarValue(intQuants.polymerAdsorption()));
        }
    }

    template<class T>
    void updateNum(const std::string& name, std::vector<T>& numbers)
    {
        const auto& simulator = this->simulator();
        const auto& eclState = simulator.vanguard().eclState();

        if (!eclState.fieldProps().has_int(name))
            return;

        const auto& numData = eclState.fieldProps().get_int(name);
        const auto& vanguard = simulator.vanguard();

        unsigned numElems = vanguard.gridView().size(/*codim=*/0);
        numbers.resize(numElems);
        for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx) {
            numbers[elemIdx] = static_cast<T>(numData[elemIdx]) - 1;
        }
    }

    void updatePvtnum_()
    {
        updateNum("PVTNUM", pvtnum_);
    }

    void updateSatnum_()
    {
        updateNum("SATNUM", satnum_);
    }

    void updateMiscnum_()
    {
        updateNum("MISCNUM", miscnum_);
    }

    void updatePlmixnum_()
    {
        updateNum("PLMIXNUM", plmixnum_);
    }

    struct PffDofData_
    {
        Ewoms::ConditionalStorage<enableEnergy, Scalar> thermalHalfTrans;
        Scalar transmissibility;
    };

    // update the prefetch friendly data object
    void updatePffDofData_()
    {
        const auto& distFn =
            [this](PffDofData_& dofData,
                   const Stencil& stencil,
                   unsigned localDofIdx)
            -> void
        {
            const auto& elementMapper = this->model().elementMapper();

            unsigned globalElemIdx = elementMapper.index(stencil.entity(localDofIdx));
            if (localDofIdx != 0) {
                unsigned globalCenterElemIdx = elementMapper.index(stencil.entity(/*dofIdx=*/0));
                dofData.transmissibility = transmissibilities_.transmissibility(globalCenterElemIdx, globalElemIdx);

                if (enableEnergy)
                    *dofData.thermalHalfTrans = transmissibilities_.thermalHalfTrans(globalCenterElemIdx, globalElemIdx);
            }
        };

        pffDofData_.update(distFn);
    }

    void readBoundaryConditions_()
    {
        nonTrivialBoundaryConditions_ = false;
        const auto& simulator = this->simulator();
        const auto& vanguard = simulator.vanguard();
        const auto& bcconfig = vanguard.eclState().getSimulationConfig().bcconfig();
        if (bcconfig.size() > 0) {
            nonTrivialBoundaryConditions_ = true;

            size_t numCartDof = vanguard.cartesianSize();
            unsigned numElems = vanguard.gridView().size(/*codim=*/0);
            std::vector<int> cartesianToCompressedElemIdx(numCartDof, -1);

            for (unsigned elemIdx = 0; elemIdx < numElems; ++elemIdx)
                cartesianToCompressedElemIdx[vanguard.cartesianIndex(elemIdx)] = elemIdx;

            massratebcXMinus_.resize(numElems, 0.0);
            massratebcX_.resize(numElems, 0.0);
            massratebcYMinus_.resize(numElems, 0.0);
            massratebcY_.resize(numElems, 0.0);
            massratebcZMinus_.resize(numElems, 0.0);
            massratebcZ_.resize(numElems, 0.0);
            freebcX_.resize(numElems, false);
            freebcXMinus_.resize(numElems, false);
            freebcY_.resize(numElems, false);
            freebcYMinus_.resize(numElems, false);
            freebcZ_.resize(numElems, false);
            freebcZMinus_.resize(numElems, false);

            for (const auto& bcface : bcconfig) {
                const auto& type = bcface.bctype;
                if (type == BCType::RATE) {
                    int compIdx = 0; // default initialize to avoid -Wmaybe-uninitialized warning

                    switch (bcface.component) {
                    case BCComponent::OIL:
                        compIdx = Indices::canonicalToActiveComponentIndex(oilCompIdx);
                        break;
                    case BCComponent::GAS:
                        compIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
                        break;
                    case BCComponent::WATER:
                        compIdx = Indices::canonicalToActiveComponentIndex(waterCompIdx);
                        break;
                    case BCComponent::SOLVENT:
                        if (!enableSolvent)
                            throw std::logic_error("solvent is disabled and you're trying to add solvent to BC");

                        compIdx = Indices::solventSaturationIdx;
                        break;
                    case BCComponent::POLYMER:
                        if (!enablePolymer)
                            throw std::logic_error("polymer is disabled and you're trying to add polymer to BC");

                        compIdx = Indices::polymerConcentrationIdx;
                        break;
                    case BCComponent::NONE:
                        throw std::logic_error("you need to specify the component when RATE type is set in BC");
                        break;
                    }

                    std::vector<RateVector>* data = nullptr;
                    switch (bcface.dir) {
                    case FaceDir::XMinus:
                        data = &massratebcXMinus_;
                        break;
                    case FaceDir::XPlus:
                        data = &massratebcX_;
                        break;
                    case FaceDir::YMinus:
                        data = &massratebcYMinus_;
                        break;
                    case FaceDir::YPlus:
                        data = &massratebcY_;
                        break;
                    case FaceDir::ZMinus:
                        data = &massratebcZMinus_;
                        break;
                    case FaceDir::ZPlus:
                        data = &massratebcZ_;
                        break;
                    }

                    const Evaluation rate = bcface.rate;
                    for (int i = bcface.i1; i <= bcface.i2; ++i) {
                        for (int j = bcface.j1; j <= bcface.j2; ++j) {
                            for (int k = bcface.k1; k <= bcface.k2; ++k) {
                                std::array<int, 3> tmp = {{i,j,k}};
                                auto elemIdx = cartesianToCompressedElemIdx[vanguard.cartesianIndex(tmp)];
                                if (elemIdx >= 0)
                                    (*data)[elemIdx][compIdx] = rate;
                            }
                        }
                    }
                } else if (type == BCType::FREE) {
                    std::vector<bool>* data = nullptr;
                    switch (bcface.dir) {
                    case FaceDir::XMinus:
                        data = &freebcXMinus_;
                        break;
                    case FaceDir::XPlus:
                        data = &freebcX_;
                        break;
                    case FaceDir::YMinus:
                        data = &freebcYMinus_;
                        break;
                    case FaceDir::YPlus:
                        data = &freebcY_;
                        break;
                    case FaceDir::ZMinus:
                        data = &freebcZMinus_;
                        break;
                    case FaceDir::ZPlus:
                        data = &freebcZ_;
                        break;
                    }

                    for (int i = bcface.i1; i <= bcface.i2; ++i) {
                        for (int j = bcface.j1; j <= bcface.j2; ++j) {
                            for (int k = bcface.k1; k <= bcface.k2; ++k) {
                                std::array<int, 3> tmp = {{i,j,k}};
                                auto elemIdx = cartesianToCompressedElemIdx[vanguard.cartesianIndex(tmp)];
                                if (elemIdx >= 0)
                                    (*data)[elemIdx] = true;
                            }
                        }
                    }

                    // TODO: either the real initial solution needs to be computed or read from the restart file
                    const auto& eclState = simulator.vanguard().eclState();
                    const auto& initconfig = eclState.getInitConfig();
                    if (initconfig.restartRequested()) {
                        throw std::logic_error("restart is not compatible with using free boundary conditions");
                    }
                } else {
                    throw std::logic_error("invalid type for BC. Use FREE or RATE");
                }
            }
        }
    }

    // this method applies the runtime constraints specified via the deck and/or command
    // line parameters for the size of the next time step.
    Scalar limitNextTimeStepSize_(Scalar dtNext) const
    {
        if (!enableExperiments)
            return dtNext;

        const auto& simulator = this->simulator();
        const auto& events = simulator.vanguard().schedule().getEvents();
        int episodeIdx = simulator.episodeIndex();

        // first thing in the morning, limit the time step size to the maximum size
        dtNext = std::min(dtNext, maxTimeStepSize_);

        Scalar remainingEpisodeTime =
            simulator.episodeStartTime() + simulator.episodeLength()
            - (simulator.startTime() + simulator.time());
        assert(remainingEpisodeTime >= 0.0);

        // if we would have a small amount of time left over in the current episode, make
        // two equal time steps instead of a big and a small one
        if (remainingEpisodeTime/2.0 < dtNext && dtNext < remainingEpisodeTime*(1.0 - 1e-5))
            // note: limiting to the maximum time step size here is probably not strictly
            // necessary, but it should not hurt and is more fool-proof
            dtNext = std::min(maxTimeStepSize_, remainingEpisodeTime/2.0);

        if (simulator.episodeStarts()) {
            // if a well event occured, respect the limit for the maximum time step after
            // that, too
            int reportStepIdx = std::max(episodeIdx, 0);
            bool wellEventOccured =
                    events.hasEvent(Ewoms::ScheduleEvents::NEW_WELL, reportStepIdx)
                    || events.hasEvent(Ewoms::ScheduleEvents::PRODUCTION_UPDATE, reportStepIdx)
                    || events.hasEvent(Ewoms::ScheduleEvents::INJECTION_UPDATE, reportStepIdx)
                    || events.hasEvent(Ewoms::ScheduleEvents::WELL_STATUS_CHANGE, reportStepIdx);
            if (episodeIdx >= 0 && wellEventOccured && maxTimeStepAfterWellEvent_ > 0)
                dtNext = std::min(dtNext, maxTimeStepAfterWellEvent_);
        }

        return dtNext;
    }

    static std::string briefDescription_;

    std::array<std::vector<Scalar>, 2> referencePorosity_;
    std::vector<Scalar> elementCenterDepth_;
    EclTransmissibility<TypeTag> transmissibilities_;

    std::shared_ptr<EclMaterialLawManager> materialLawManager_;
    std::shared_ptr<EclThermalLawManager> thermalLawManager_;

    EclThresholdPressure<TypeTag> thresholdPressures_;

    std::vector<int> pvtnum_;
    std::vector<unsigned short> satnum_;
    std::vector<unsigned short> miscnum_;
    std::vector<unsigned short> plmixnum_;

    std::vector<unsigned short> rockTableIdx_;
    std::vector<RockParams> rockParams_;

    std::vector<Scalar> maxPolymerAdsorption_;

    std::vector<InitialFluidState> initialFluidStates_;

    std::vector<Scalar> polymerConcentration_;
    // polymer molecular weight
    std::vector<Scalar> polymerMoleWeight_;
    std::vector<Scalar> solventSaturation_;

    std::vector<bool> dRsDtOnlyFreeGas_; // apply the DRSDT rate limit only to cells that exhibit free gas
    std::vector<Scalar> lastRs_;
    std::vector<Scalar> maxDRs_;
    std::vector<Scalar> lastRv_;
    std::vector<Scalar> maxDRv_;
    constexpr static Scalar freeGasMinSaturation_ = 1e-7;
    std::vector<Scalar> maxOilSaturation_;
    std::vector<Scalar> maxWaterSaturation_;
    std::vector<Scalar> overburdenPressure_;
    std::vector<Scalar> minOilPressure_;

    std::vector<TabulatedTwoDFunction> rockCompPoroMultWc_;
    std::vector<TabulatedTwoDFunction> rockCompTransMultWc_;
    std::vector<TabulatedFunction> rockCompPoroMult_;
    std::vector<TabulatedFunction> rockCompTransMult_;

    bool enableDriftCompensation_;
    GlobalEqVector drift_;

    EclWellModel wellModel_;
    bool enableAquifers_;
    EclAquiferModel aquiferModel_;

    bool enableEclOutput_;
    std::unique_ptr<EclWriterType> eclWriter_;

    PffGridVector<GridView, Stencil, PffDofData_, DofMapper> pffDofData_;
    TracerModel tracerModel_;

    bool nonTrivialBoundaryConditions_;
    std::vector<bool> freebcX_;
    std::vector<bool> freebcXMinus_;
    std::vector<bool> freebcY_;
    std::vector<bool> freebcYMinus_;
    std::vector<bool> freebcZ_;
    std::vector<bool> freebcZMinus_;

    std::vector<RateVector> massratebcX_;
    std::vector<RateVector> massratebcXMinus_;
    std::vector<RateVector> massratebcY_;
    std::vector<RateVector> massratebcYMinus_;
    std::vector<RateVector> massratebcZ_;
    std::vector<RateVector> massratebcZMinus_;

    // time stepping parameters
    bool enableTuning_;
    Scalar initialTimeStepSize_;
    Scalar maxTimeStepAfterWellEvent_;
    Scalar maxTimeStepSize_;
    Scalar restartShrinkFactor_;
    unsigned maxFails_;
    Scalar minTimeStepSize_;
};

template <class TypeTag>
std::string EclProblem<TypeTag>::briefDescription_;

} // namespace Ewoms

#endif
