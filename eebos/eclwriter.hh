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
 * \copydoc Ewoms::EclWriter
 */
#ifndef EWOMS_ECL_WRITER_HH
#define EWOMS_ECL_WRITER_HH

#include "collecttoiorank.hh"
#include "ecloutputblackoilmodule.hh"

#include <ewoms/numerics/models/blackoil/blackoilmodel.hh>
#include <ewoms/numerics/discretizations/ecfv/ecfvdiscretization.hh>
#include <ewoms/numerics/io/baseoutputwriter.hh>
#include <ewoms/common/parallel/tasklets.hh>

#include <eebos/nncsorter.hh>

#include <ewoms/eclio/output/eclipseio.hh>
#include <ewoms/eclio/output/restartvalue.hh>
#include <ewoms/eclio/output/summary.hh>
#include <ewoms/eclio/parser/units/unitsystem.hh>

#include <ewoms/eclsimulators/utils/parallelrestart.hh>
#include <ewoms/eclgrids/gridhelpers.hh>
#include <ewoms/eclgrids/utility/cartesiantocompressed.hh>

#include <ewoms/common/valgrind.hh>
#include <ewoms/common/exceptions.hh>

#include <ewoms/eclio/opmlog/opmlog.hh>

#include <list>
#include <utility>
#include <string>
#include <chrono>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

BEGIN_PROPERTIES

NEW_PROP_TAG(EnableEclOutput);
NEW_PROP_TAG(EnableAsyncEclOutput);
NEW_PROP_TAG(EclOutputDoublePrecision);
NEW_PROP_TAG(EquilGrid);

END_PROPERTIES

namespace Ewoms {

template <class TypeTag>
class EclWriter;

template <class TypeTag>
class EclOutputBlackOilModule;

template <class TypeTag>
class EclTransmissibility;

/*!
 * \brief Detect whether two cells are direct vertical neighbours.
 *
 * I.e. have the same i and j index and all cartesian cells between them
 * along the vertical column are inactive.
 *
 * \tparam CM The type of the cartesian index mapper.
 * \param cartMapper The mapper onto cartesian indices.
 * \param cartesianToActive The mapping of cartesian indices to active indices.
 * \param smallGlobalIndex The cartesian cell index of the cell with smaller index
 * \param largeGlobalIndex The cartesian cell index of the cell with larger index
 * \return True if the cells have the same i and j indices and all cartesian cells
 *         between them are inactive.
 */
inline
bool directVerticalNeighbors(const std::array<int, 3>& cartDims,
                             const std::unordered_map<int,int>& cartesianToActive,
                             int smallGlobalIndex, int largeGlobalIndex)
{
    assert(smallGlobalIndex <= largeGlobalIndex);
    std::array<int, 3> ijk1, ijk2;
    auto globalToIjk = [cartDims](int gc) {
                           std::array<int, 3> ijk;
                           ijk[0] = gc % cartDims[0];
                           gc /= cartDims[0];
                           ijk[1] = gc % cartDims[1];
                           ijk[2] = gc / cartDims[1];
                           return ijk;
                       };
    ijk1 = globalToIjk(smallGlobalIndex);
    ijk2 = globalToIjk(largeGlobalIndex);
    assert(ijk2[2]>=ijk1[2]);

    if ( ijk1[0] == ijk2[0] && ijk1[1] == ijk2[1] && (ijk2[2] - ijk1[2]) > 1)
    {
        assert((largeGlobalIndex-smallGlobalIndex)%(cartDims[0]*cartDims[1])==0);
        for ( int gi = smallGlobalIndex + cartDims[0] * cartDims[1]; gi < largeGlobalIndex;
              gi += cartDims[0] * cartDims[1] )
        {
            if ( cartesianToActive.find( gi ) != cartesianToActive.end() )
            {
                return false;
            }
        }
        return true;
    } else
        return false;
}

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Collects necessary output values and pass it to ewoms-eclio.
 *
 * Caveats:
 * - For this class to do do anything meaningful, you will have to
 *   have the eWoms module ewoms-eclio.
 * - The only DUNE grid which is currently supported is Dune::CpGrid
 *   from the eWoms module "ewoms-eclgrids". Using another grid won't
 *   fail at compile time but you will provoke a fatal exception as
 *   soon as you try to write an ECL output file.
 * - This class requires to use the black oil model with the element
 *   centered finite volume discretization.
 */
template <class TypeTag>
class EclWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Vanguard) Vanguard;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, EquilGrid) EquilGrid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;

    typedef CollectDataToIORank<Vanguard> CollectDataToIORankType;

    typedef std::vector<Scalar> ScalarBuffer;

    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };

public:
    static void registerParameters()
    {
        EclOutputBlackOilModule<TypeTag>::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, bool, EnableAsyncEclOutput,
                             "Write the ECL-formated results in a non-blocking way (i.e., using a separate thread).");
    }

    // The Simulator object should preferably have been const - the
    // only reason that is not the case is due to the SummaryState
    // object owned deep down by the vanguard.
    EclWriter(Simulator& simulator)
        : simulator_(simulator)
        , collectToIORank_(simulator_.vanguard())
        , eclOutputModule_(simulator, collectToIORank_)
    {
        if (collectToIORank_.isIORank()) {
            eclIO_.reset(new Ewoms::EclipseIO(simulator_.vanguard().eclState(),
                                            Ewoms::UgGridHelpers::createEclipseGrid(globalGrid(), simulator_.vanguard().eclState().getInputGrid()),
                                            simulator_.vanguard().schedule(),
                                            simulator_.vanguard().summaryConfig()));
        }

        // create output thread if enabled and rank is I/O rank
        // async output is enabled by default if pthread are enabled
        bool enableAsyncOutput = EWOMS_GET_PARAM(TypeTag, bool, EnableAsyncEclOutput);
        int numWorkerThreads = 0;
        if (enableAsyncOutput && collectToIORank_.isIORank())
            numWorkerThreads = 1;
        taskletRunner_.reset(new TaskletRunner(numWorkerThreads));
    }

    ~EclWriter()
    { }

    const Ewoms::EclipseIO& eclIO() const
    {
        assert(eclIO_);
        return *eclIO_;
    }

    const EquilGrid& globalGrid() const
    {
        return simulator_.vanguard().equilGrid();
    }

    void writeInit()
    {
        if (collectToIORank_.isIORank()) {
            std::map<std::string, std::vector<int> > integerVectors;
            if (collectToIORank_.isParallel())
                integerVectors.emplace("MPI_RANK", collectToIORank_.globalRanks());
            auto cartMap = Ewoms::cartesianToCompressed(globalGrid().size(0),
                                                      Ewoms::UgGridHelpers::globalCell(globalGrid()));
            eclIO_->writeInitial(computeTrans_(cartMap), integerVectors, exportNncStructure_(cartMap));
        }
    }

    /*!
     * \brief collect and pass data and pass it to eclIO writer
     */

    void evalSummaryState(bool isSubStep)
    {
        int reportStepNum = simulator_.episodeIndex() + 1;
        /*
          The summary data is not evaluated for timestep 0, that is
          implemented with a:

             if (time_step == 0)
                 return;

          check somewhere in the summary code. When the summary code was
          split in separate methods Summary::eval() and
          Summary::add_timestep() it was necessary to pull this test out
          here to ensure that the well and group related keywords in the
          restart file, like XWEL and XGRP were "correct" also in the
          initial report step.

          "Correct" in this context means unchanged behavior, might very
          well be more correct to actually remove this if test.
        */
        if (reportStepNum == 0)
            return;

        Scalar curTime = simulator_.time() + simulator_.timeStepSize();
        Scalar totalCpuTime =
            simulator_.executionTimer().realTimeElapsed() +
            simulator_.setupTimer().realTimeElapsed() +
            simulator_.vanguard().externalSetupTime();

        Ewoms::data::Wells localWellData = simulator_.problem().wellModel().wellData();

        const auto& gridView = simulator_.vanguard().gridView();
        int numElements = gridView.size(/*codim=*/0);
        bool log = collectToIORank_.isIORank();
        eclOutputModule_.allocBuffers(numElements, reportStepNum, isSubStep, log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);
        ElementIterator elemIt = gridView.template begin</*codim=*/0>();
        const ElementIterator& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            eclOutputModule_.processElement(elemCtx);
        }

        if (collectToIORank_.isParallel())
            collectToIORank_.collect({}, eclOutputModule_.getBlockData(), localWellData);

        std::map<std::string, double> miscSummaryData;
        std::map<std::string, std::vector<double>> regionData;
        eclOutputModule_.outputFipLog(miscSummaryData, regionData, isSubStep);

        bool forceDisableProdOutput = false;
        bool forceDisableInjOutput = false;
        bool forceDisableCumOutput = false;
        eclOutputModule_.outputProdLog(reportStepNum, isSubStep, forceDisableProdOutput);
        eclOutputModule_.outputInjLog(reportStepNum, isSubStep, forceDisableInjOutput);
        eclOutputModule_.outputCumLog(reportStepNum, isSubStep, forceDisableCumOutput);

        std::vector<char> buffer;
        if (collectToIORank_.isIORank()) {
            const auto& summary = eclIO_->summary();
            const auto& eclState = simulator_.vanguard().eclState();

            // Add TCPU
            if (totalCpuTime != 0.0)
                miscSummaryData["TCPU"] = totalCpuTime;

            const Ewoms::data::Wells& wellData = collectToIORank_.isParallel() ? collectToIORank_.globalWellData() : localWellData;

            const std::map<std::pair<std::string, int>, double>& blockData
                = collectToIORank_.isParallel()
                ? collectToIORank_.globalBlockData()
                : eclOutputModule_.getBlockData();

            summary.eval(summaryState(),
                         reportStepNum,
                         curTime,
                         eclState,
                         schedule(),
                         wellData,
                         miscSummaryData,
                         regionData,
                         blockData);
            buffer = summaryState().serialize();
        }

        if (collectToIORank_.isParallel()) {
#ifdef HAVE_MPI
            unsigned long buffer_size = buffer.size();
            MPI_Bcast(&buffer_size, 1, MPI_UNSIGNED_LONG, collectToIORank_.ioRank, MPI_COMM_WORLD);
            if (!collectToIORank_.isIORank())
                buffer.resize( buffer_size );

            MPI_Bcast(buffer.data(), buffer_size, MPI_CHAR, collectToIORank_.ioRank, MPI_COMM_WORLD);
            if (!collectToIORank_.isIORank()) {
                Ewoms::SummaryState& st = summaryState();
                st.deserialize(buffer);
            }
#endif
        }
    }

    void writeOutput(bool isSubStep)
    {
        Scalar curTime = simulator_.time() + simulator_.timeStepSize();
        Scalar nextStepSize = simulator_.problem().nextTimeStepSize();

        // output using eclWriter if enabled
        Ewoms::data::Wells localWellData = simulator_.problem().wellModel().wellData();

        int reportStepNum = simulator_.episodeIndex() + 1;
        const auto& gridView = simulator_.vanguard().gridView();
        int numElements = gridView.size(/*codim=*/0);
        bool log = collectToIORank_.isIORank();
        eclOutputModule_.allocBuffers(numElements, reportStepNum, isSubStep, log, /*isRestart*/ false);

        ElementContext elemCtx(simulator_);
        ElementIterator elemIt = gridView.template begin</*codim=*/0>();
        const ElementIterator& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const Element& elem = *elemIt;
            elemCtx.updatePrimaryStencil(elem);
            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            eclOutputModule_.processElement(elemCtx);
        }
        eclOutputModule_.outputErrorLog();

        // collect all data to I/O rank and assign to sol
        Ewoms::data::Solution localCellData = {};
        if (!isSubStep)
            eclOutputModule_.assignToSolution(localCellData);

        // add cell data to perforations for Rft output
        if (!isSubStep)
            eclOutputModule_.addRftDataToWells(localWellData, reportStepNum);

        if (collectToIORank_.isParallel())
            collectToIORank_.collect(localCellData, eclOutputModule_.getBlockData(), localWellData);

        if (collectToIORank_.isIORank()) {
            const auto& eclState = simulator_.vanguard().eclState();
            const auto& simConfig = eclState.getSimulationConfig();

            bool enableDoublePrecisionOutput = EWOMS_GET_PARAM(TypeTag, bool, EclOutputDoublePrecision);
            const Ewoms::data::Solution& cellData = collectToIORank_.isParallel() ? collectToIORank_.globalCellData() : localCellData;
            const Ewoms::data::Wells& wellData = collectToIORank_.isParallel() ? collectToIORank_.globalWellData() : localWellData;
            Ewoms::RestartValue restartValue(cellData, wellData);

            if (simConfig.useThresholdPressure())
                restartValue.addExtra("THRESHPR", Ewoms::UnitSystem::measure::pressure, simulator_.problem().thresholdPressure().data());

            // Add suggested next timestep to extra data.
            if (!isSubStep)
                restartValue.addExtra("OPMEXTRA", std::vector<double>(1, nextStepSize));

            // first, create a tasklet to write the data for the current time step to disk
            auto eclWriteTasklet = std::make_shared<EclWriteTasklet>(summaryState(),
                                                                     *eclIO_,
                                                                     reportStepNum,
                                                                     isSubStep,
                                                                     curTime,
                                                                     restartValue,
                                                                     enableDoublePrecisionOutput);

            // then, make sure that the previous I/O request has been completed and the
            // number of incomplete tasklets does not increase between time steps
            taskletRunner_->barrier();

            // finally, start a new output writing job
            taskletRunner_->dispatch(eclWriteTasklet);
        }
    }

    void beginRestart()
    {
        bool enableHysteresis = simulator_.problem().materialLawManager()->enableHysteresis();
        bool enableSwatinit = simulator_.vanguard().eclState().fieldProps().has_double("SWATINIT");
        std::vector<Ewoms::RestartKey> solutionKeys{
            {"PRESSURE", Ewoms::UnitSystem::measure::pressure},
            {"SWAT", Ewoms::UnitSystem::measure::identity, static_cast<bool>(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))},
            {"SGAS", Ewoms::UnitSystem::measure::identity, static_cast<bool>(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))},
            {"TEMP" , Ewoms::UnitSystem::measure::temperature, enableEnergy},
            {"SSOLVENT" , Ewoms::UnitSystem::measure::identity, enableSolvent},
            {"RS", Ewoms::UnitSystem::measure::gas_oil_ratio, FluidSystem::enableDissolvedGas()},
            {"RV", Ewoms::UnitSystem::measure::oil_gas_ratio, FluidSystem::enableVaporizedOil()},
            {"SOMAX", Ewoms::UnitSystem::measure::identity, simulator_.problem().vapparsActive()},
            {"PCSWM_OW", Ewoms::UnitSystem::measure::identity, enableHysteresis},
            {"KRNSW_OW", Ewoms::UnitSystem::measure::identity, enableHysteresis},
            {"PCSWM_GO", Ewoms::UnitSystem::measure::identity, enableHysteresis},
            {"KRNSW_GO", Ewoms::UnitSystem::measure::identity, enableHysteresis},
            {"PPCW", Ewoms::UnitSystem::measure::pressure, enableSwatinit}
        };

        const auto& inputThpres = eclState().getSimulationConfig().getThresholdPressure();
        std::vector<Ewoms::RestartKey> extraKeys = {{"OPMEXTRA", Ewoms::UnitSystem::measure::identity, false},
                                                  {"THRESHPR", Ewoms::UnitSystem::measure::pressure, inputThpres.active()}};

        // The episodeIndex is rewined one back before beginRestart is called
        // and can not be used here.
        // We just ask the initconfig directly to be sure that we use the correct
        // index.
        const auto& initconfig = simulator_.vanguard().eclState().getInitConfig();
        int restartStepIdx = initconfig.getRestartStep();

        const auto& gridView = simulator_.vanguard().gridView();
        unsigned numElements = gridView.size(/*codim=*/0);
        eclOutputModule_.allocBuffers(numElements, restartStepIdx, /*isSubStep=*/false, /*log=*/false, /*isRestart*/ true);

        {
            Ewoms::SummaryState& summaryState = simulator_.vanguard().summaryState();
            auto restartValues = loadParallelRestart(eclIO_.get(), summaryState, solutionKeys, extraKeys,
                                                     gridView.grid().comm());

            for (unsigned elemIdx = 0; elemIdx < numElements; ++elemIdx) {
                unsigned globalIdx = collectToIORank_.localIdxToGlobalIdx(elemIdx);
                eclOutputModule_.setRestart(restartValues.solution, elemIdx, globalIdx);
            }

            if (inputThpres.active()) {
                Simulator& mutableSimulator = const_cast<Simulator&>(simulator_);
                auto& thpres = mutableSimulator.problem().thresholdPressure();
                const auto& thpresValues = restartValues.getExtra("THRESHPR");
                thpres.setFromRestart(thpresValues);
            }
            restartTimeStepSize_ = restartValues.getExtra("OPMEXTRA")[0];

            // initialize the well model from restart values
            simulator_.problem().wellModel().initFromRestartFile(restartValues);

            if (!restartValues.aquifer.empty())
                simulator_.problem().mutableAquiferModel().initFromRestart(restartValues.aquifer);
        }
    }

    void endRestart()
    {}

    const EclOutputBlackOilModule<TypeTag>& eclOutputModule() const
    { return eclOutputModule_; }

    Scalar restartTimeStepSize() const
    { return restartTimeStepSize_; }

private:
    static bool enableEclOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, EnableEclOutput); }

    Ewoms::data::Solution computeTrans_(const std::unordered_map<int,int>& cartesianToActive) const
    {
        const auto& cartMapper = simulator_.vanguard().equilCartesianIndexMapper();
        const auto& cartDims = cartMapper.cartesianDimensions();
        const int globalSize = cartDims[0]*cartDims[1]*cartDims[2];

        Ewoms::data::CellData tranx = {Ewoms::UnitSystem::measure::transmissibility, std::vector<double>(globalSize), Ewoms::data::TargetType::INIT};
        Ewoms::data::CellData trany = {Ewoms::UnitSystem::measure::transmissibility, std::vector<double>(globalSize), Ewoms::data::TargetType::INIT};
        Ewoms::data::CellData tranz = {Ewoms::UnitSystem::measure::transmissibility, std::vector<double>(globalSize), Ewoms::data::TargetType::INIT};

        for (size_t i = 0; i < tranx.data.size(); ++i) {
            tranx.data[0] = 0.0;
            trany.data[0] = 0.0;
            tranz.data[0] = 0.0;
        }

        typedef typename EquilGrid :: LeafGridView  GlobalGridView;
        const GlobalGridView& globalGridView = globalGrid().leafGridView();
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView> ElementMapper;
        ElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());
#else
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView, Dune::MCMGElementLayout> ElementMapper;
        ElementMapper globalElemMapper(globalGridView);
#endif

        const EclTransmissibility<TypeTag>* globalTrans;

        if (!collectToIORank_.isParallel())
        {
            // in the sequential case we must use the transmissibilites defined by
            // the problem. (because in the sequential case, the grid manager does
            // not compute "global" transmissibilities for performance reasons. in
            // the parallel case, the problem's transmissibilities can't be used
            // because this object refers to the distributed grid and we need the
            // sequential version here.)
            globalTrans = &simulator_.problem().eclTransmissibilities();
        }
        else
        {
            globalTrans = &(simulator_.vanguard().globalTransmissibility());
        }

        auto elemIt = globalGridView.template begin</*codim=*/0>();
        const auto& elemEndIt = globalGridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            const auto& elem = *elemIt;

            auto isIt = globalGridView.ibegin(elem);
            const auto& isEndIt = globalGridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                const auto& is = *isIt;

                if (!is.neighbor())
                    continue; // intersection is on the domain boundary

                unsigned c1 = globalElemMapper.index(is.inside());
                unsigned c2 = globalElemMapper.index(is.outside());

                if (c1 > c2)
                    continue; // we only need to handle each connection once, thank you.

                // Ordering of compressed and uncompressed index should be the same
                const int cartIdx1 = cartMapper.cartesianIndex( c1 );
                const int cartIdx2 = cartMapper.cartesianIndex( c2 );
                // Ordering of compressed and uncompressed index should be the same
                assert(cartIdx1 <= cartIdx2);
                int gc1 = std::min(cartIdx1, cartIdx2);
                int gc2 = std::max(cartIdx1, cartIdx2);

                if (gc2 - gc1 == 1) {
                    tranx.data[gc1] = globalTrans->transmissibility(c1, c2);
                    continue; // skip other if clauses as they are false, last one needs some computation
                }

                if (gc2 - gc1 == cartDims[0]) {
                    trany.data[gc1] = globalTrans->transmissibility(c1, c2);
                    continue; // skipt next if clause as it needs some computation
                }

                if ( gc2 - gc1 == cartDims[0]*cartDims[1] ||
                     directVerticalNeighbors(cartDims, cartesianToActive, gc1, gc2))
                    tranz.data[gc1] = globalTrans->transmissibility(c1, c2);
            }
        }

        return {{"TRANX", tranx},
                {"TRANY", trany},
                {"TRANZ", tranz}};
    }

    Ewoms::NNC exportNncStructure_(const std::unordered_map<int,int>& cartesianToActive) const
    {
        std::size_t nx = eclState().getInputGrid().getNX();
        std::size_t ny = eclState().getInputGrid().getNY();
        auto nncData = sortNncAndApplyEditnnc(eclState().getInputNNC().data(),
                                              eclState().getInputEDITNNC().data());
        const auto& unitSystem = simulator_.vanguard().deck().getActiveUnitSystem();
        std::vector<Ewoms::NNCdata> outputNnc;
        std::size_t index = 0;

        for( const auto& entry : nncData ) {
            // test whether NNC is not a neighboring connection
            // cell2>=cell1 holds due to sortNncAndApplyEditnnc
            assert( entry.cell2 >= entry.cell1 );
            auto cellDiff = entry.cell2 - entry.cell1;

            if (cellDiff != 1 && cellDiff != nx && cellDiff != nx*ny) {
                auto tt = unitSystem.from_si(Ewoms::UnitSystem::measure::transmissibility, entry.trans);
                // Eclipse ignores NNCs (with EDITNNC applied) that are small. Seems like the threshold is 1.0e-6
                if ( tt >= 1.0e-6 )
                    outputNnc.emplace_back(entry.cell1, entry.cell2, entry.trans);
            }
            ++index;
        }

        auto nncCompare =  []( const Ewoms::NNCdata& nnc1, const Ewoms::NNCdata& nnc2){
                               return nnc1.cell1 < nnc2.cell1 ||
                                      ( nnc1.cell1 == nnc2.cell1 && nnc1.cell2 < nnc2.cell2);};
        // Sort the nncData values from the deck as they need to be
        // Checked when writing NNC transmissibilities from the simulation.
        std::sort(nncData.begin(), nncData.end(), nncCompare);

        typedef typename EquilGrid :: LeafGridView  GlobalGridView;
        const GlobalGridView& globalGridView = globalGrid().leafGridView();
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView> ElementMapper;
        ElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());

#else
        typedef Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView, Dune::MCMGElementLayout> ElementMapper;
        ElementMapper globalElemMapper(globalGridView);
#endif

        const EclTransmissibility<TypeTag>* globalTrans;
        if (!collectToIORank_.isParallel()) {
            // in the sequential case we must use the transmissibilites defined by
            // the problem. (because in the sequential case, the grid manager does
            // not compute "global" transmissibilities for performance reasons. in
            // the parallel case, the problem's transmissibilities can't be used
            // because this object refers to the distributed grid and we need the
            // sequential version here.)
            globalTrans = &simulator_.problem().eclTransmissibilities();
        }
        else
        {
            globalTrans = &(simulator_.vanguard().globalTransmissibility());
        }

        // Cartesian index mapper for the serial I/O grid
        const auto& equilCartMapper =  simulator_.vanguard().equilCartesianIndexMapper();

        const auto& cartDims = simulator_.vanguard().cartesianIndexMapper().cartesianDimensions();
        auto elemIt = globalGridView.template begin</*codim=*/0>();
        const auto& elemEndIt = globalGridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++ elemIt) {
            const auto& elem = *elemIt;

            auto isIt = globalGridView.ibegin(elem);
            const auto& isEndIt = globalGridView.iend(elem);
            for (; isIt != isEndIt; ++ isIt) {
                const auto& is = *isIt;

                if (!is.neighbor())
                    continue; // intersection is on the domain boundary

                unsigned c1 = globalElemMapper.index(is.inside());
                unsigned c2 = globalElemMapper.index(is.outside());

                if (c1 > c2)
                    continue; // we only need to handle each connection once, thank you.

                std::size_t cc1 = equilCartMapper.cartesianIndex( c1 ); //globalIOGrid_.globalCell()[c1];
                std::size_t cc2 = equilCartMapper.cartesianIndex( c2 ); //globalIOGrid_.globalCell()[c2];

                if ( cc2 < cc1 )
                    std::swap(cc1, cc2);

                auto cellDiff = cc2 - cc1;

                if (cellDiff != 1 &&
                    cellDiff != nx &&
                    cellDiff != nx*ny &&
                    ! directVerticalNeighbors(cartDims, cartesianToActive, cc1, cc2)) {
                    // We need to check whether an NNC for this face was also specified
                    // via the NNC keyword in the deck (i.e. in the first origNncSize entries.
                    auto t = globalTrans->transmissibility(c1, c2);
                    auto candidate = std::lower_bound(nncData.begin(), nncData.end(), Ewoms::NNCdata(cc1, cc2, 0.0), nncCompare);

                    while ( candidate != nncData.end() && candidate->cell1 == cc1
                         && candidate->cell2 == cc2) {
                        t -= candidate->trans;
                        ++candidate;
                    }
                    // eclipse ignores NNCs with zero transmissibility (different threshold than for NNC
                    // with corresponding EDITNNC above). In addition we do set small transmissibilties
                    // to zero when setting up the simulator. These will be ignored here, too.
                    auto tt = unitSystem.from_si(Ewoms::UnitSystem::measure::transmissibility, std::abs(t));
                    if ( tt > 1e-12 )
                        outputNnc.push_back({cc1, cc2, t});
                }
            }
        }
        Ewoms::NNC ret;
        for(const auto& nncItem: outputNnc)
            ret.addNNC(nncItem.cell1, nncItem.cell2, nncItem.trans);
        return ret;
    }

    struct EclWriteTasklet
        : public TaskletInterface
    {
        Ewoms::SummaryState summaryState_;
        Ewoms::EclipseIO& eclIO_;
        int reportStepNum_;
        bool isSubStep_;
        double secondsElapsed_;
        Ewoms::RestartValue restartValue_;
        bool writeDoublePrecision_;

        explicit EclWriteTasklet(const Ewoms::SummaryState& summaryState,
                                 Ewoms::EclipseIO& eclIO,
                                 int reportStepNum,
                                 bool isSubStep,
                                 double secondsElapsed,
                                 Ewoms::RestartValue restartValue,
                                 bool writeDoublePrecision)
            : summaryState_(summaryState)
            , eclIO_(eclIO)
            , reportStepNum_(reportStepNum)
            , isSubStep_(isSubStep)
            , secondsElapsed_(secondsElapsed)
            , restartValue_(restartValue)
            , writeDoublePrecision_(writeDoublePrecision)
        { }

        // callback to eclIO serial writeTimeStep method
        void run()
        {
            eclIO_.writeTimeStep(summaryState_,
                                 reportStepNum_,
                                 isSubStep_,
                                 secondsElapsed_,
                                 restartValue_,
                                 writeDoublePrecision_);
        }
    };

    const Ewoms::EclipseState& eclState() const
    { return simulator_.vanguard().eclState(); }

    Ewoms::SummaryState& summaryState()
    { return simulator_.vanguard().summaryState(); }

    const Ewoms::Schedule& schedule() const
    { return simulator_.vanguard().schedule(); }

    Simulator& simulator_;
    CollectDataToIORankType collectToIORank_;
    EclOutputBlackOilModule<TypeTag> eclOutputModule_;
    std::unique_ptr<Ewoms::EclipseIO> eclIO_;
    std::unique_ptr<TaskletRunner> taskletRunner_;
    Scalar restartTimeStepSize_;

};
} // namespace Ewoms

#endif
