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
#ifndef EWOMS_BLACKOILWELLMODEL_HH
#define EWOMS_BLACKOILWELLMODEL_HH

#include <ewoms/eclio/opmlog/opmlog.hh>

#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <ewoms/common/optional.hh>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <stddef.h>

#include <ewoms/eclio/parser/eclipsestate/runspec.hh>

#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wellteststate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/guiderate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/group.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/gconsale.hh>

#include <ewoms/eclsimulators/timestepping/simulatorreport.hh>
#include <ewoms/eclsimulators/wells/perforationdata.hh>
#include <ewoms/eclsimulators/wells/vfpinjproperties.hh>
#include <ewoms/eclsimulators/wells/vfpprodproperties.hh>
#include <ewoms/eclsimulators/eflow/countglobalcells.hh>
#include <ewoms/eclsimulators/wells/wellstatefullyimplicitblackoil.hh>
#include <ewoms/eclsimulators/wells/rateconverter.hh>
#include <ewoms/eclsimulators/wells/wellinterface.hh>
#include <ewoms/eclsimulators/wells/standardwell.hh>
#include <ewoms/eclsimulators/wells/multisegmentwell.hh>
#include <ewoms/eclsimulators/wells/wellgrouphelpers.hh>
#include <ewoms/eclsimulators/wells/wellprodindexcalculator.hh>
#include <ewoms/eclsimulators/wells/parallelwellinfo.hh>
#include <ewoms/eclsimulators/timestepping/gatherconvergencereport.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmatrix.hh>

#include <ewoms/common/densead/math.hh>

#include <ewoms/eclsimulators/utils/deferredlogger.hh>

BEGIN_PROPERTIES

NEW_PROP_TAG(EnableTerminalOutput);

END_PROPERTIES

namespace Ewoms {

        /// Class for handling the blackoil well model.
        template<typename TypeTag>
        class BlackoilWellModel : public Ewoms::BaseAuxiliaryModule<TypeTag>
        {
        public:
            // ---------      Types      ---------
            typedef WellStateFullyImplicitBlackoil WellState;
            typedef BlackoilModelParameters<TypeTag> ModelParameters;

            using Grid = GET_PROP_TYPE(TypeTag, Grid);
            using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
            using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
            using Indices = GET_PROP_TYPE(TypeTag, Indices);
            using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
            using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
            using RateVector = GET_PROP_TYPE(TypeTag, RateVector);
            using GlobalEqVector = GET_PROP_TYPE(TypeTag, GlobalEqVector);
            using SparseMatrixAdapter = GET_PROP_TYPE(TypeTag, SparseMatrixAdapter);

            typedef typename Ewoms::BaseAuxiliaryModule<TypeTag>::NeighborSet NeighborSet;

            static const int numEq = Indices::numEq;
            static const int solventSaturationIdx = Indices::solventSaturationIdx;

            // TODO: where we should put these types, WellInterface or Well Model?
            // or there is some other strategy, like TypeTag
            typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
            typedef Dune::BlockVector<VectorBlockType> BVector;

#if  DUNE_VERSION_NEWER_REV(DUNE_ISTL, 2 , 5, 1)
            // 3x3 matrix block inversion was unstable from at least 2.3 until and
            // including 2.5.0
            typedef Dune::FieldMatrix<Scalar, numEq, numEq > MatrixBlockType;
#else
            typedef Dune::FieldMatrix<Scalar, numEq, numEq > MatrixBlockType;
#endif

            typedef Ewoms::BlackOilPolymerModule<TypeTag> PolymerModule;

            // For the conversion between the surface volume rate and resrevoir voidage rate
            using RateConverterType = RateConverter::
                SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;

            BlackoilWellModel(Simulator& eebosSimulator);

            void init();

            /////////////
            // <eWoms auxiliary module stuff>
            /////////////
            unsigned numDofs() const
            // No extra dofs are inserted for wells. (we use a Schur complement.)
            { return 0; }

            void addNeighbors(std::vector<NeighborSet>& neighbors) const;

            void applyInitial()
            {}

            void linearize(SparseMatrixAdapter& jacobian, GlobalEqVector& res);

            void postSolve(GlobalEqVector& deltaX)
            {
                recoverWellSolutionAndUpdateWellState(deltaX);
            }

            /////////////
            // </ eWoms auxiliary module stuff>
            /////////////

            template <class Restarter>
            void deserialize(Restarter& /* res */)
            {
                // TODO (?)
            }

            /*!
             * \brief This method writes the complete state of the well
             *        to the harddisk.
             */
            template <class Restarter>
            void serialize(Restarter& /* res*/)
            {
                // TODO (?)
            }

            void beginEpisode()
            {
                beginReportStep(eebosSimulator_.episodeIndex());
            }

            void beginTimeStep();

            void beginIteration()
            {
                assemble(eebosSimulator_.model().newtonMethod().numIterations(),
                         eebosSimulator_.timeStepSize());
            }

            void endIteration()
            { }

            void endTimeStep()
            {
                timeStepSucceeded(eebosSimulator_.time(), eebosSimulator_.timeStepSize());
            }

            void endEpisode()
            {
                endReportStep();
            }

            template <class Context>
            void computeTotalRatesForDof(RateVector& rate,
                                         const Context& context,
                                         unsigned spaceIdx,
                                         unsigned timeIdx) const;

            using WellInterfacePtr = std::shared_ptr<WellInterface<TypeTag> >;
            WellInterfacePtr well(const std::string& wellName) const;

            void initFromRestartFile(const RestartValue& restartValues);

            Ewoms::data::GroupAndNetworkValues
            groupAndNetworkData(const int reportStepIdx, const Ewoms::Schedule& sched) const
            {
                auto grp_nwrk_values = ::Ewoms::data::GroupAndNetworkValues{};

                this->assignGroupValues(reportStepIdx, sched, grp_nwrk_values.groupData);
                this->assignNodeValues(grp_nwrk_values.nodeData);

                return grp_nwrk_values;
            }

            Ewoms::data::Wells wellData() const
            {
                auto wsrpt = well_state_.report(phase_usage_, Ewoms::UgGridHelpers::globalCell(grid()));

                this->assignWellGuideRates(wsrpt);
                this->assignShutConnections(wsrpt);

                return wsrpt;
            }

            // substract Binv(D)rw from r;
            void apply( BVector& r) const;

            // subtract B*inv(D)*C * x from A*x
            void apply(const BVector& x, BVector& Ax) const;

#if HAVE_CUDA || HAVE_OPENCL
            // accumulate the contributions of all Wells in the WellContributions object
            void getWellContributions(WellContributions& x) const;
#endif

            // apply well model with scaling of alpha
            void applyScaleAdd(const Scalar alpha, const BVector& x, BVector& Ax) const;

            // Check if well equations is converged.
            ConvergenceReport getWellConvergence(const std::vector<Scalar>& B_avg, const bool checkGroupConvergence = false) const;

            // return the internal well state, ignore the passed one.
            // Used by the legacy code to make it compatible with the legacy well models.
            const WellState& wellState(const WellState& well_state EWOMS_UNUSED) const;

            // return the internal well state
            const WellState& wellState() const;

            const SimulatorReportSingle& lastReport() const;

            void addWellContributions(SparseMatrixAdapter& jacobian) const
            {
                for ( const auto& well: well_container_ ) {
                    well->addWellContributions(jacobian);
                }
            }

            // called at the beginning of a report step
            void beginReportStep(const int time_step);

            /// Return true if any well has a THP constraint.
            bool hasTHPConstraints() const;

            /// Shut down any single well, but only if it is in prediction mode.
            /// Returns true if the well was actually found and shut.
            bool forceShutWellByNameIfPredictionMode(const std::string& wellname, const double simulation_time);

        protected:
            Simulator& eebosSimulator_;

            std::vector< Well > wells_ecl_;
            std::vector< std::vector<PerforationData> > well_perf_data_;
            std::vector< WellProdIndexCalculator > prod_index_calc_;

            std::vector< ParallelWellInfo > parallel_well_info_;
            std::vector< ParallelWellInfo* > local_parallel_well_info_;

            bool wells_active_;

            // a vector of all the wells.
            std::vector<WellInterfacePtr > well_container_;

            // map from logically cartesian cell indices to compressed ones
            std::vector<int> cartesian_to_compressed_;

            std::vector<bool> is_cell_perforated_;

            std::function<bool(const Well&)> is_shut_or_defunct_;

            void initializeWellProdIndCalculators();
            void initializeWellPerfData();
            void initializeWellState(const int           timeStepIdx,
                                     const int           globalNumWells,
                                     const SummaryState& summaryState);

            // create the well container
            std::vector<WellInterfacePtr > createWellContainer(const int time_step);

            WellInterfacePtr
            createWellPointer(const int wellID,
                              const int time_step) const;

            template <typename WellType>
            std::unique_ptr<WellType>
            createTypedWellPointer(const int wellID,
                                   const int time_step) const;

            WellInterfacePtr createWellForWellTest(const std::string& well_name, const int report_step, Ewoms::DeferredLogger& deferred_logger) const;

            WellState well_state_;
            WellState previous_well_state_;
            WellState well_state_nupcol_;

            const ModelParameters param_;
            bool terminal_output_;
            bool has_solvent_;
            bool has_zFraction_;
            bool has_polymer_;
            std::vector<int> pvt_region_idx_;
            PhaseUsage phase_usage_;
            size_t global_num_cells_;
            // the number of the cells in the local grid
            size_t local_num_cells_;
            double gravity_;
            std::vector<double> depth_;
            bool initial_step_;
            bool report_step_starts_;
            bool glift_debug = false;
            bool alternative_well_rate_init_;

            Ewoms::optional<int> last_run_wellpi_{};

            std::unique_ptr<RateConverterType> rateConverter_;
            std::unique_ptr<VFPProperties<VFPInjProperties,VFPProdProperties>> vfp_properties_;

            SimulatorReportSingle last_report_;

            WellTestState wellTestState_;
            std::unique_ptr<GuideRate> guideRate_;

            std::map<std::string, double> node_pressures_; // Storing network pressures for output.

            // used to better efficiency of calcuation
            mutable BVector scaleAddRes_;

            const Grid& grid() const
            { return eebosSimulator_.vanguard().grid(); }

            const EclipseState& eclState() const
            { return eebosSimulator_.vanguard().eclState(); }

            const Schedule& schedule() const
            { return eebosSimulator_.vanguard().schedule(); }

            void gliftDebug(
                const std::string &msg,
                Ewoms::DeferredLogger& deferred_logger) const;

            /// \brief Get the wells of our partition that are not shut.
            /// \param timeStepIdx The index of the time step.
            /// \param[out] globalNumWells the number of wells globally.
            std::vector< Well > getLocalNonshutWells(const int timeStepIdx,
                                                     int& globalNumWells) const;

            /// \brief Create the parallel well information
            /// \param localWells The local wells from ECL schedule
            std::vector< ParallelWellInfo* >
            createLocalParallelWellInfo(const std::vector<Well>& localWells);

            // compute the well fluxes and assemble them in to the reservoir equations as source terms
            // and in the well equations.
            void assemble(const int iterationIdx,
                          const double dt);

            // called at the end of a time step
            void timeStepSucceeded(const double& simulationTime, const double dt);

            // called at the end of a report step
            void endReportStep();

            // using the solution x to recover the solution xw for wells and applying
            // xw to update Well State
            void recoverWellSolutionAndUpdateWellState(const BVector& x);

            void updateWellControls(Ewoms::DeferredLogger& deferred_logger, const bool checkGroupControls);

            void updateAndCommunicateGroupData();
            void updateNetworkPressures();

            // setting the well_solutions_ based on well_state.
            void updatePrimaryVariables(Ewoms::DeferredLogger& deferred_logger);

            void setupCartesianToCompressed_(const int* global_cell, int local_num__cells);

            void computeRepRadiusPerfLength(const Grid& grid, Ewoms::DeferredLogger& deferred_logger);

            void computeAverageFormationFactor(std::vector<Scalar>& B_avg) const;

            // Calculating well potentials for each well
            void computeWellPotentials(std::vector<double>& well_potentials, const int reportStepIdx, Ewoms::DeferredLogger& deferred_logger);

            const std::vector<double>& wellPerfEfficiencyFactors() const;

            void calculateEfficiencyFactors(const int reportStepIdx);

            void calculateProductivityIndexValues(DeferredLogger& deferred_logger);

            // it should be able to go to prepareTimeStep(), however, the updateWellControls() and initPrimaryVariablesEvaluation()
            // makes it a little more difficult. unless we introduce if (iterationIdx != 0) to avoid doing the above functions
            // twice at the beginning of the time step
            /// Calculating the explict quantities used in the well calculation. By explicit, we mean they are cacluated
            /// at the beginning of the time step and no derivatives are included in these quantities
            void calculateExplicitQuantities(Ewoms::DeferredLogger& deferred_logger) const;

            SimulatorReportSingle solveWellEq(const std::vector<Scalar>& B_avg, const double dt, Ewoms::DeferredLogger& deferred_logger);

            void initPrimaryVariablesEvaluation() const;

            // The number of components in the model.
            int numComponents() const;

            int numLocalWells() const;

            int numPhases() const;

            void assembleWellEq(const std::vector<Scalar>& B_avg, const double dt, Ewoms::DeferredLogger& deferred_logger);

            // some preparation work, mostly related to group control and RESV,
            // at the beginning of each time step (Not report step)
            void prepareTimeStep(Ewoms::DeferredLogger& deferred_logger);

            void extractLegacyCellPvtRegionIndex_();

            void extractLegacyDepth_();

            /// return true if wells are available in the reservoir
            bool wellsActive() const;

            void setWellsActive(const bool wells_active);

            /// return true if wells are available on this process
            bool localWellsActive() const;

            /// upate the wellTestState related to economic limits
            void updateWellTestState(const double& simulationTime, WellTestState& wellTestState) const;

            void updatePerforationIntensiveQuantities();

            void wellTesting(const int timeStepIdx, const double simulationTime, Ewoms::DeferredLogger& deferred_logger);

            // convert well data from ewoms-eclio to well state from ewoms-eclsimulators
            void wellsToState( const data::Wells& wells,
                               const data::GroupAndNetworkValues& grpNwrkValues,
                               const PhaseUsage& phases,
                               const bool handle_ms_well,
                               WellStateFullyImplicitBlackoil& state ) const;

            // whether there exists any multisegment well open on this process
            bool anyMSWellOpenLocal() const;

            const Well& getWellEcl(const std::string& well_name) const;

            void updateGroupIndividualControls(Ewoms::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);
            void updateGroupIndividualControl(const Group& group, Ewoms::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);
            bool checkGroupConstraints(const Group& group, Ewoms::DeferredLogger& deferred_logger) const;
            Group::ProductionCMode checkGroupProductionConstraints(const Group& group, Ewoms::DeferredLogger& deferred_logger) const;
            Group::InjectionCMode checkGroupInjectionConstraints(const Group& group, const Phase& phase) const;
            void checkGconsaleLimits(const Group& group, WellState& well_state, Ewoms::DeferredLogger& deferred_logger ) const;

            void updateGroupHigherControls(Ewoms::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);
            void checkGroupHigherConstraints(const Group& group, Ewoms::DeferredLogger& deferred_logger, std::set<std::string>& switched_groups);

            void actionOnBrokenConstraints(const Group& group, const Group::ExceedAction& exceed_action, const Group::ProductionCMode& newControl, Ewoms::DeferredLogger& deferred_logger);

            void actionOnBrokenConstraints(const Group& group, const Group::InjectionCMode& newControl, const Phase& topUpPhase, Ewoms::DeferredLogger& deferred_logger);

            WellInterfacePtr getWell(const std::string& well_name) const;

            void updateWsolvent(const Group& group, const Schedule& schedule, const int reportStepIdx, const WellStateFullyImplicitBlackoil& wellState);

            void setWsolvent(const Group& group, const Schedule& schedule, const int reportStepIdx, double wsolvent);

            void runWellPIScaling(const int timeStepIdx, DeferredLogger& local_deferredLogger);

            void assignWellGuideRates(data::Wells& wsrpt) const;
            void assignShutConnections(data::Wells& wsrpt) const;
            void assignGroupValues(const int                               reportStepIdx,
                                   const Schedule&                         sched,
                                   std::map<std::string, data::GroupData>& gvalues) const;

            void assignNodeValues(std::map<std::string, data::NodeData>& gvalues) const;

            std::unordered_map<std::string, data::GroupGuideRates>
            calculateAllGroupGuiderates(const int reportStepIdx, const Schedule& sched) const;

            void assignGroupControl(const Group& group, data::GroupData& gdata) const;
            data::GuideRateValue getGuideRateValues(const Well& well) const;
            data::GuideRateValue getGuideRateValues(const Group& group) const;
            void getGuideRateValues(const GuideRate::RateVector& qs,
                                    const bool                   is_inj,
                                    const std::string&           wgname,
                                    data::GuideRateValue&        grval) const;

            void assignGroupGuideRates(const Group& group,
                                       const std::unordered_map<std::string, data::GroupGuideRates>& groupGuideRates,
                                       data::GroupData& gdata) const;
        };

} // namespace Ewoms

#include "blackoilwellmodel_impl.hh"
#endif
