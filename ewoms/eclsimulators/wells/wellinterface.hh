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
#ifndef EWOMS_WELLINTERFACE_HH
#define EWOMS_WELLINTERFACE_HH

#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/errormacros.hh>
#include <ewoms/eclio/exceptions.hh>

#include <ewoms/eclio/parser/eclipsestate/schedule/well/well.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wellteststate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/guiderate.hh>

#include <ewoms/eclsimulators/deprecated/props/blackoilphases.hh>
#include <ewoms/eclsimulators/timestepping/simulatorreport.hh>

#include <ewoms/eclsimulators/wells/rateconverter.hh>
#include <ewoms/eclsimulators/wells/vfpproperties.hh>
#include <ewoms/eclsimulators/wells/wellhelpers.hh>
#include <ewoms/eclsimulators/wells/wellgrouphelpers.hh>
#include <ewoms/eclsimulators/wells/wellprodindexcalculator.hh>
#include <ewoms/eclsimulators/wells/wellstatefullyimplicitblackoil.hh>
#include <ewoms/eclsimulators/eflow/blackoilmodelparameters.hh>

#include <ewoms/eclsimulators/timestepping/convergencereport.hh>
#include <ewoms/eclsimulators/utils/deferredlogger.hh>

#include<dune/common/fmatrix.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/matrixmatrix.hh>

#include <ewoms/common/densead/math.hh>
#include <ewoms/common/densead/evaluation.hh>

#include <string>
#include <memory>
#include <vector>
#include <cassert>

namespace Ewoms
{

    template<typename TypeTag>
    class WellInterface
    {
    public:

        using WellState = WellStateFullyImplicitBlackoil;

        typedef BlackoilModelParameters<TypeTag> ModelParameters;

        static const int Water = BlackoilPhases::Aqua;
        static const int Oil = BlackoilPhases::Liquid;
        static const int Gas = BlackoilPhases::Vapour;

        using Grid = GET_PROP_TYPE(TypeTag, Grid);
        using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
        using FluidSystem = GET_PROP_TYPE(TypeTag, FluidSystem);
        using Indices = GET_PROP_TYPE(TypeTag, Indices);
        using IntensiveQuantities = GET_PROP_TYPE(TypeTag, IntensiveQuantities);
        using MaterialLaw = GET_PROP_TYPE(TypeTag, MaterialLaw);
        using SparseMatrixAdapter = GET_PROP_TYPE(TypeTag, SparseMatrixAdapter);
        using RateVector = GET_PROP_TYPE(TypeTag, RateVector);

        static const int numEq = Indices::numEq;
        typedef double Scalar;

        typedef Dune::FieldVector<Scalar, numEq    > VectorBlockType;
        typedef Dune::FieldMatrix<Scalar, numEq, numEq > MatrixBlockType;
        typedef Dune::BlockVector<VectorBlockType> BVector;
        typedef DenseAd::Evaluation<double, /*size=*/numEq> Eval;

        static const bool has_solvent = GET_PROP_VALUE(TypeTag, EnableSolvent);
        static const bool has_zFraction = GET_PROP_VALUE(TypeTag, EnableSsaSolvent);
        static const bool has_polymer = GET_PROP_VALUE(TypeTag, EnablePolymer);
        static const bool has_energy = GET_PROP_VALUE(TypeTag, EnableEnergy);
        static const bool has_temperature = GET_PROP_VALUE(TypeTag, EnableTemperature);
        // flag for polymer molecular weight related
        static const bool has_polymermw = GET_PROP_VALUE(TypeTag, EnablePolymerMW);
        static const bool has_foam = GET_PROP_VALUE(TypeTag, EnableFoam);
        static const bool has_brine = GET_PROP_VALUE(TypeTag, EnableBrine);
        static const int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static const int contiZfracEqIdx = Indices::contiZfracEqIdx;
        static const int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
        // index for the polymer molecular weight continuity equation
        static const int contiPolymerMWEqIdx = Indices::contiPolymerMWEqIdx;
        static const int contiFoamEqIdx = Indices::contiFoamEqIdx;
        static const int contiBrineEqIdx = Indices::contiBrineEqIdx;

        // For the conversion between the surface volume rate and reservoir voidage rate
        using RateConverterType = RateConverter::
        SurfaceToReservoirVoidage<FluidSystem, std::vector<int> >;
        static const bool compositionSwitchEnabled = Indices::gasEnabled;
        using FluidState = Ewoms::BlackOilFluidState<Eval,
                                                   FluidSystem,
                                                   has_temperature,
                                                   has_energy,
                                                   compositionSwitchEnabled,
                                                   has_brine,
                                                   Indices::numPhases >;
        /// Constructor
        WellInterface(const Well& well,
                      const ParallelWellInfo& pw_info,
                      const int time_step,
                      const ModelParameters& param,
                      const RateConverterType& rate_converter,
                      const int pvtRegionIdx,
                      const int num_components,
                      const int num_phases,
                      const int index_of_well,
                      const int first_perf_index,
                      const std::vector<PerforationData>& perf_data);

        /// Virutal destructor
        virtual ~WellInterface() {}

        /// Well name.
        const std::string& name() const;

        /// True if the well is an injector.
        bool isInjector() const;

        /// True if the well is a producer.
        bool isProducer() const;

        /// Index of well in the wells struct and wellState
        int indexOfWell() const;

        /// Well cells.
        const std::vector<int>& cells() const {return well_cells_; }

        void setVFPProperties(const VFPProperties<VFPInjProperties,VFPProdProperties>* vfp_properties_arg);

        void setGuideRate(const GuideRate* guide_rate_arg);

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const int num_cells,
                          const std::vector< Scalar >& B_avg);

        virtual void initPrimaryVariablesEvaluation() const = 0;

        virtual ConvergenceReport getWellConvergence(const WellState& well_state, const std::vector<double>& B_avg, Ewoms::DeferredLogger& deferred_logger, const bool relax_tolerance = false) const = 0;

        virtual void solveEqAndUpdateWellState(WellState& well_state, Ewoms::DeferredLogger& deferred_logger) = 0;

        virtual void assembleWellEq(const Simulator& eebosSimulator,
                                    const std::vector<Scalar>& B_avg,
                                    const double dt,
                                    WellState& well_state,
                                    Ewoms::DeferredLogger& deferred_logger
                                    ) = 0;

        virtual void maybeDoGasLiftOptimization (
            WellState& well_state,
            const Simulator& eebosSimulator,
            DeferredLogger& deferred_logger
        ) const = 0;

        void updateWellTestState(const WellState& well_state,
                                 const double& simulationTime,
                                 const bool& writeMessageToOpmLog,
                                 WellTestState& wellTestState,
                                 Ewoms::DeferredLogger& deferred_logger) const;

        void setWellEfficiencyFactor(const double efficiency_factor);

        void computeRepRadiusPerfLength(const Grid& grid, const std::vector<int>& cartesian_to_compressed, Ewoms::DeferredLogger& deferred_logger);

        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        virtual void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                                           WellState& well_state,
                                                           Ewoms::DeferredLogger& deferred_logger) const = 0;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const = 0;

        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const = 0;

        // TODO: before we decide to put more information under mutable, this function is not const
        virtual void computeWellPotentials(const Simulator& eebosSimulator,
                                           const std::vector<Scalar>& B_avg,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials,
                                           Ewoms::DeferredLogger& deferred_logger) = 0;

        virtual void updateWellStateWithTarget(const Simulator& eebos_simulator,
                                               WellState& well_state,
                                               Ewoms::DeferredLogger& deferred_logger) const = 0;

        enum class IndividualOrGroup { Individual, Group, Both };
        bool updateWellControl(const Simulator& eebos_simulator,
                               const IndividualOrGroup iog,
                               WellState& well_state,
                               Ewoms::DeferredLogger& deferred_logger) /* const */;

        virtual void updatePrimaryVariables(const WellState& well_state, Ewoms::DeferredLogger& deferred_logger) const = 0;

        virtual void calculateExplicitQuantities(const Simulator& eebosSimulator,
                                                 const WellState& well_state,
                                                 Ewoms::DeferredLogger& deferred_logger) = 0; // should be const?

        virtual void updateProductivityIndex(const Simulator& eebosSimulator,
                                             const WellProdIndexCalculator& wellPICalc,
                                             WellState& well_state,
                                             DeferredLogger& deferred_logger) const = 0;

        /// \brief Wether the Jacobian will also have well contributions in it.
        virtual bool jacobianContainsWellContributions() const
        {
            return false;
        }

        // updating the voidage rates in well_state when requested
        void calculateReservoirRates(WellState& well_state) const;

        // Add well contributions to matrix
        virtual void addWellContributions(SparseMatrixAdapter&) const = 0;

        void addCellRates(RateVector& rates, int cellIdx) const;

        Scalar volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const;

        template <class EvalWell>
        Eval restrictEval(const EvalWell& in) const
        {
            Eval out = 0.0;
            out.setValue(in.value());
            for(int eqIdx = 0; eqIdx < numEq;++eqIdx) {
                out.setDerivative(eqIdx, in.derivative(eqIdx));
            }
            return out;
        }

        void closeCompletions(WellTestState& wellTestState);

        const Well& wellEcl() const;

        // TODO: theoretically, it should be a const function
        // Simulator is not const is because that assembleWellEq is non-const Simulator
        void wellTesting(const Simulator& simulator, const std::vector<double>& B_avg,
                         const double simulation_time, const int report_step,
                         const WellTestConfig::Reason testing_reason,
                         /* const */ WellState& well_state, WellTestState& welltest_state,
                         Ewoms::DeferredLogger& deferred_logger);

        void updatePerforatedCell(std::vector<bool>& is_cell_perforated);

        void checkWellOperability(const Simulator& eebos_simulator, const WellState& well_state, Ewoms::DeferredLogger& deferred_logger);

        // check whether the well is operable under the current reservoir condition
        // mostly related to BHP limit and THP limit
        void updateWellOperability(const Simulator& eebos_simulator,
                                   const WellState& well_state,
                                   Ewoms::DeferredLogger& deferred_logger);

        // whether the well is operable
        bool isOperable() const;

        /// Returns true if the well has one or more THP limits/constraints.
        bool wellHasTHPConstraints(const SummaryState& summaryState) const;

        /// Returns true if the well is currently in prediction mode (i.e. not history mode).
        bool underPredictionMode() const;

        // update perforation water throughput based on solved water rate
        virtual void updateWaterThroughput(const double dt, WellState& well_state) const = 0;

        /// Compute well rates based on current reservoir conditions and well variables.
        /// Used in updateWellStateRates().
        virtual std::vector<double> computeCurrentWellRates(const Simulator& eebosSimulator,
                                                            DeferredLogger& deferred_logger) const = 0;

        /// Modify the well_state's rates if there is only one nonzero rate.
        /// If so, that rate is kept as is, but the others are set proportionally
        /// to the rates returned by computeCurrentWellRates().
        void updateWellStateRates(const Simulator& eebosSimulator,
                                  WellState& well_state,
                                  DeferredLogger& deferred_logger) const;

        void stopWell() {
            wellIsStopped_ = true;
        }
        void openWell() {
            wellIsStopped_ = false;
        }

        bool wellIsStopped() const {
            return wellIsStopped_;
        }

        void setWsolvent(const double wsolvent);

        void setDynamicThpLimit(const double thp_limit);

    protected:

        // to indicate a invalid completion
        static const int INVALIDCOMPLETION = INT_MAX;

        Well well_ecl_;

        const ParallelWellInfo& parallel_well_info_;

        const int current_step_;

        // simulation parameters
        const ModelParameters& param_;

        // number of the perforations for this well
        int number_of_perforations_;

        // well index for each perforation
        std::vector<double> well_index_;

        // depth for each perforation
        std::vector<double> perf_depth_;

        // reference depth for the BHP
        double ref_depth_;

        double well_efficiency_factor_;

        // cell index for each well perforation
        std::vector<int> well_cells_;

        // saturation table nubmer for each well perforation
        std::vector<int> saturation_table_number_;

        // representative radius of the perforations, used in shear calculation
        std::vector<double> perf_rep_radius_;

        // length of the perforations, use in shear calculation
        std::vector<double> perf_length_;

        // well bore diameter
        std::vector<double> bore_diameters_;

        /*
         *  completions_ contains the mapping from completion id to connection indices
         *  {
         *      2 : [ConnectionIndex, ConnectionIndex],
         *      1 : [ConnectionIndex, ConnectionIndex, ConnectionIndex],
         *      5 : [ConnectionIndex],
         *      7 : [ConnectionIndex]
         *      ...
         *   }
         *   The integer IDs correspond to the COMPLETION id given by the COMPLUMP keyword.
         *   When there is no COMPLUMP keyword used, a default completion number will be assigned
         *   based on the order of the declaration of the connections.
         *   Since the connections not OPEN is not included in the Wells, so they will not be considered
         *   in this mapping relation.
         */
        std::map<int, std::vector<int>> completions_;

        const PhaseUsage* phase_usage_;

        bool getAllowCrossFlow() const;

        const VFPProperties<VFPInjProperties,VFPProdProperties>* vfp_properties_;

        const GuideRate* guide_rate_;

        double gravity_;

        // For the conversion between the surface volume rate and resrevoir voidage rate
        const RateConverterType& rateConverter_;

        // The pvt region of the well. We assume
        // We assume a well to not penetrate more than one pvt region.
        const int pvtRegionIdx_;

        const int num_components_;

        // number of phases
        int number_of_phases_;

        // the index of well in Wells struct
        int index_of_well_;

        // record the index of the first perforation
        // of states of individual well.
        int first_perf_;

        const std::vector<PerforationData>* perf_data_;

        std::vector<RateVector> connectionRates_;

        bool wellIsStopped_;

        double wsolvent_;

        Ewoms::optional<double> dynamic_thp_limit_;

        std::vector< Scalar > B_avg_;

        // the vectors used to describe the inflow performance relationship (IPR)
        // Q = IPR_A - BHP * IPR_B
        // TODO: it minght need to go to WellInterface, let us implement it in StandardWell first
        // it is only updated and used for producers for now
        mutable std::vector<double> ipr_a_;
        mutable std::vector<double> ipr_b_;

        bool changed_to_stopped_this_step_ = false;

        const PhaseUsage& phaseUsage() const;

        int eflowPhaseToEebosCompIdx( const int phaseIdx ) const;

        int eebosCompIdxToEFlowCompIdx( const unsigned compIdx ) const;

        double wsolvent() const;

        double wpolymer() const;

        double wfoam() const;

        double wsalt() const;

        bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                                 const WellState& well_state,
                                 Ewoms::DeferredLogger& deferred_logger) const;

        double getTHPConstraint(const SummaryState& summaryState) const;

        // Component fractions for each phase for the well
        const std::vector<double>& compFrac() const;

        double mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const;

        struct RatioLimitCheckReport;

        void checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                              const WellState& well_state,
                              RatioLimitCheckReport& report) const;

        void checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                              const WellState& well_state,
                              RatioLimitCheckReport& report) const;

        void checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                              const WellState& well_state,
                              RatioLimitCheckReport& report) const;

        void checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                                  const WellState& well_state,
                                  RatioLimitCheckReport& report,
                                  Ewoms::DeferredLogger& deferred_logger) const;

        template <typename RatioFunc>
        bool checkMaxRatioLimitWell(const WellState& well_state,
                                    const double max_ratio_limit,
                                    const RatioFunc& ratioFunc) const;

        template <typename RatioFunc>
        void checkMaxRatioLimitCompletions(const WellState& well_state,
                                           const double max_ratio_limit,
                                           const RatioFunc& ratioFunc,
                                           RatioLimitCheckReport& report) const;

        double scalingFactor(const int comp_idx) const;

        // whether a well is specified with a non-zero and valid VFP table number
        bool isVFPActive(Ewoms::DeferredLogger& deferred_logger) const;

        struct OperabilityStatus;

        OperabilityStatus operability_status_;

        // check whether the well is operable under BHP limit with current reservoir condition
        virtual void checkOperabilityUnderBHPLimitProducer(const WellState& well_state, const Simulator& eebos_simulator, Ewoms::DeferredLogger& deferred_logger) =0;

        // check whether the well is operable under THP limit with current reservoir condition
        virtual void checkOperabilityUnderTHPLimitProducer(const Simulator& eebos_simulator, const WellState& well_state, Ewoms::DeferredLogger& deferred_logger) =0;

        virtual void updateIPR(const Simulator& eebos_simulator, Ewoms::DeferredLogger& deferred_logger) const=0;

        void wellTestingEconomic(const Simulator& simulator, const std::vector<double>& B_avg,
                                 const double simulation_time, const WellState& well_state,
                                 WellTestState& welltest_state, Ewoms::DeferredLogger& deferred_logger);

        void wellTestingPhysical(const Simulator& simulator, const std::vector<double>& B_avg,
                                 const double simulation_time, const int report_step,
                                 WellState& well_state, WellTestState& welltest_state, Ewoms::DeferredLogger& deferred_logger);

        virtual void assembleWellEqWithoutIteration(const Simulator& eebosSimulator,
                                                    const double dt,
                                                    const Well::InjectionControls& inj_controls,
                                                    const Well::ProductionControls& prod_controls,
                                                    WellState& well_state,
                                                    Ewoms::DeferredLogger& deferred_logger) = 0;

        // iterate well equations with the specified control until converged
        virtual bool iterateWellEqWithControl(const Simulator& eebosSimulator,
                                              const std::vector<double>& B_avg,
                                              const double dt,
                                              const Well::InjectionControls& inj_controls,
                                              const Well::ProductionControls& prod_controls,
                                              WellState& well_state,
                                              Ewoms::DeferredLogger& deferred_logger) = 0;

        bool iterateWellEquations(const Simulator& eebosSimulator,
                                  const std::vector<double>& B_avg,
                                  const double dt,
                                  WellState& well_state,
                                  Ewoms::DeferredLogger& deferred_logger);

        void updateWellTestStateEconomic(const WellState& well_state,
                                         const double simulation_time,
                                         const bool write_message_to_opmlog,
                                         WellTestState& well_test_state,
                                         Ewoms::DeferredLogger& deferred_logger) const;

        void updateWellTestStatePhysical(const WellState& well_state,
                                         const double simulation_time,
                                         const bool write_message_to_opmlog,
                                         WellTestState& well_test_state,
                                         Ewoms::DeferredLogger& deferred_logger) const;

        void solveWellForTesting(const Simulator& eebosSimulator, WellState& well_state,
                                 const std::vector<double>& B_avg,
                                 Ewoms::DeferredLogger& deferred_logger);

        void initCompletions();

        bool checkConstraints(WellState& well_state,
                              const Schedule& schedule,
                              const SummaryState& summaryState,
                              DeferredLogger& deferred_logger) const;

        bool checkIndividualConstraints(WellState& well_state,
                                        const SummaryState& summaryState) const;

        bool checkGroupConstraints(WellState& well_state,
                                   const Schedule& schedule,
                                   const SummaryState& summaryState,
                                   DeferredLogger& deferred_logger) const;

        std::pair<bool, double> checkGroupConstraintsProd(const Group& group,
                                       const WellState& well_state,
                                       const double efficiencyFactor,
                                       const Schedule& schedule,
                                       const SummaryState& summaryState,
                                       DeferredLogger& deferred_logger) const;

        std::pair<bool, double> checkGroupConstraintsInj(const Group& group,
                                      const WellState& well_state,
                                      const double efficiencyFactor,
                                      const Schedule& schedule,
                                      const SummaryState& summaryState,
                                      DeferredLogger& deferred_logger) const;

        template <class EvalWell>
        void getGroupInjectionControl(const Group& group,
                                      const WellState& well_state,
                                      const Ewoms::Schedule& schedule,
                                      const SummaryState& summaryState,
                                      const InjectorType& injectorType,
                                      const EvalWell& bhp,
                                      const EvalWell& injection_rate,
                                      EvalWell& control_eq,
                                      double efficiencyFactor);

        template <class EvalWell>
        void getGroupProductionControl(const Group& group,
                                       const WellState& well_state,
                                       const Ewoms::Schedule& schedule,
                                       const SummaryState& summaryState,
                                       const EvalWell& bhp,
                                       const std::vector<EvalWell>& rates,
                                       EvalWell& control_eq,
                                       double efficiencyFactor);

        template <class EvalWell, class BhpFromThpFunc>
        void assembleControlEqInj(const WellState& well_state,
                                  const Ewoms::Schedule& schedule,
                                  const SummaryState& summaryState,
                                  const Well::InjectionControls& controls,
                                  const EvalWell& bhp,
                                  const EvalWell& injection_rate,
                                  BhpFromThpFunc bhp_from_thp,
                                  EvalWell& control_eq,
                                  Ewoms::DeferredLogger& deferred_logger);

        template <class EvalWell, class BhpFromThpFunc>
        void assembleControlEqProd(const WellState& well_state,
                                   const Ewoms::Schedule& schedule,
                                   const SummaryState& summaryState,
                                   const Well::ProductionControls& controls,
                                   const EvalWell& bhp,
                                   const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                                   BhpFromThpFunc bhp_from_thp,
                                   EvalWell& control_eq,
                                   Ewoms::DeferredLogger& deferred_logger);
    };

    // definition of the struct OperabilityStatus
    template<typename TypeTag>
    struct
    WellInterface<TypeTag>::
    OperabilityStatus {
        bool isOperable() const {
            if (!operable_under_only_bhp_limit) {
                return false;
            } else {
                return ( (isOperableUnderBHPLimit() || isOperableUnderTHPLimit()) );
            }
        }

        bool isOperableUnderBHPLimit() const {
            return operable_under_only_bhp_limit && obey_thp_limit_under_bhp_limit;
        }

        bool isOperableUnderTHPLimit() const {
            return can_obtain_bhp_with_thp_limit && obey_bhp_limit_with_thp_limit;
        }

        void reset() {
            operable_under_only_bhp_limit = true;
            obey_thp_limit_under_bhp_limit = true;
            can_obtain_bhp_with_thp_limit = true;
            obey_bhp_limit_with_thp_limit = true;
        }

        // whether the well can be operated under bhp limit
        // without considering other limits.
        // if it is false, then the well is not operable for sure.
        bool operable_under_only_bhp_limit = true;
        // if the well can be operated under bhp limit, will it obey(not violate)
        // the thp limit when operated under bhp limit
        bool obey_thp_limit_under_bhp_limit = true;
        // whether the well operate under the thp limit only
        bool can_obtain_bhp_with_thp_limit = true;
        // whether the well obey bhp limit when operated under thp limit
        bool obey_bhp_limit_with_thp_limit = true;

    };

    template<typename TypeTag>
    struct
    WellInterface<TypeTag>::
    RatioLimitCheckReport{
        bool ratio_limit_violated = false;
        int worst_offending_completion = INVALIDCOMPLETION;
        double violation_extent = 0.0;
    };

    const std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };

}

#include "wellinterface_impl.hh"

#endif // EWOMS_WELLINTERFACE_HH
