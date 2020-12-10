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
#ifndef EWOMS_WELLGROUPHELPERS_HH
#define EWOMS_WELLGROUPHELPERS_HH

#include <ewoms/eclio/parser/eclipsestate/schedule/group/group.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/guiderate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclsimulators/utils/deferredlogger.hh>
#include <ewoms/eclsimulators/utils/deferredloggingerrorhelpers.hh>
#include <ewoms/eclsimulators/wells/vfpprodproperties.hh>
#include <ewoms/eclsimulators/wells/wellstatefullyimplicitblackoil.hh>

#include <algorithm>
#include <cassert>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

namespace Ewoms
{

namespace WellGroupHelpers
{

    void setCmodeGroup(const Group& group,
                       const Schedule& schedule,
                       const SummaryState& summaryState,
                       const int reportStepIdx,
                       WellStateFullyImplicitBlackoil& wellState);

    void accumulateGroupEfficiencyFactor(const Group& group,
                                         const Schedule& schedule,
                                         const int reportStepIdx,
                                         double& factor);

    double sumWellPhaseRates(const std::vector<double>& rates,
                             const Group& group,
                             const Schedule& schedule,
                             const WellStateFullyImplicitBlackoil& wellState,
                             const int reportStepIdx,
                             const int phasePos,
                             const bool injector);

    double sumWellRates(const Group& group,
                        const Schedule& schedule,
                        const WellStateFullyImplicitBlackoil& wellState,
                        const int reportStepIdx,
                        const int phasePos,
                        const bool injector);

    double sumWellResRates(const Group& group,
                           const Schedule& schedule,
                           const WellStateFullyImplicitBlackoil& wellState,
                           const int reportStepIdx,
                           const int phasePos,
                           const bool injector);

    double sumSolventRates(const Group& group,
                           const Schedule& schedule,
                           const WellStateFullyImplicitBlackoil& wellState,
                           const int reportStepIdx,
                           const bool injector);

    void updateGroupTargetReduction(const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    const bool isInjector,
                                    const PhaseUsage& pu,
                                    const GuideRate& guide_rate,
                                    const WellStateFullyImplicitBlackoil& wellStateNupcol,
                                    WellStateFullyImplicitBlackoil& wellState,
                                    std::vector<double>& groupTargetReduction);

    template <class Comm>
    void updateGuideRateForGroups(const Group& group,
                                  const Schedule& schedule,
                                  const PhaseUsage& pu,
                                  const int reportStepIdx,
                                  const double& simTime,
                                  const bool isInjector,
                                  WellStateFullyImplicitBlackoil& wellState,
                                  const Comm& comm,
                                  GuideRate* guideRate,
                                  std::vector<double>& pot)
    {
        const int np = pu.num_phases;
        for (const std::string& groupName : group.groups()) {
            std::vector<double> thisPot(np, 0.0);
            const Group& groupTmp = schedule.getGroup(groupName, reportStepIdx);
            updateGuideRateForGroups(
                groupTmp, schedule, pu, reportStepIdx, simTime, isInjector, wellState, comm, guideRate, thisPot);

            const auto gefac = groupTmp.getGroupEfficiencyFactor();

            // accumulate group contribution from sub group unconditionally
            if (isInjector) {
                const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
                for (Phase phase : all) {
                    int phasePos;
                    if (phase == Phase::GAS && pu.phase_used[BlackoilPhases::Vapour])
                        phasePos = pu.phase_pos[BlackoilPhases::Vapour];
                    else if (phase == Phase::OIL && pu.phase_used[BlackoilPhases::Liquid])
                        phasePos = pu.phase_pos[BlackoilPhases::Liquid];
                    else if (phase == Phase::WATER && pu.phase_used[BlackoilPhases::Aqua])
                        phasePos = pu.phase_pos[BlackoilPhases::Aqua];
                    else
                        continue;

                    pot[phasePos] += gefac * thisPot[phasePos];
                }
            } else {
                const Group::ProductionCMode& currentGroupControl = wellState.currentProductionGroupControl(groupName);
                if (currentGroupControl != Group::ProductionCMode::FLD
                    && currentGroupControl != Group::ProductionCMode::NONE) {
                    continue;
                }
                for (int phase = 0; phase < np; phase++) {
                    pot[phase] += gefac * thisPot[phase];
                }
            }
        }
        for (const std::string& wellName : group.wells()) {
            const auto& wellTmp = schedule.getWell(wellName, reportStepIdx);
            const auto wefac = wellTmp.getEfficiencyFactor();

            if (wellTmp.isProducer() && isInjector)
                continue;

            if (wellTmp.isInjector() && !isInjector)
                continue;

            if (wellTmp.getStatus() == Well::Status::SHUT)
                continue;
            const auto& end = wellState.wellMap().end();
            const auto& it = wellState.wellMap().find(wellName);
            if (it == end) // the well is not found
                continue;

            int well_index = it->second[0];

            if (! wellState.wellIsOwned(well_index, wellName) ) // Only sum once
            {
                continue;
            }

            const auto wellrate_index = well_index * wellState.numPhases();
            // add contribution from wells unconditionally
            for (int phase = 0; phase < np; phase++) {
                pot[phase] += wefac * wellState.wellPotentials()[wellrate_index + phase];
            }
        }

        double oilPot = 0.0;
        if (pu.phase_used[BlackoilPhases::Liquid])
            oilPot = pot[pu.phase_pos[BlackoilPhases::Liquid]];

        double gasPot = 0.0;
        if (pu.phase_used[BlackoilPhases::Vapour])
            gasPot = pot[pu.phase_pos[BlackoilPhases::Vapour]];

        double waterPot = 0.0;
        if (pu.phase_used[BlackoilPhases::Aqua])
            waterPot = pot[pu.phase_pos[BlackoilPhases::Aqua]];

        oilPot = comm.sum(oilPot);
        gasPot = comm.sum(gasPot);
        waterPot = comm.sum(waterPot);

        if (isInjector) {
            wellState.setCurrentGroupInjectionPotentials(group.name(), pot);
        } else {
            guideRate->compute(group.name(), reportStepIdx, simTime, oilPot, gasPot, waterPot);
        }
    }

    template <class Comm>
    void updateGuideRatesForWells(const Schedule& schedule,
                                  const PhaseUsage& pu,
                                  const int reportStepIdx,
                                  const double& simTime,
                                  const WellStateFullyImplicitBlackoil& wellState,
                                  const Comm& comm,
                                  GuideRate* guideRate)
    {

        const auto& end = wellState.wellMap().end();
        for (const auto& well : schedule.getWells(reportStepIdx)) {
            double oilpot = 0.0;
            double gaspot = 0.0;
            double waterpot = 0.0;

            const auto& it = wellState.wellMap().find(well.name());
            if (it != end && wellState.wellIsOwned(it->second[0], well.name()))
            {
                // the well is found and owned
                int well_index = it->second[0];

                const auto wpot = wellState.wellPotentials().data() + well_index * wellState.numPhases();
                if (pu.phase_used[BlackoilPhases::Liquid] > 0)
                    oilpot = wpot[pu.phase_pos[BlackoilPhases::Liquid]];

                if (pu.phase_used[BlackoilPhases::Vapour] > 0)
                    gaspot = wpot[pu.phase_pos[BlackoilPhases::Vapour]];

                if (pu.phase_used[BlackoilPhases::Aqua] > 0)
                    waterpot = wpot[pu.phase_pos[BlackoilPhases::Aqua]];
            }
            oilpot = comm.sum(oilpot);
            gaspot = comm.sum(gaspot);
            waterpot = comm.sum(waterpot);
            guideRate->compute(well.name(), reportStepIdx, simTime, oilpot, gaspot, waterpot);
        }
    }

    void updateVREPForGroups(const Group& group,
                             const Schedule& schedule,
                             const int reportStepIdx,
                             const WellStateFullyImplicitBlackoil& wellStateNupcol,
                             WellStateFullyImplicitBlackoil& wellState);

    void updateReservoirRatesInjectionGroups(const Group& group,
                                             const Schedule& schedule,
                                             const int reportStepIdx,
                                             const WellStateFullyImplicitBlackoil& wellStateNupcol,
                                             WellStateFullyImplicitBlackoil& wellState);

    void updateWellRates(const Group& group,
                         const Schedule& schedule,
                         const int reportStepIdx,
                         const WellStateFullyImplicitBlackoil& wellStateNupcol,
                         WellStateFullyImplicitBlackoil& wellState);

    void updateGroupProductionRates(const Group& group,
                                    const Schedule& schedule,
                                    const int reportStepIdx,
                                    const WellStateFullyImplicitBlackoil& wellStateNupcol,
                                    WellStateFullyImplicitBlackoil& wellState);

    void updateREINForGroups(const Group& group,
                             const Schedule& schedule,
                             const int reportStepIdx,
                             const PhaseUsage& pu,
                             const SummaryState& st,
                             const WellStateFullyImplicitBlackoil& wellStateNupcol,
                             WellStateFullyImplicitBlackoil& wellState);

    std::map<std::string, double>
    computeNetworkPressures(const Ewoms::Network::ExtNetwork& network,
                            const WellStateFullyImplicitBlackoil& well_state,
                            const VFPProdProperties& vfp_prod_props,
                            const Schedule& schedule,
                            const int report_time_step);

    GuideRate::RateVector
    getRateVector(const WellStateFullyImplicitBlackoil& well_state, const PhaseUsage& pu, const std::string& name);

    GuideRate::RateVector
    getProductionGroupRateVector(const WellStateFullyImplicitBlackoil& well_state, const PhaseUsage& pu, const std::string& group_name);

    double getGuideRate(const std::string& name,
                        const Schedule& schedule,
                        const WellStateFullyImplicitBlackoil& wellState,
                        const int reportStepIdx,
                        const GuideRate* guideRate,
                        const GuideRateModel::Target target,
                        const PhaseUsage& pu);

    double getGuideRateInj(const std::string& name,
                           const Schedule& schedule,
                           const WellStateFullyImplicitBlackoil& wellState,
                           const int reportStepIdx,
                           const GuideRate* guideRate,
                           const GuideRateModel::Target target,
                           const Phase& injectionPhase,
                           const PhaseUsage& pu);

    int groupControlledWells(const Schedule& schedule,
                             const WellStateFullyImplicitBlackoil& well_state,
                             const int report_step,
                             const std::string& group_name,
                             const std::string& always_included_child);

    class FractionCalculator
    {
    public:
        FractionCalculator(const Schedule& schedule,
                           const WellStateFullyImplicitBlackoil& well_state,
                           const int report_step,
                           const GuideRate* guide_rate,
                           const GuideRateModel::Target target,
                           const PhaseUsage& pu);
        double fraction(const std::string& name, const std::string& control_group_name, const bool always_include_this);
        double localFraction(const std::string& name, const std::string& always_included_child);

    private:
        std::string parent(const std::string& name);
        double guideRateSum(const Group& group, const std::string& always_included_child);
        double guideRate(const std::string& name, const std::string& always_included_child);
        int groupControlledWells(const std::string& group_name, const std::string& always_included_child);
        GuideRate::RateVector getGroupRateVector(const std::string& group_name);
        const Schedule& schedule_;
        const WellStateFullyImplicitBlackoil& well_state_;
        int report_step_;
        const GuideRate* guide_rate_;
        GuideRateModel::Target target_;
        PhaseUsage pu_;
    };

    double fractionFromGuideRates(const std::string& name,
                                  const std::string& controlGroupName,
                                  const Schedule& schedule,
                                  const WellStateFullyImplicitBlackoil& wellState,
                                  const int reportStepIdx,
                                  const GuideRate* guideRate,
                                  const GuideRateModel::Target target,
                                  const PhaseUsage& pu,
                                  const bool alwaysIncludeThis = false);

    double fractionFromInjectionPotentials(const std::string& name,
                                           const std::string& controlGroupName,
                                           const Schedule& schedule,
                                           const WellStateFullyImplicitBlackoil& wellState,
                                           const int reportStepIdx,
                                           const GuideRate* guideRate,
                                           const GuideRateModel::Target target,
                                           const PhaseUsage& pu,
                                           const Phase& injectionPhase,
                                           const bool alwaysIncludeThis = false);

    std::pair<bool, double> checkGroupConstraintsInj(const std::string& name,
                                                     const std::string& parent,
                                                     const Group& group,
                                                     const WellStateFullyImplicitBlackoil& wellState,
                                                     const int reportStepIdx,
                                                     const GuideRate* guideRate,
                                                     const double* rates,
                                                     Phase injectionPhase,
                                                     const PhaseUsage& pu,
                                                     const double efficiencyFactor,
                                                     const Schedule& schedule,
                                                     const SummaryState& summaryState,
                                                     const std::vector<double>& resv_coeff,
                                                     DeferredLogger& deferred_logger);

    std::vector<std::string> groupChainTopBot(const std::string& bottom,
                                              const std::string& top,
                                              const Schedule& schedule,
                                              const int report_step);

    std::pair<bool, double> checkGroupConstraintsProd(const std::string& name,
                                                      const std::string& parent,
                                                      const Group& group,
                                                      const WellStateFullyImplicitBlackoil& wellState,
                                                      const int reportStepIdx,
                                                      const GuideRate* guideRate,
                                                      const double* rates,
                                                      const PhaseUsage& pu,
                                                      const double efficiencyFactor,
                                                      const Schedule& schedule,
                                                      const SummaryState& summaryState,
                                                      const std::vector<double>& resv_coeff,
                                                      DeferredLogger& deferred_logger);

} // namespace WellGroupHelpers

} // namespace Ewoms

#endif
