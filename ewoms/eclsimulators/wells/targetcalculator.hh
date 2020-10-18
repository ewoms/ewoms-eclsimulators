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
#ifndef EWOMS_TARGETCALCULATOR_HH
#define EWOMS_TARGETCALCULATOR_HH

#include <ewoms/eclio/parser/eclipsestate/schedule/group/group.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/group/guiderate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclsimulators/utils/deferredlogger.hh>
#include <ewoms/eclsimulators/utils/deferredloggingerrorhelpers.hh>
#include <ewoms/eclsimulators/wells/wellstatefullyimplicitblackoil.hh>

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <vector>

namespace Ewoms
{

namespace WellGroupHelpers
{

    /// Based on a group control mode, extract or calculate rates, and
    /// provide other conveniences.
    class TargetCalculator
    {
    public:
        TargetCalculator(const Group::ProductionCMode cmode,
                         const PhaseUsage& pu,
                         const std::vector<double>& resv_coeff,
                         const double group_grat_target_from_sales)
            : cmode_(cmode)
            , pu_(pu)
            , resv_coeff_(resv_coeff)
            , group_grat_target_from_sales_(group_grat_target_from_sales)
        {

        }

        template <typename RateVec>
        auto calcModeRateFromRates(const RateVec& rates) const
        {
            // ElemType is just the plain element type of the rates container,
            // without any reference, const or volatile modifiers.
            using ElemType = std::remove_cv_t<std::remove_reference_t<decltype(rates[0])>>;
            switch (cmode_) {
            case Group::ProductionCMode::ORAT: {
                assert(pu_.phase_used[BlackoilPhases::Liquid]);
                const int pos = pu_.phase_pos[BlackoilPhases::Liquid];
                return rates[pos];
            }
            case Group::ProductionCMode::WRAT: {
                assert(pu_.phase_used[BlackoilPhases::Aqua]);
                const int pos = pu_.phase_pos[BlackoilPhases::Aqua];
                return rates[pos];
            }
            case Group::ProductionCMode::GRAT: {
                assert(pu_.phase_used[BlackoilPhases::Vapour]);
                const int pos = pu_.phase_pos[BlackoilPhases::Vapour];
                return rates[pos];
            }
            case Group::ProductionCMode::LRAT: {
                assert(pu_.phase_used[BlackoilPhases::Liquid]);
                assert(pu_.phase_used[BlackoilPhases::Aqua]);
                const int opos = pu_.phase_pos[BlackoilPhases::Liquid];
                const int wpos = pu_.phase_pos[BlackoilPhases::Aqua];
                return rates[opos] + rates[wpos];
            }
            case Group::ProductionCMode::RESV: {
                assert(pu_.phase_used[BlackoilPhases::Liquid]);
                assert(pu_.phase_used[BlackoilPhases::Aqua]);
                assert(pu_.phase_used[BlackoilPhases::Vapour]);
                ElemType mode_rate = zero<ElemType>();
                for (int phase = 0; phase < pu_.num_phases; ++phase) {
                    mode_rate += rates[phase] * resv_coeff_[phase];
                }
                return mode_rate;
            }
            default:
                // Should never be here.
                assert(false);
                return zero<ElemType>();
            }
        }

        double groupTarget(const Group::ProductionControls ctrl) const
        {
            switch (cmode_) {
            case Group::ProductionCMode::ORAT:
                return ctrl.oil_target;
            case Group::ProductionCMode::WRAT:
                return ctrl.water_target;
            case Group::ProductionCMode::GRAT:
            {
                // gas target may have been adjusted by GCONSALE
                if ( group_grat_target_from_sales_ > 0)
                    return group_grat_target_from_sales_;

                return ctrl.gas_target;
            }
            case Group::ProductionCMode::LRAT:
                return ctrl.liquid_target;
            case Group::ProductionCMode::RESV:
                return ctrl.resv_target;
            default:
                // Should never be here.
                assert(false);
                return 0.0;
            }
        }

        GuideRateModel::Target guideTargetMode() const
        {
            switch (cmode_) {
            case Group::ProductionCMode::ORAT:
                return GuideRateModel::Target::OIL;
            case Group::ProductionCMode::WRAT:
                return GuideRateModel::Target::WAT;
            case Group::ProductionCMode::GRAT:
                return GuideRateModel::Target::GAS;
            case Group::ProductionCMode::LRAT:
                return GuideRateModel::Target::LIQ;
            case Group::ProductionCMode::RESV:
                return GuideRateModel::Target::RES;
            default:
                // Should never be here.
                assert(false);
                return GuideRateModel::Target::NONE;
            }
        }

    private:
        template <typename ElemType>
        static ElemType zero()
        {
            // This is for Evaluation types.
            ElemType x;
            x = 0.0;
            return x;
        }
        Group::ProductionCMode cmode_;
        const PhaseUsage& pu_;
        const std::vector<double>& resv_coeff_;
        const double group_grat_target_from_sales_;
    };

} // namespace WellGroupHelpers

} // namespace Ewoms

#endif
