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
#ifndef EWOMS_SIM_FIBO_DETAILS_HH
#define EWOMS_SIM_FIBO_DETAILS_HH

#include <utility>
#include <algorithm>
#include <locale>

#include <ewoms/eclio/parser/eclipsestate/schedule/events.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/well.hh>

namespace Ewoms
{
    namespace SimFIBODetails {
        typedef std::unordered_map<std::string, Well > WellMap;

        inline WellMap
        mapWells(const std::vector< Well >& wells)
        {
            WellMap wmap;

            for (const auto& w : wells)
            {
                wmap.insert(std::make_pair(w.name(), w));
            }

            return wmap;
        }

        inline void
        historyRates(const PhaseUsage&               pu,
                     const Well::ProductionControls& p,
                     std::vector<double>&            rates)
        {
            assert (! p.prediction_mode);
            assert (rates.size() ==
                    std::vector<double>::size_type(pu.num_phases));

            if (pu.phase_used[ BlackoilPhases::Aqua ]) {
                const std::vector<double>::size_type
                    i = pu.phase_pos[ BlackoilPhases::Aqua ];

                rates[i] = p.water_rate;
            }

            if (pu.phase_used[ BlackoilPhases::Liquid ]) {
                const std::vector<double>::size_type
                    i = pu.phase_pos[ BlackoilPhases::Liquid ];

                rates[i] = p.oil_rate;
            }

            if (pu.phase_used[ BlackoilPhases::Vapour ]) {
                const std::vector<double>::size_type
                    i = pu.phase_pos[ BlackoilPhases::Vapour ];

                rates[i] = p.gas_rate;
            }
        }
    } // namespace SimFIBODetails
} // namespace Ewoms

#endif // EWOMS_SIM_FIBO_DETAILS_HH
