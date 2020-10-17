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
#include <config.h>

#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/python/python.hh>
#include <ewoms/eclio/parser/units/units.hh>
#include <ewoms/eclsimulators/wells/vfpproperties.hh>
#include <ewoms/eclsimulators/wells/vfpinjproperties.hh>
#include <ewoms/eclsimulators/wells/vfpprodproperties.hh>

#include <iostream>
#include <iomanip>

using namespace Ewoms;

struct Setup
{
    using VFP = VFPProperties<VFPInjProperties, VFPProdProperties>;
    std::unique_ptr<const EclipseState> ecl_state;
    std::shared_ptr<Python> python;
    std::unique_ptr<const Schedule> schedule;
    std::unique_ptr<SummaryState> summary_state;
    std::unique_ptr<VFP> vfp_properties;

    Setup(const std::string& file)
    {
        Parser parser;
        auto deck = parser.parseFile(file);
        ecl_state.reset(new EclipseState(deck) );
        {
          const TableManager table( deck );
          const Runspec runspec(deck);
          python = std::make_shared<Python>();
          schedule.reset( new Schedule(deck, *ecl_state, python));
          summary_state.reset( new SummaryState(std::chrono::system_clock::from_time_t(schedule->getStartTime())));
        }
        const int step = 0;
        vfp_properties = std::make_unique<VFP>(schedule->getVFPInjTables(step), schedule->getVFPProdTables(step));
    };
};

double computeBhp(const VFPProdTable& table,
                  const double flo,
                  const double thp,
                  const double wfr,
                  const double gfr,
                  const double alq)
{

    // First, find the values to interpolate between.
    // Assuming positive flo here!
    assert(flo > 0.0);
    auto flo_i = detail::findInterpData(flo, table.getFloAxis());
    auto thp_i = detail::findInterpData(thp, table.getTHPAxis()); // assume constant
    auto wfr_i = detail::findInterpData(wfr, table.getWFRAxis());
    auto gfr_i = detail::findInterpData(gfr, table.getGFRAxis());
    auto alq_i = detail::findInterpData(alq, table.getALQAxis()); //assume constant

    return detail::interpolate(table, flo_i, thp_i, wfr_i, gfr_i, alq_i).value;
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        return EXIT_FAILURE;
    }
    Setup setup(argv[1]);

//    const int table_id = 1;
    const int table_id = 4;
    const double wct = 0.0;
    const double gor = 35.2743;
    const double alq = 0.0;
    const int n = 51;
    const double m3pd = unit::cubic(unit::meter)/unit::day;
    const double rate_min = 20.0 * m3pd;
    const double rate_max = 2000.0 * m3pd;
    // const double thp = 32.1744 * unit::barsa;
    // const double thp = 10.0 * unit::barsa;
    const double thp_min = 10.0 * unit::barsa;
    const double thp_max = 35.0 * unit::barsa;
    std::vector<double> rates(n);
    std::vector<double> thps(n);
    for (int ii = 0; ii < n; ++ii) {
        const double q = double(ii) / double(n-1);
        rates[ii] = (1.0 - q) * rate_min + q * rate_max;
        thps[ii] = (1.0 - q) * thp_min + q * thp_max;
    }

    const VFPProdTable& table = *(setup.vfp_properties->getProd()->getTable(table_id));
    std::cout.precision(12);
    for (double rate : rates) {
        for (double thp : thps) {
            const double bhp = computeBhp(table, rate, thp, wct, gor, alq);
            std::cout //<< std::setw(18) << unit::convert::to(rate, m3pd)
                      //<< std::setw(18) << unit::convert::to(thp, unit::barsa)
                      << std::setw(18) << unit::convert::to(bhp, unit::barsa) << '\n';
        }
    }

    return EXIT_SUCCESS;
}
