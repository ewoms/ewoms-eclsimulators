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
#ifndef EWOMS_SIMULATORREPORT_HH
#define EWOMS_SIMULATORREPORT_HH
#include <cassert>
#include <iosfwd>
#include <vector>

namespace Ewoms
{

    /// A struct for returning timing data from a simulator to its caller.
    struct SimulatorReportSingle
    {
        double pressure_time;
        double transport_time;
        double total_time;
        double solver_time;
        double assemble_time;
        double pre_post_time;
        double assemble_time_well;
        double linear_solve_setup_time;
        double linear_solve_time;
        double update_time;
        double output_write_time;

        unsigned int total_well_iterations;
        unsigned int total_linearizations;
        unsigned int total_newton_iterations;
        unsigned int total_linear_iterations;

        bool converged;
        int exit_status;

        double global_time;
        double timestep_length;

        /// Default constructor initializing all times to 0.0.
        SimulatorReportSingle();
        /// Increment this report's times by those in sr.
        void operator+=(const SimulatorReportSingle& sr);
        /// Print a report suitable for a single simulation step.
        void reportStep(std::ostringstream& os) const;
        /// Print a report suitable for the end of a fully implicit case, leaving out the pressure/transport time.
        void reportFullyImplicit(std::ostream& os, const SimulatorReportSingle* failedReport = nullptr) const;
    };

    struct SimulatorReport
    {
        SimulatorReportSingle success;
        SimulatorReportSingle failure;
        std::vector<SimulatorReportSingle> stepreports;

        void operator+=(const SimulatorReportSingle& sr);
        void operator+=(const SimulatorReport& sr);
        void reportFullyImplicit(std::ostream& os) const;
        void fullReports(std::ostream& os) const;
    };

    } // namespace Ewoms

#endif // EWOMS_SIMULATORREPORT_HH
