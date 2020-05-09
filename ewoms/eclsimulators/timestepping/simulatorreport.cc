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

#include "config.h"

#include <ewoms/eclsimulators/timestepping/simulatorreport.hh>
#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/parser/units/units.hh>

#include <iomanip>
#include <ostream>
#include <sstream>

namespace Ewoms
{
    SimulatorReportSingle::SimulatorReportSingle()
        : pressure_time(0.0),
          transport_time(0.0),
          total_time(0.0),
          solver_time(0.0),
          assemble_time(0.0),
          linear_solve_setup_time(0.0),
          linear_solve_time(0.0),
          update_time(0.0),
          output_write_time(0.0),
          total_well_iterations(0),
          total_linearizations( 0 ),
          total_newton_iterations( 0 ),
          total_linear_iterations( 0 ),
          converged(false),
          exit_status(EXIT_SUCCESS),
          global_time(0),
          timestep_length(0.0)
    {
    }

    void SimulatorReportSingle::operator+=(const SimulatorReportSingle& sr)
    {
        pressure_time += sr.pressure_time;
        transport_time += sr.transport_time;
        linear_solve_setup_time += sr.linear_solve_setup_time;
        linear_solve_time += sr.linear_solve_time;
        solver_time += sr.solver_time;
        assemble_time += sr.assemble_time;
        update_time += sr.update_time;
        output_write_time += sr.output_write_time;
        total_time += sr.total_time;
        total_well_iterations += sr.total_well_iterations;
        total_linearizations += sr.total_linearizations;
        total_newton_iterations += sr.total_newton_iterations;
        total_linear_iterations += sr.total_linear_iterations;
        global_time = sr.global_time; // It makes no sense adding time points, so = not += here.
    }

    void SimulatorReportSingle::reportStep(std::ostringstream& ss) const
    {
        ss << "Time step summary: ";
        if (total_well_iterations != 0) {
            ss << "well its = " << std::setw(2) << total_well_iterations << ", ";
        }
        ss << "newton its = " << std::setw(2) << total_newton_iterations << ", "
           << "linearizations = "  << std::setw(2) << total_linearizations
           << " ("  << std::fixed << std::setprecision(3) << std::setw(6) << assemble_time << " sec), "
           << "linear its = " << std::setw(3) << total_linear_iterations
           << " ("  << std::fixed << std::setprecision(3) << std::setw(6) << linear_solve_time << " sec)";
    }

    void SimulatorReportSingle::reportFullyImplicit(std::ostream& os, const SimulatorReportSingle* failureReport) const
    {
        double t = total_time;
        os << "Total time (seconds):         " << t;
        os << std::endl;

        t = solver_time + (failureReport ? failureReport->solver_time : 0.0);
        os << "Solver time (seconds):        " << t;
        os << std::endl;

        if (assemble_time > 0.0 || linear_solve_time > 0.0) {
            t = assemble_time + (failureReport ? failureReport->assemble_time : 0.0);
            os << " Assembly time (seconds):     " << t;
            if (failureReport) {
                os << " (Failed: " << failureReport->assemble_time << "; "
                   << 100*failureReport->assemble_time/t << "%)";
            }
            os << std::endl;

            t = linear_solve_time + (failureReport ? failureReport->linear_solve_time : 0.0);
            os << " Linear solve time (seconds): " << t;
            if (failureReport) {
                os << " (Failed: " << failureReport->linear_solve_time << "; "
                   << 100*failureReport->linear_solve_time/t << "%)";
            }
            os << std::endl;

            t = linear_solve_setup_time + (failureReport ? failureReport->linear_solve_setup_time : 0.0);
            os << " Linear solve setup time (seconds): " << t;
            if (failureReport) {
                os << " (Failed: " << failureReport->linear_solve_setup_time << "; "
                   << 100*failureReport->linear_solve_setup_time/t << "%)";
            }
            os << std::endl;

            t = update_time + (failureReport ? failureReport->update_time : 0.0);
            os << " Update time (seconds):       " << t;
            if (failureReport) {
                os << " (Failed: " << failureReport->update_time << "; "
                   << 100*failureReport->update_time/t << "%)";
            }
            os << std::endl;

            t = output_write_time + (failureReport ? failureReport->output_write_time : 0.0);
            os << " Output write time (seconds): " << t;
            os << std::endl;

        }

        int n = total_well_iterations + (failureReport ? failureReport->total_well_iterations : 0);
        os << "Overall Well Iterations:      " << n;
        if (failureReport) {
            os << " (Failed: " << failureReport->total_well_iterations << "; "
               << 100.0*failureReport->total_well_iterations/n << "%)";
        }
        os << std::endl;

        n = total_linearizations + (failureReport ? failureReport->total_linearizations : 0);
        os << "Overall Linearizations:       " << n;
        if (failureReport) {
            os << " (Failed: " << failureReport->total_linearizations << "; "
               << 100.0*failureReport->total_linearizations/n << "%)";
        }
        os << std::endl;

        n = total_newton_iterations + (failureReport ? failureReport->total_newton_iterations : 0);
        os << "Overall Newton Iterations:    " << n;
        if (failureReport) {
            os << " (Failed: " << failureReport->total_newton_iterations << "; "
               << 100.0*failureReport->total_newton_iterations/n << "%)";
        }
        os << std::endl;

        n = total_linear_iterations + (failureReport ? failureReport->total_linear_iterations : 0);
        os << "Overall Linear Iterations:    " << n;
        if (failureReport) {
            os << " (Failed: " << failureReport->total_linear_iterations << "; "
               << 100.0*failureReport->total_linear_iterations/n << "%)";
        }
        os << std::endl;
    }

    void SimulatorReport::operator+=(const SimulatorReportSingle& sr)
    {
        if (sr.converged) {
            success += sr;
        } else {
            failure += sr;
        }
        stepreports.push_back(sr);
    }

    void SimulatorReport::operator+=(const SimulatorReport& sr)
    {
        success += sr.success;
        failure += sr.failure;
        stepreports.insert(stepreports.end(), sr.stepreports.begin(), sr.stepreports.end());
    }

    void SimulatorReport::reportFullyImplicit(std::ostream& os) const
    {
        success.reportFullyImplicit(os, &failure);
    }

    void SimulatorReport::fullReports(std::ostream& os) const
    {
        os << "  Time(day)  TStep(day)  Assembly    LSolve    LSetup    Update    Output WellIt Lins NewtIt LinIt Conv\n";
        for (size_t i = 0; i < this->stepreports.size(); ++i) {
            const SimulatorReportSingle& sr = this->stepreports[i];
            os.precision(10);
            os << std::defaultfloat;
            os << std::setw(11) << unit::convert::to(sr.global_time, unit::day) << " ";
            os << std::setw(11) << unit::convert::to(sr.timestep_length, unit::day) << " ";
            os.precision(4);
            os << std::fixed;
            os << std::setw(9) << sr.assemble_time << " ";
            os << std::setw(9) << sr.linear_solve_setup_time << " ";
            os << std::setw(9) << sr.linear_solve_time << " ";
            os << std::setw(9) << sr.update_time << " ";
            os << std::setw(9) << sr.output_write_time << " ";
            os.precision(6);
            os << std::defaultfloat;
            os << std::setw(6) << sr.total_well_iterations << " ";
            os << std::setw(4) << sr.total_linearizations << " ";
            os << std::setw(6) << sr.total_newton_iterations << " ";
            os << std::setw(5) << sr.total_linear_iterations << " ";
            os << std::setw(4) << sr.converged << "\n";
        }
    }

} // namespace Ewoms
