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

#ifndef EWOMS_SIMULATORTIMER_HH
#define EWOMS_SIMULATORTIMER_HH

#include <ewoms/eclio/parser/eclipsestate/schedule/timemap.hh>
#include <ewoms/eclsimulators/timestepping/simulatortimerinterface.hh>

#include <iosfwd>
#include <vector>

namespace Ewoms
{
    class SimulatorTimer : public SimulatorTimerInterface
    {
    public:
        // use default implementation of these methods
        using SimulatorTimerInterface::currentDateTime;
        using SimulatorTimerInterface::currentPosixTime;

        /// Default constructor.
        SimulatorTimer();

        /// Use the SimulatorTimer as a shim around ewoms-eclio's Ewoms::TimeMap
        void init(const TimeMap& timeMap, size_t report_step = 0);

        /// Whether the current step is the first step.
        bool initialStep() const;

        /// Total number of steps.
        int numSteps() const;

        /// Current step number. This is the number of timesteps that
        /// has been completed from the start of the run. The time
        /// after initialization but before the simulation has started
        /// is timestep number zero.
        int currentStepNum() const;

        /// Set current step number.
        void setCurrentStepNum(int step);

        /// Current step length. This is the length of the step
        /// the simulator will take in the next iteration.
        ///
        /// @note if done(), it is an error to call currentStepLength().
        double currentStepLength() const;

        /// Previous step length. This is the length of the step that
        /// was taken to arrive at this time.
        ///
        /// @note if no increments have been done (i.e. the timer is
        /// still in its constructed state and currentStepNum() == 0),
        /// it is an error to call stepLengthTaken().
        double stepLengthTaken () const;

        /// Time elapsed since the start of the simulation until the
        /// beginning of the current time step [s].
        double simulationTimeElapsed() const;

        /// Total time.
        double totalTime() const;

        /// Return start date of simulation
        boost::posix_time::ptime startDateTime() const;

        /// Return current date.
        boost::posix_time::ptime currentDateTime() const;

        /// Set total time.
        /// This is primarily intended for multi-epoch schedules,
        /// where a timer for a given epoch does not have
        /// access to later timesteps.
        void setTotalTime(double time);

        /// Print a report with current and total time etc.
        /// Note: if done(), it is an error to call report().
        void report(std::ostream& os) const;

        /// advance time by currentStepLength
        SimulatorTimer& operator++();

        /// advance time by currentStepLength
        void advance() { this->operator++(); }

        /// Return true if op++() has been called numSteps() times.
        bool done() const;

        /// Always return false. Timestep failures is handled in the
        /// substepTimer
        bool lastStepFailed() const {return false;}

        /// return copy of object
        virtual std::unique_ptr< SimulatorTimerInterface > clone() const;

    private:
        std::vector<double> timesteps_;
        int current_step_;
        double current_time_;
        double total_time_;
        boost::gregorian::date start_date_;
    };

} // namespace Ewoms

#endif // EWOMS_SIMULATORTIMER_HH