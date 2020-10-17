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
#include "config.h"
#include <ewoms/eclsimulators/timestepping/simulatortimer.hh>
#include <ewoms/eclio/utility/parameters/parametergroup.hh>
#include <ewoms/eclio/parser/units/units.hh>
#include <ostream>
#include <numeric>

namespace Ewoms
{

    /// Default constructor.
    SimulatorTimer::SimulatorTimer()
        : current_step_(0),
          current_time_(0.0),
          start_date_(2012,1,1)    // A really arbitrary default starting value?!
    {
    }

    /// Use the SimulatorTimer as a shim around ewoms-eclio's Ewoms::TimeMap
    void SimulatorTimer::init(const TimeMap& timeMap, size_t report_step)
    {
        total_time_ = timeMap.getTotalTime();
        timesteps_.resize(timeMap.numTimesteps());
        for ( size_t i = 0; i < timeMap.numTimesteps(); ++i ) {
            timesteps_[i] = timeMap.getTimeStepLength(i);
        }

        setCurrentStepNum(report_step);
        start_date_ = boost::posix_time::from_time_t( timeMap.getStartTime(0)).date();
    }

    /// Whether the current step is the first step.
    bool SimulatorTimer::initialStep() const
    {
        return (current_step_ == 0);
    }

    /// Total number of steps.
    int SimulatorTimer::numSteps() const
    {
        return timesteps_.size();
    }

    /// Current step number.
    int SimulatorTimer::currentStepNum() const
    {
        return current_step_;
    }

    /// Set current step number.
    void SimulatorTimer::setCurrentStepNum(int step)
    {
        current_step_ = step;
        current_time_ = std::accumulate(timesteps_.begin(), timesteps_.begin() + step, 0.0);
    }

    /// Current step length.
    double SimulatorTimer::currentStepLength() const
    {
        assert(!done());
        return timesteps_[current_step_];
    }

    double SimulatorTimer::stepLengthTaken() const
    {
        assert(current_step_ > 0);
        return timesteps_[current_step_ - 1];
    }

    /// time elapsed since the start of the simulation [s].
    double SimulatorTimer::simulationTimeElapsed() const
    {
        return current_time_;
    }

    boost::posix_time::ptime SimulatorTimer::startDateTime() const
    {
        return boost::posix_time::ptime(start_date_);
    }

    boost::posix_time::ptime SimulatorTimer::currentDateTime() const
    {
        // Boost uses only 32 bit long for seconds, but 64 bit for milliseconds.
        // As a workaround for very large times we just use milliseconds.
        // The cast is necessary because boost::posix_time::milliseconds requires
        // an integer argument.
        return startDateTime() + boost::posix_time::milliseconds(static_cast<long long>(simulationTimeElapsed() / Ewoms::prefix::milli));
    }

    /// Total time.
    double SimulatorTimer::totalTime() const
    {
        return total_time_;
    }

    /// Set total time.
    /// This is primarily intended for multi-epoch schedules,
    /// where a timer for a given epoch does not have
    /// access to later timesteps.
    void SimulatorTimer::setTotalTime(double time)
    {
        total_time_ = time;
    }

    /// Print a report with current and total time etc.
    void SimulatorTimer::report(std::ostream& os) const
    {
        os << "\n\n---------------    Simulation step number " << currentStepNum() << "    ---------------"
           << "\n      Current time (days)     " << Ewoms::unit::convert::to(simulationTimeElapsed(), Ewoms::unit::day)
           << "\n      Current stepsize (days) " << Ewoms::unit::convert::to(currentStepLength(), Ewoms::unit::day)
           << "\n      Total time (days)       " << Ewoms::unit::convert::to(totalTime(), Ewoms::unit::day)
           << "\n" << std::endl;
    }

    /// Next step.
    SimulatorTimer& SimulatorTimer::operator++()
    {
        assert(!done());
        current_time_ += timesteps_[current_step_];
        ++current_step_;
        return *this;
    }

    /// Return true if op++() has been called numSteps() times.
    bool SimulatorTimer::done() const
    {
        return int(timesteps_.size()) == current_step_;
    }

    /// return copy of object
    std::unique_ptr< SimulatorTimerInterface >
    SimulatorTimer::clone() const
    {
       return std::make_unique<SimulatorTimer>(*this);
    }

} // namespace Ewoms
