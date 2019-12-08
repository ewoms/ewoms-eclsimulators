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

#ifndef EWOMS_CONVERGENCEREPORT_HH
#define EWOMS_CONVERGENCEREPORT_HH

#include <cassert>
#include <numeric>
#include <string>
#include <vector>

namespace Ewoms
{

    /// Represents the convergence status of the whole simulator, to
    /// make it possible to query and store the reasons for
    /// convergence failures.
    class ConvergenceReport
    {
    public:

        // ----------- Types -----------

        enum Status { AllGood            = 0,
                      ReservoirFailed    = 1 << 0,
                      WellFailed         = 1 << 1 };
        enum struct Severity { None       = 0,
                               Normal     = 1,
                               TooLarge   = 2,
                               NotANumber = 3 };
        class ReservoirFailure
        {
        public:
            enum struct Type { Invalid, MassBalance, Cnv };
            ReservoirFailure(Type t, Severity s, int phase)
                : type_(t), severity_(s), phase_(phase)
            {
            }
            Type type() const { return type_; }
            Severity severity() const { return severity_; }
            int phase() const { return phase_; }
        private:
            Type type_;
            Severity severity_;
            int phase_;
        };
        class WellFailure
        {
        public:
            enum struct Type { Invalid, MassBalance, Pressure, ControlBHP, ControlTHP, ControlRate };
            WellFailure(Type t, Severity s, int phase, const std::string& well_name)
                : type_(t), severity_(s), phase_(phase), well_name_(well_name)
            {
            }
            Type type() const { return type_; }
            Severity severity() const { return severity_; }
            int phase() const { return phase_; }
            const std::string& wellName() const { return well_name_; }
        private:
            Type type_;
            Severity severity_;
            int phase_;
            std::string well_name_;
        };

        // ----------- Mutating member functions -----------

        ConvergenceReport()
            : status_{AllGood}
            , res_failures_{}
            , well_failures_{}
        {
        }

        void clear()
        {
            status_ = AllGood;
            res_failures_.clear();
            well_failures_.clear();
        }

        void setReservoirFailed(const ReservoirFailure& rf)
        {
            status_ = static_cast<Status>(status_ | ReservoirFailed);
            res_failures_.push_back(rf);
        }

        void setWellFailed(const WellFailure& wf)
        {
            status_ = static_cast<Status>(status_ | WellFailed);
            well_failures_.push_back(wf);
        }

        ConvergenceReport& operator+=(const ConvergenceReport& other)
        {
            status_ = static_cast<Status>(status_ | other.status_);
            res_failures_.insert(res_failures_.end(), other.res_failures_.begin(), other.res_failures_.end());
            well_failures_.insert(well_failures_.end(), other.well_failures_.begin(), other.well_failures_.end());
            assert(reservoirFailed() != res_failures_.empty());
            assert(wellFailed() != well_failures_.empty());
            return *this;
        }

        // ----------- Const member functions (queries) -----------

        bool converged() const
        {
            return status_ == AllGood;
        }

        bool reservoirFailed() const
        {
            return status_ & ReservoirFailed;
        }

        bool wellFailed() const
        {
            return status_ & WellFailed;
        }

        const std::vector<ReservoirFailure>& reservoirFailures() const
        {
            return res_failures_;
        }

        const std::vector<WellFailure>& wellFailures() const
        {
            return well_failures_;
        }

        Severity severityOfWorstFailure() const
        {
            // A function to get the worst of two severities.
            auto smax = [](Severity s1, Severity s2) {
                return s1 < s2 ? s2 : s1;
            };
            auto s = Severity::None;
            for (const auto f : res_failures_) {
                s = smax(s, f.severity());
            }
            for (const auto f : well_failures_) {
                s = smax(s, f.severity());
            }
            return s;
        }

    private:

        // ----------- Member variables -----------
        Status status_;
        std::vector<ReservoirFailure> res_failures_;
        std::vector<WellFailure> well_failures_;
    };

} // namespace Ewoms

#endif // EWOMS_CONVERGENCEREPORT_HH
