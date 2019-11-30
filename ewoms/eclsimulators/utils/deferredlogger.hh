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

#ifndef EWOMS_DEFERREDLOGGER_HH
#define EWOMS_DEFERREDLOGGER_HH

#include <ewoms/eclio/opmlog/opmlog.hh>

#include <string>
#include <vector>

namespace Ewoms
{
    /** This class implements a deferred logger:
     * 1) messages can be pushed back to a vector
     * 2) a call to logMessages adds the messages to OpmLog backends
     * */

    class DeferredLogger
    {
    public:

        struct Message
        {
            int64_t flag;
            std::string tag;
            std::string text;
        };

        void info(const std::string& tag, const std::string& message);
        void warning(const std::string& tag, const std::string& message);
        void error(const std::string& tag, const std::string& message);
        void problem(const std::string& tag, const std::string& message);
        void bug(const std::string& tag, const std::string& message);
        void debug(const std::string& tag, const std::string& message);
        void note(const std::string& tag, const std::string& message);

        void info(const std::string& message);
        void warning(const std::string& message);
        void error(const std::string& message);
        void problem(const std::string& message);
        void bug(const std::string& message);
        void debug(const std::string& message);
        void note(const std::string& message);

        /// Log all messages to the OpmLog backends,
        /// and clear the message container.
        void logMessages();

        /// Clear the message container without logging them.
        void clearMessages();

    private:
        std::vector<Message> messages_;
        friend Ewoms::DeferredLogger gatherDeferredLogger(const Ewoms::DeferredLogger& local_deferredlogger);
    };

} // namespace Ewoms

#endif // EWOMS_DEFERREDLOGGER_HH
