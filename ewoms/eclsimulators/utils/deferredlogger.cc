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

#include <ewoms/eclsimulators/utils/deferredlogger.hh>

namespace Ewoms
{

    void DeferredLogger::info(const std::string& tag, const std::string& message)
    {
        messages_.push_back({Log::MessageType::Info, tag, message});
    }
    void DeferredLogger::warning(const std::string& tag, const std::string& message)
    {
        messages_.push_back({Log::MessageType::Warning, tag, message});
    }
    void DeferredLogger::error(const std::string& tag, const std::string& message)
    {
        messages_.push_back({Log::MessageType::Error, tag, message});
    }
    void DeferredLogger::problem(const std::string& tag, const std::string& message)
    {
        messages_.push_back({Log::MessageType::Problem, tag, message});
    }
    void DeferredLogger::bug(const std::string& tag, const std::string& message)
    {
        messages_.push_back({Log::MessageType::Bug, tag, message});
    }
    void DeferredLogger::debug(const std::string& tag, const std::string& message)
    {
        messages_.push_back({Log::MessageType::Debug, tag, message});
    }
    void DeferredLogger::note(const std::string& tag, const std::string& message)
    {
        messages_.push_back({Log::MessageType::Note, tag, message});
    }

    void DeferredLogger::info(const std::string& message)
    {
        messages_.push_back({Log::MessageType::Info, "", message});
    }
    void DeferredLogger::warning(const std::string& message)
    {
        messages_.push_back({Log::MessageType::Warning, "", message});
    }
    void DeferredLogger::error(const std::string& message)
    {
        messages_.push_back({Log::MessageType::Error, "", message});
    }
    void DeferredLogger::problem(const std::string& message)
    {
        messages_.push_back({Log::MessageType::Problem, "", message});
    }
    void DeferredLogger::bug(const std::string& message)
    {
        messages_.push_back({Log::MessageType::Bug, "", message});
    }
    void DeferredLogger::debug(const std::string& message)
    {
        messages_.push_back({Log::MessageType::Debug, "", message});
    }
    void DeferredLogger::note(const std::string& message)
    {
        messages_.push_back({Log::MessageType::Note, "", message});
    }

    void DeferredLogger::logMessages()
    {
        for (const auto& m : messages_) {
            OpmLog::addTaggedMessage(m.flag, m.tag, m.text);
        }
        messages_.clear();
    }

    void DeferredLogger::clearMessages()
    {
        messages_.clear();
    }

} // namespace Ewoms
