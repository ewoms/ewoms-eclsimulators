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
#ifndef EWOMS_READDECK_HH
#define EWOMS_READDECK_HH

#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/opmlog/eclipseprtlog.hh>
#include <ewoms/eclio/opmlog/logutil.hh>

#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/summaryconfig/summaryconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/messagelimits.hh>

#include <memory>
#include <string>

namespace Ewoms
{
enum class FileOutputMode {
    //! \brief No output to files.
    OUTPUT_NONE = 0,
    //! \brief Output only to log files, no eclipse output.
    OUTPUT_LOG_ONLY = 1,
    //! \brief Output to all files.
    OUTPUT_ALL = 3
};

// Setup the OpmLog backends
FileOutputMode setupLogging(int mpi_rank_, const std::string& deck_filename, const std::string& cmdline_output_dir, const std::string& cmdline_output, bool output_cout_, const std::string& stdout_log_id);

void setupMessageLimiter(const Ewoms::MessageLimits msgLimits,  const std::string& stdout_log_id);

/// \brief Reads the deck and creates all necessary objects if needed
///
/// If pointers already contains objects then they are used otherwise they are created and can be used outside later.
void readDeck(int rank, std::string& deckFilename, std::unique_ptr<Ewoms::Deck>& deck, std::unique_ptr<Ewoms::EclipseState>& eclipseState,
              std::unique_ptr<Ewoms::Schedule>& schedule, std::unique_ptr<Ewoms::SummaryConfig>& summaryConfig,
              std::unique_ptr<ErrorGuard> errorGuard, std::unique_ptr<ParseContext> parseContext,
              bool initFromRestart, bool checkDeck);
} // end namespace Ewoms

#endif // EWOMS_READDECK_HH
