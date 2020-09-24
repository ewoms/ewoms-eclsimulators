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

#include <ewoms/eclsimulators/utils/parallelserialization.hh>

#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/dynamicstate.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/action/astnode.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/sicd.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/valve.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqactive.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqastnode.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/udq/udqconfig.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wlist.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wlistmanager.hh>
#include <ewoms/eclio/parser/eclipsestate/summaryconfig/summaryconfig.hh>

#include <eebos/eclmpiserializer.hh>

#include <dune/common/parallel/mpihelper.hh>

namespace Ewoms {

void eclStateBroadcast(EclipseState& eclState, Schedule& schedule,
                       SummaryConfig& summaryConfig)
{
    Ewoms::EclMpiSerializer ser(Dune::MPIHelper::getCollectiveCommunication());
    ser.broadcast(eclState);
    ser.broadcast(schedule);
    ser.broadcast(summaryConfig);
}

}
