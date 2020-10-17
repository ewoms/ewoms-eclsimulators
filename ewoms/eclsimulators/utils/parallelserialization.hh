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
#ifndef PARALLEL_SERIALIZATION_HH
#define PARALLEL_SERIALIZATION_HH

namespace Ewoms {

class EclipseState;
class Schedule;
class SummaryConfig;

/*! \brief Broadcasts an eclipse state from root node in parallel runs.
 *! \param eclState EclipseState to broadcast
 *! \param schedule Schedule to broadcast
 *! \param summaryConfig SummaryConfig to broadcast
*/
void eclStateBroadcast(EclipseState& eclState, Schedule& schedule,
                       SummaryConfig& summaryConfig);

} // end namespace Ewoms

#endif // PARALLEL_SERIALIZATION_HH
