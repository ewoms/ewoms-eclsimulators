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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#ifndef EWOMS_EBOS_NNCSORTER_HH
#define EWOMS_EBOS_NNCSORTER_HH

#include <ewoms/eclio/parser/eclipsestate/grid/nnc.hh>

#include <vector>

namespace Ewoms
{
/// \brief Scale NNC data wit informtion form EDITNNC and sort it.
/// \param nncData The NNC data as provided by the deck.
/// \param editnncData The EDITNNC data as provided by the deck.
/// \return A lexicographically sorted vector of the scaled NNC data.
///         For each entry entry.cell1<entry.cell2 will hold for convenience.
std::vector<Ewoms::NNCdata> sortNncAndApplyEditnnc(const std::vector<Ewoms::NNCdata>& nncData, std::vector<Ewoms::NNCdata> editnncData,
                                                 bool log = true);
}
#endif
