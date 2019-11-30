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
#ifndef EFLOW_OILWATER_POLYMER_INJECTIVITY_HH
#define EFLOW_OILWATER_POLYMER_INJECTIVITY_HH

#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>

namespace Ewoms {
  // void eflowOilWaterPolymerInjectivitySetDeck(Deck& deck, EclipseState& eclState);
int eflowOilWaterPolymerInjectivityMain(int argc, char** argv, bool outputCout, bool outputFiles);
}

#endif // EFLOW_OILWATER_POLYMER_INJECTIVITY_HH
