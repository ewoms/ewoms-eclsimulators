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
/*!
 * \file
 *
 * \brief The function prototypes required to start the solvent variant of eebos
 */
#ifndef EEBOS_SSA_SOLVENT_HH
#define EEBOS_SSA_SOLVENT_HH

#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/parser/parsecontext.hh>
#include <ewoms/eclio/parser/errorguard.hh>

namespace Ewoms {
void eebosSsaSolventSetDeck(Ewoms::Deck* deck,
                        Ewoms::ParseContext* parseContext,
                        Ewoms::ErrorGuard* errorGuard,
                        double externalSetupTime);

int eebosSsaSolventMain(int argc, char** argv);
}

#endif
