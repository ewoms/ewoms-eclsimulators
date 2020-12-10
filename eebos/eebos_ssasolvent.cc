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
 * \brief A general-purpose simulator for ECL decks using the black-oil model.
 */
#include "config.h"

#include "eebos.hh"

#include <ewoms/numerics/utils/start.hh>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EebosSsaSolventTypeTag, INHERITS_FROM(EebosTypeTag));

// enable the polymer extension of the black oil model
SET_BOOL_PROP(EebosSsaSolventTypeTag, EnableSsaSolvent, true);

END_PROPERTIES

namespace Ewoms {

void eebosSsaSolventSetDeck(Ewoms::Deck* deck,
                      Ewoms::ParseContext* parseContext,
                      Ewoms::ErrorGuard* errorGuard,
                      double externalSetupTime)
{
    using ProblemTypeTag = TTAG(EebosSsaSolventTypeTag);
    using Vanguard = GET_PROP_TYPE(ProblemTypeTag, Vanguard);

    Vanguard::setExternalSetupTime(externalSetupTime);
    Vanguard::setExternalParseContext(parseContext);
    Vanguard::setExternalErrorGuard(errorGuard);
    Vanguard::setExternalDeck(deck);
}

int eebosSsaSolventMain(int argc, char **argv)
{
    using ProblemTypeTag = TTAG(EebosSsaSolventTypeTag);
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}

}
