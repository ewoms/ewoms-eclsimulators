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
/*!
 * \file
 *
 * \brief A general-purpose simulator for ECL decks using the black-oil model.
 */
#include "config.h"

#include "eebos.hh"

#include <ewoms/numerics/utils/start.hh>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EebosOilWaterTypeTag, INHERITS_FROM(EebosTypeTag));

//! The indices indices which only enable oil and water
SET_PROP(EebosOilWaterTypeTag, Indices)
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    typedef typename GET_PROP_TYPE(TTAG(EebosTypeTag), FluidSystem) FluidSystem;

public:
    typedef Ewoms::BlackOilTwoPhaseIndices<GET_PROP_VALUE(TypeTag, EnableSolvent),
                                         GET_PROP_VALUE(TypeTag, EnablePolymer),
                                         GET_PROP_VALUE(TypeTag, EnableEnergy),
                                         GET_PROP_VALUE(TypeTag, EnableFoam),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::gasCompIdx> type;
};

END_PROPERTIES

namespace Ewoms {

void ebosOilWaterSetDeck(Ewoms::Deck* deck,
                         Ewoms::ParseContext* parseContext,
                         Ewoms::ErrorGuard* errorGuard,
                         double externalSetupTime)
{
    typedef TTAG(EebosOilWaterTypeTag) ProblemTypeTag;
    typedef GET_PROP_TYPE(ProblemTypeTag, Vanguard) Vanguard;

    Vanguard::setExternalSetupTime(externalSetupTime);
    Vanguard::setExternalParseContext(parseContext);
    Vanguard::setExternalErrorGuard(errorGuard);
    Vanguard::setExternalDeck(deck);
}

int ebosOilWaterMain(int argc, char **argv)
{
    typedef TTAG(EebosOilWaterTypeTag) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}

}
