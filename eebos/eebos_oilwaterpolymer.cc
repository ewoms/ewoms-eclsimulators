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

NEW_TYPE_TAG(EebosOilWaterPolymerTypeTag, INHERITS_FROM(EebosTypeTag));

// enable the polymer extension of the black oil model
SET_BOOL_PROP(EebosOilWaterPolymerTypeTag, EnablePolymer, true);

//! The indices indices which only enable oil and water
SET_PROP(EebosOilWaterPolymerTypeTag, Indices)
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using FluidSystem = GET_PROP_TYPE(TTAG(EebosTypeTag), FluidSystem);

public:
    typedef Ewoms::BlackOilTwoPhaseIndices<GET_PROP_VALUE(TypeTag, EnableSolvent),
                                         GET_PROP_VALUE(TypeTag, EnablePolymer),
                                         GET_PROP_VALUE(TypeTag, EnableEnergy),
                                         GET_PROP_VALUE(TypeTag, EnableFoam),
                                         GET_PROP_VALUE(TypeTag, EnableBrine),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::gasCompIdx> type;
};

END_PROPERTIES

namespace Ewoms {

void eebosOilWaterPolymerSetDeck(Ewoms::Deck* deck,
                                Ewoms::ParseContext* parseContext,
                                Ewoms::ErrorGuard* errorGuard,
                                double externalSetupTime)
{
    using ProblemTypeTag = TTAG(EebosOilWaterPolymerTypeTag);
    using Vanguard = GET_PROP_TYPE(ProblemTypeTag, Vanguard);

    Vanguard::setExternalSetupTime(externalSetupTime);
    Vanguard::setExternalParseContext(parseContext);
    Vanguard::setExternalErrorGuard(errorGuard);
    Vanguard::setExternalDeck(deck);
}

int eebosOilWaterPolymerMain(int argc, char **argv)
{
    typedef TTAG(EebosOilWaterPolymerTypeTag) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}

}
