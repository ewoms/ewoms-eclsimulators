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
 * \brief The main function for the stand alone gas-oil variant of eebos.
 *
 * This only calls the eebosGasWaterMain() function.
 */
#include "config.h"

#include "eebos_gaswater.hh"

#include "eebos.hh"
#include "starteebos.hh"

BEGIN_PROPERTIES

NEW_TYPE_TAG(EebosGasWaterTypeTag, INHERITS_FROM(EebosTypeTag));

//! The indices indices which only enable oil and water
SET_PROP(EebosGasWaterTypeTag, Indices)
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    using FluidSystem = GET_PROP_TYPE(TTAG(EebosTypeTag), FluidSystem);

public:
    typedef Ewoms::BlackOilTwoPhaseIndices<GET_PROP_VALUE(TypeTag, EnableSolvent),
                                         GET_PROP_VALUE(TypeTag, EnableSsaSolvent),
                                         GET_PROP_VALUE(TypeTag, EnablePolymer),
                                         GET_PROP_VALUE(TypeTag, EnableEnergy),
                                         GET_PROP_VALUE(TypeTag, EnableFoam),
                                         GET_PROP_VALUE(TypeTag, EnableBrine),
                                         /*PVOffset=*/0,
                                         /*disabledCompIdx=*/FluidSystem::oilCompIdx> type;
};

END_PROPERTIES

namespace Ewoms {

void eebosGasWaterSetDeck(std::unique_ptr<Ewoms::Deck> deck,
                         std::unique_ptr<Ewoms::ParseContext> parseContext,
                         std::unique_ptr<Ewoms::ErrorGuard> errorGuard,
                         double externalSetupTime)
{
    using ProblemTypeTag = TTAG(EebosGasWaterTypeTag);
    using Vanguard = GET_PROP_TYPE(ProblemTypeTag, Vanguard);

    Vanguard::setExternalSetupTime(externalSetupTime);
    Vanguard::setExternalParseContext(std::move(parseContext));
    Vanguard::setExternalErrorGuard(std::move(errorGuard));
    Vanguard::setExternalDeck(std::move(deck));
}

int eebosGasWaterMain(int argc, char **argv)
{
    using ProblemTypeTag = TTAG(EebosGasWaterTypeTag);
    return Ewoms::startEebos<ProblemTypeTag>(argc, argv);
}

}
