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

namespace Ewoms {

bool eebosBlackOilDeckFileNameIsSet(int argc, char** argv)
{
    typedef TTAG(EebosTypeTag) ProblemTypeTag;

    // use the ewoms parameter machinery and the blackoil vanguard to handle the grunt of
    // the work
    EWOMS_RESET_PARAMS_(ProblemTypeTag);
    Ewoms::setupParameters_<ProblemTypeTag>(argc,
                                          const_cast<const char**>(argv),
                                          /*doRegistration=*/true,
                                          /*allowUnused=*/true,
                                          /*handleHelp=*/false);
    bool result = EWOMS_PARAM_IS_SET(ProblemTypeTag, std::string, EclDeckFileName);
    EWOMS_RESET_PARAMS_(ProblemTypeTag);

    return result;
}

std::string eebosBlackOilGetDeckFileName(int argc, char** argv)
{
    typedef TTAG(EebosTypeTag) ProblemTypeTag;
    typedef GET_PROP_TYPE(ProblemTypeTag, Vanguard) Vanguard;

    // use the ewoms parameter machinery and the blackoil vanguard to handle the grunt of
    // the work
    EWOMS_RESET_PARAMS_(ProblemTypeTag);
    Ewoms::setupParameters_<ProblemTypeTag>(argc,
                                          const_cast<const char**>(argv),
                                          /*doRegistration=*/true,
                                          /*allowUnused=*/true,
                                          /*handleHelp=*/false);
    std::string rawDeckFileName = EWOMS_GET_PARAM(ProblemTypeTag, std::string, EclDeckFileName);
    std::string result = Vanguard::canonicalDeckPath(rawDeckFileName).string();
    EWOMS_RESET_PARAMS_(ProblemTypeTag);

    return result;
}

std::unique_ptr<Ewoms::ParseContext> eebosBlackOilCreateParseContext(int argc, char** argv)
{
    typedef TTAG(EebosTypeTag) ProblemTypeTag;
    typedef GET_PROP_TYPE(ProblemTypeTag, Vanguard) Vanguard;

    // use the ewoms parameter machinery and the blackoil vanguard to handle the grunt of
    // the work
    EWOMS_RESET_PARAMS_(ProblemTypeTag);
    Ewoms::setupParameters_<ProblemTypeTag>(argc,
                                          const_cast<const char**>(argv),
                                          /*doRegistration=*/true,
                                          /*allowUnused=*/true,
                                          /*handleHelp=*/false);
    std::unique_ptr<Ewoms::ParseContext> result = Vanguard::createParseContext();
    EWOMS_RESET_PARAMS_(ProblemTypeTag);

    return result;
}

void eebosBlackOilSetDeck(Ewoms::Deck* deck,
                         Ewoms::ParseContext* parseContext,
                         Ewoms::ErrorGuard* errorGuard,
                         double externalSetupTime)
{
    typedef TTAG(EebosTypeTag) ProblemTypeTag;
    typedef GET_PROP_TYPE(ProblemTypeTag, Vanguard) Vanguard;

    Vanguard::setExternalSetupTime(externalSetupTime);
    Vanguard::setExternalParseContext(parseContext);
    Vanguard::setExternalErrorGuard(errorGuard);
    Vanguard::setExternalDeck(deck);
}

int eebosBlackOilMain(int argc, char **argv)
{
    typedef TTAG(EebosTypeTag) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}

}
