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

#if HAVE_DUNE_FEM
#warning "eebos is not compatible with dune-fem. Disabling."
#undef HAVE_DUNE_FEM
#elif defined(HAVE_DUNE_FEM)
#undef HAVE_DUNE_FEM
#endif

#include "eclproblem.hh"

#include <ewoms/numerics/utils/start.hh>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EebosPlainTypeTag, INHERITS_FROM(BlackOilModel, EclBaseProblem));

END_PROPERTIES

namespace Ewoms {
namespace CO2DefaultTables {
#include <ewoms/material/components/co2tables.inc.cc>
}}

int main(int argc, char **argv)
{
    using ProblemTypeTag = TTAG(EebosPlainTypeTag);
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
