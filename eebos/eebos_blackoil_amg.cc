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
 * \brief An eebos variant for blackoil cases that uses an AMG-based linear solver.
 *
 * AMG based linear solvers have some theoretical advantages, in particular for parallel
 * simulations on large grids. That said, they require substantially more CPU time to set
 * up, so the performance is usually slower that with plain-old stabilized BiCG and we
 * thus do not use AMG by default.
 */
#include "config.h"

#include "eebos.hh"
#include "eclamgbackend.hh"

#include <ewoms/numerics/utils/start.hh>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EebosAmgTypeTag, INHERITS_FROM(EebosTypeTag));

// use the eebos optimized linear solver backend
SET_TAG_PROP(EebosAmgTypeTag, LinearSolverSplice, EclAmgLinearSolver);

// reducing the residual in the linear solver by a factor of 20 is enough for us
SET_SCALAR_PROP(EebosAmgTypeTag, LinearSolverTolerance, 0.05);

// use the AMG-based linear solver

END_PROPERTIES

namespace Ewoms {
namespace CO2DefaultTables {
#include <ewoms/material/components/co2tables.inc.cc>
}}

int main(int argc, char **argv)
{
    typedef TTAG(EebosAmgTypeTag) ProblemTypeTag;
    return Ewoms::start<ProblemTypeTag>(argc, argv);
}
