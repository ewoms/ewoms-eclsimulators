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
#include "config.h"

#define NVERBOSE  // Suppress own messages when throw()ing

#define BOOST_TEST_MODULE RelpermDiagnostics

#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION / 100000 == 1 && BOOST_VERSION / 100 % 1000 < 71
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/opmlog/counterlog.hh>

#include <ewoms/eclgrids/unstructuredgrid.h>
#include <ewoms/eclgrids/cart_grid.h>
#include <ewoms/eclgrids/gridmanager.hh>

#include <ewoms/eclsimulators/deprecated/props/satfunc/relpermdiagnostics.hh>
#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/deck/deck.hh>

BOOST_AUTO_TEST_SUITE ()

BOOST_AUTO_TEST_CASE(diagnosis)
{
    using namespace Ewoms;
    Parser parser;

    Ewoms::Deck deck = parser.parseFile("../tests/relpermDiagnostics.DATA");
    EclipseState eclState(deck);
    GridManager gm(eclState.getInputGrid());
    const UnstructuredGrid& grid = *gm.c_grid();
    std::shared_ptr<CounterLog> counterLog = std::make_shared<CounterLog>(Log::DefaultMessageTypes);
    OpmLog::addBackend( "COUNTERLOG" , counterLog );
    RelpermDiagnostics diagnostics;
    diagnostics.diagnosis(eclState, grid);
    BOOST_CHECK(counterLog->numMessages(Log::MessageType::Warning) > 1);
}
BOOST_AUTO_TEST_SUITE_END()
