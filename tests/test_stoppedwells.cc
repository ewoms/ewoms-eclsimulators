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

#include <config.h>
#include <chrono>

#define BOOST_TEST_MODULE StoppedWellsTests

#include <boost/test/unit_test.hpp>

#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/schedule.hh>

using namespace Ewoms;

BOOST_AUTO_TEST_CASE(TestStoppedWells)
{
    const std::string filename = "wells_stopped.data";
    Ewoms::Parser parser;
    Ewoms::Deck deck(parser.parseFile(filename));
    Ewoms::EclipseState eclipseState(deck);
    const Schedule sched(deck, eclipseState);

    // Both wells are open in the first schedule step
    {
        auto wells = sched.getWells(0);
        BOOST_CHECK(wells[0].getStatus() == Ewoms::Well::Status::OPEN);
        BOOST_CHECK(wells[1].getStatus() == Ewoms::Well::Status::OPEN);
    }

    // The injector is stopped
    {
        auto wells = sched.getWells(1);
        BOOST_CHECK(wells[0].getStatus() == Ewoms::Well::Status::STOP);
        BOOST_CHECK(wells[1].getStatus() == Ewoms::Well::Status::OPEN);
    }
}
