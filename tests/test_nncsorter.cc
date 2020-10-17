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

#define BOOST_TEST_MODULE NNCSortTest
#include <boost/test/unit_test.hpp>
#include <eebos/nncsorter.hh>
#include <ewoms/eclio/parser/eclipsestate/grid/nnc.hh>

#include <cmath>

BOOST_AUTO_TEST_CASE(Test1) {
    std::vector<Ewoms::NNCdata> nncDataIn =
        { {9, 8, 10.0 }, { 1, 2, 3.0 }, { 3, 4, 2.0 }, { 2, 1, 5.0 } };
    std::vector<Ewoms::NNCdata> editnncData =
        { {20, 5, .001}, { 2, 1, .1}, {3, 4, .001}, {0, 0, 0.0}, {4, 3, 2.0} };
    std::vector<Ewoms::NNCdata> nncDataOut1 =
        { { 1, 2, 0.3 }, { 1, 2, 0.5 }, { 3, 4, 0.004 }, { 8, 9, 10.0 } };
    std::vector<Ewoms::NNCdata> nncDataOut2 =
        { { 1, 2, 0.5 }, { 1, 2, 0.3 }, { 3, 4, 0.4 }, { 8, 9, 10.0 } };

    auto nncDataProcessed = Ewoms::sortNncAndApplyEditnnc(nncDataIn, editnncData);
    BOOST_CHECK(nncDataProcessed.size() == nncDataOut1.size());
    auto expectedNnc1 = nncDataOut1.begin();
    auto expectedNnc2 = nncDataOut2.begin();

    for(const auto& entry: nncDataProcessed) {
        BOOST_CHECK( entry.cell1 == expectedNnc1->cell1);
        BOOST_CHECK( entry.cell2 == expectedNnc1->cell2);
        BOOST_CHECK( std::abs(entry.trans - expectedNnc1->trans) / std::abs(entry.trans) < 1e-5 ||
                std::abs(entry.trans - expectedNnc2->trans) / std::abs(entry.trans) < 1e-5);
        ++expectedNnc1;
        ++expectedNnc2;
    }
}
