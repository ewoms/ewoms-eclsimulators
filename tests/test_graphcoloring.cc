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
#include <config.h>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/graph.hh>

#include <ewoms/eclsimulators/linalg/graphcoloring.hh>

#define BOOST_TEST_MODULE GraphColoringTest
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

///! \brief check that all indices are represented in the new ordering.
void checkAllIndices(const std::vector<std::size_t>& ordering)
{
    std::vector<int> counters(ordering.size(), 0);
    for(auto index: ordering)
    {
        ++counters[index];
    }

    for(auto count: counters)
    {
        BOOST_CHECK(count==1);
    }
}

BOOST_AUTO_TEST_CASE(TestWelschPowell)
{
    using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>>;
    using Graph = Dune::Amg::MatrixGraph<Matrix>;
    int N = 10;
    Matrix matrix(N*N, N*N, 5, 0.4, Matrix::implicit);
    for( int j = 0; j < N; j++)
    {
        for(int i = 0; i < N; i++)
        {
            auto index = j*10+i;
            matrix.entry(index,index) = 1;

            if ( i > 0 )
            {
                matrix.entry(index,index-1) = 1;
            }
            if ( i  < N - 1)
            {
                matrix.entry(index,index+1) = 1;
            }

            if ( j > 0 )
            {
                matrix.entry(index,index-N) = 1;
            }
            if ( j  < N - 1)
            {
                matrix.entry(index,index+N) = 1;
            }

        }
    }
    matrix.compress();

    Graph graph(matrix);
    auto colorsTuple = Ewoms::colorVerticesWelshPowell(graph);
    const auto& colors = std::get<0>(colorsTuple);
    const auto& verticesPerColor = std::get<2>(colorsTuple);
    auto noColors = std::get<1>(colorsTuple);
    auto firstCornerColor = colors[0];
    BOOST_CHECK(noColors == 2);

    // Check for checkerboard coloring

    for( int j = 0, index = 0; j < N; j++)
    {
        auto expectedColor = firstCornerColor;

        for(int i = 0; i < N; i++)
        {
            BOOST_CHECK(colors[index]==expectedColor);
            index++;
            expectedColor = (expectedColor + 1) % 2;
        }
        firstCornerColor=(firstCornerColor + 1) % 2;
    }
    auto newOrder = Ewoms::reorderVerticesPreserving(colors, noColors, verticesPerColor,
                                                   graph);
    std::vector<std::size_t> colorIndex(noColors, 0);
    std::partial_sum(verticesPerColor.begin(),
                    verticesPerColor.begin()+verticesPerColor.size()-1,
                    colorIndex.begin()+1);

    for (auto vertex : graph)
    {
        BOOST_CHECK(colorIndex[colors[vertex]]++ == newOrder[vertex]);
    }
    checkAllIndices(newOrder);
    newOrder = Ewoms::reorderVerticesSpheres(colors, noColors, verticesPerColor,
                                           graph, 0);
    checkAllIndices(newOrder);
}
