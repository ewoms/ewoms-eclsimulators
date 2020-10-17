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

#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE ParallelIstlInformation
#include <boost/test/unit_test.hpp>
#include "duneistltesthelpers.hh"
#include <ewoms/eclsimulators/linalg/parallelistlinformation.hh>
#include <functional>
#ifdef HAVE_DUNE_ISTL

template<typename T>
void runSumMaxMinTest(const T offset)
{
    const int N=100;
    int start, end, istart, iend;
    std::tie(start,istart,iend,end) = computeRegions(N);
    Ewoms::ParallelISTLInformation comm(MPI_COMM_WORLD);
    auto mat = create1DLaplacian(*comm.indexSet(), N, start, end, istart, iend);
    std::vector<T> x(end-start);
    assert(comm.indexSet()->size()==x.size());
    for(auto it=comm.indexSet()->begin(), itend=comm.indexSet()->end(); it!=itend; ++it)
        x[it->local()]=it->global()+offset;
    auto containers = std::make_tuple(x, x, x, x, x);
    auto operators  = std::make_tuple(Ewoms::Reduction::makeGlobalSumFunctor<T>(),
                                      Ewoms::Reduction::makeGlobalMaxFunctor<T>(),
                                      Ewoms::Reduction::makeGlobalMinFunctor<T>(),
                                      Ewoms::Reduction::makeInnerProductFunctor<T>(),
                                      Ewoms::Reduction::makeLInfinityNormFunctor<T>());
    auto values     = std::tuple<T,T,T,T,T>(0,0,100000, 0, 0);
    auto oldvalues  = values;
    start = offset;
    end   = start+N;
    comm.computeReduction(containers,operators,values);
    BOOST_CHECK(std::get<0>(values)==std::get<0>(oldvalues)+((N-1+2*offset)*N)/2);
    BOOST_CHECK(std::get<1>(values)==std::max(N+offset-1, std::get<1>(oldvalues)));
    BOOST_CHECK(std::get<2>(values)==std::min(offset, std::get<2>(oldvalues)));
    BOOST_CHECK(std::get<3>(values)==((end-1)*end*(2*end-1)-(start-1)*start*(2*start-1))/6+std::get<3>(oldvalues));
    // Must avoid std::abs() directly to prevent ambiguity with unsigned integers.
    Ewoms::Reduction::detail::MaxAbsFunctor<T> maxabsfunc;
    BOOST_CHECK(std::get<4>(values)==maxabsfunc(offset, N+offset-1));
}

BOOST_AUTO_TEST_CASE(tupleReductionTestInt)
{
    runSumMaxMinTest<int>(-200);
    runSumMaxMinTest<int>(0);
    runSumMaxMinTest<int>(20);
    runSumMaxMinTest<int>(-20);
}

BOOST_AUTO_TEST_CASE(tupleReductionTestUnsignedInt)
{
    runSumMaxMinTest<std::size_t>(0);
    runSumMaxMinTest<std::size_t>(20);
}
BOOST_AUTO_TEST_CASE(tupleReductionTestFloat)
{
    runSumMaxMinTest<float>(-200);
    runSumMaxMinTest<float>(0);
    runSumMaxMinTest<float>(20);
    runSumMaxMinTest<float>(-20);
}

BOOST_AUTO_TEST_CASE(singleContainerReductionTest)
{
    int N=100;
    int start, end, istart, iend;
    std::tie(start,istart,iend,end) = computeRegions(N);
    Ewoms::ParallelISTLInformation comm(MPI_COMM_WORLD);
    auto mat = create1DLaplacian(*comm.indexSet(), N, start, end, istart, iend);
    std::vector<int> x(end-start);
    assert(comm.indexSet()->size()==x.size());
    for(auto it=comm.indexSet()->begin(), itend=comm.indexSet()->end(); it!=itend; ++it)
        x[it->local()]=it->global();
    int value = 1;
    int oldvalue = value;
    comm.computeReduction(x,Ewoms::Reduction::makeGlobalSumFunctor<int>(),value);
    BOOST_CHECK(value==oldvalue+((N-1)*N)/2);
}
#endif
