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
#include<config.h>
#include<ewoms/eclsimulators/wells/parallelwellinfo.hh>
#include<vector>
#include<string>
#include<tuple>
#include<ostream>

#define BOOST_TEST_MODULE ParallelWellInfo
#include <boost/test/unit_test.hpp>
class MPIError {
public:
  /** @brief Constructor. */
  MPIError(std::string s, int e) : errorstring(s), errorcode(e){}
  /** @brief The error string. */
  std::string errorstring;
  /** @brief The mpi error code. */
  int errorcode;
};

#ifdef HAVE_MPI
void mpiErrorHandler(MPI_Comm *, int *errCode, ...){
  char *errString=new char[MPI_MAX_ERROR_STRING];
  int errLength;
  MPI_Error_string(*errCode, errString, &errLength);
  std::string s(errString, errLength);
  std::cerr << "An MPI Error ocurred:"<<std::endl<<s<<std::endl;
  delete[] errString;
  throw MPIError(s, *errCode);
}
#endif

struct MPIFixture
{
    MPIFixture()
    {
#if HAVE_MPI
    int m_argc = boost::unit_test::framework::master_test_suite().argc;
    char** m_argv = boost::unit_test::framework::master_test_suite().argv;
    helper = &Dune::MPIHelper::instance(m_argc, m_argv);
#ifdef MPI_2
    MPI_Comm_create_errhandler(mpiErrorHandler, &handler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler);
#else
        MPI_Errhandler_create(mpiErrorHandler, &handler);
        MPI_Errhandler_set(MPI_COMM_WORLD, handler);
#endif
#endif
    }
    ~MPIFixture()
    {
#if HAVE_MPI
        MPI_Finalize();
#endif
    }
    Dune::MPIHelper* helper;
#if HAVE_MPI
    MPI_Errhandler handler;
#endif
};

BOOST_GLOBAL_FIXTURE(MPIFixture);

// Needed for BOOST_CHECK_EQUAL_COLLECTIONS
namespace std
{
std::ostream& operator<<(std::ostream& os, const std::pair<std::string, bool>& p)
{
    return os << "{" << p.first << " "<< p.second << "}";
}
}
namespace Ewoms
{
std::ostream& operator<<(std::ostream& os, const Ewoms::ParallelWellInfo& w)
{
    return os << "{" << w.name() << " "<< w.hasLocalCells() << " "<<
        w.isOwner() << "}";
}
}

BOOST_AUTO_TEST_CASE(ParallelWellComparison)
{
    int argc = 0;
    char** argv = nullptr;
    const auto& helper = Dune::MPIHelper::instance(argc, argv);
    std::vector<std::pair<std::string,bool>> pairs;
    if (helper.rank() == 0)
        pairs = {{"Test1", true},{"Test2", true}, {"Test1", false} };
    else
        pairs = {{"Test1", false},{"Test2", true}, {"Test1", true} };

    std::vector<Ewoms::ParallelWellInfo> well_info;
    well_info.assign(pairs.begin(), pairs.end());

    BOOST_CHECK_EQUAL_COLLECTIONS(pairs.begin(), pairs.end(),
                                  well_info.begin(), well_info.end());

    BOOST_CHECK_EQUAL_COLLECTIONS(well_info.begin(), well_info.end(),
                                  pairs.begin(), pairs.end());

    BOOST_CHECK(well_info[0] < pairs[1]);
    BOOST_CHECK(pairs[0] != well_info[1]);
    BOOST_CHECK(pairs[0] < well_info[1]);
    BOOST_CHECK(well_info[0] == pairs[0]);

    BOOST_CHECK(well_info[0] != well_info[1]);

    Ewoms::ParallelWellInfo well0, well1;

    BOOST_CHECK(well0 == well1);
#if HAVE_MPI
    BOOST_CHECK(well0.communication()==helper.getLocalCommunicator());
#endif
    Ewoms::ParallelWellInfo well2("Test", false);
    std::pair<std::string, bool> pwell={"Test", true};
    BOOST_CHECK(well2 < pwell);
    Ewoms::ParallelWellInfo well3("Test", true);
    BOOST_CHECK(! (well3 < pwell));
    pwell.second = false;
    BOOST_CHECK(! (well3 < pwell));

    if (helper.rank() == 0)
        BOOST_CHECK(well_info[0].communication().size()==1);

#if HAVE_MPI
    Ewoms::ParallelWellInfo::Communication comm{MPI_COMM_WORLD};

    BOOST_CHECK(well_info[1].communication().size() == comm.size());

    if (helper.rank() > 0)
    {
        BOOST_CHECK(well_info[2].communication().size() == comm.size()-1);
    }
#endif

}

BOOST_AUTO_TEST_CASE(CommunicateAboveBelowSelf)
{
    auto comm = Dune::MPIHelper::getLocalCommunicator();
    Ewoms::CommunicateAboveBelow commAboveBelow{ comm };
    for(std::size_t count=0; count < 2; ++count)
    {
        std::vector<int> eclIndex = {0, 1, 2, 3, 7 , 8, 10, 11};
        std::vector<double> current(eclIndex.size());
        std::transform(eclIndex.begin(), eclIndex.end(), current.begin(),
                       [](double v){ return 1+10.0*v;});
        commAboveBelow.beginReset();
        for (std::size_t i = 0; i < current.size(); ++i)
        {
            if (i==0)
                commAboveBelow.pushBackEclIndex(-1, eclIndex[i]);
            else
                commAboveBelow.pushBackEclIndex(eclIndex[i-1], eclIndex[i]);
        }
        commAboveBelow.endReset();
        auto above = commAboveBelow.communicateAbove(-10.0, current.data(), current.size());
        BOOST_CHECK(above[0]==-10.0);
        BOOST_CHECK(above.size() == current.size());
        auto a = above.begin()+1;
        std::for_each(current.begin(), current.begin() + (current.size()-1),
                      [&a](double v){ BOOST_CHECK(*(a++) == v);});
        auto below = commAboveBelow.communicateBelow(-10.0, current.data(), current.size());
        BOOST_CHECK(below.back() == -10.0);
        BOOST_CHECK(below.size() == current.size());
        auto b = below.begin();
        std::for_each(current.begin()+1, current.end(),
                      [&b](double v){ BOOST_CHECK(*(b++) == v);});
    }
}

BOOST_AUTO_TEST_CASE(CommunicateAboveBelowSelf1)
{
    auto comm = Dune::MPIHelper::getLocalCommunicator();
    Ewoms::CommunicateAboveBelow commAboveBelow{ comm };
    for(std::size_t count=0; count < 2; ++count)
    {
        std::vector<int> eclIndex = {0};
        std::vector<double> current(eclIndex.size());
        std::transform(eclIndex.begin(), eclIndex.end(), current.begin(),
                       [](double v){ return 1+10.0*v;});
        commAboveBelow.beginReset();
        for (std::size_t i = 0; i < current.size(); ++i)
        {
            if (i==0)
                commAboveBelow.pushBackEclIndex(-1, eclIndex[i]);
            else
                commAboveBelow.pushBackEclIndex(eclIndex[i-1], eclIndex[i]);
        }
        commAboveBelow.endReset();
        auto above = commAboveBelow.communicateAbove(-10.0, current.data(), current.size());
        BOOST_CHECK(above[0]==-10.0);
        BOOST_CHECK(above.size() == current.size());
        auto a = above.begin()+1;
        std::for_each(current.begin(), current.begin() + (current.size()-1),
                      [&a](double v){ BOOST_CHECK(*(a++) == v);});
        auto below = commAboveBelow.communicateBelow(-10.0, current.data(), current.size());
        BOOST_CHECK(below.back() == -10.0);
        BOOST_CHECK(below.size() == current.size());
        auto b = below.begin();
        std::for_each(current.begin()+1, current.end(),
                      [&b](double v){ BOOST_CHECK(*(b++) == v);});
    }
}

BOOST_AUTO_TEST_CASE(CommunicateAboveBelowParallel)
{
    using MPIComm = typename Dune::MPIHelper::MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
    using Communication = Dune::Communication<MPIComm>;
#else
    using Communication = Dune::CollectiveCommunication<MPIComm>;
#endif
    auto comm = Communication(Dune::MPIHelper::getCommunicator());

    Ewoms::CommunicateAboveBelow commAboveBelow{ comm };
    for(std::size_t count=0; count < 2; ++count)
    {
        std::vector<int> globalEclIndex = {0, 1, 2, 3, 7 , 8, 10, 11};
        auto oldSize = globalEclIndex.size();
        std::size_t globalSize = 3 * comm.size();
        auto lastIndex = globalEclIndex.back();
        globalEclIndex.resize(globalSize);
        if ( globalSize > oldSize)
        {
            ++lastIndex;
            for(auto entry = globalEclIndex.begin() + oldSize;
                entry != globalEclIndex.end(); ++entry, ++lastIndex)
            {
                *entry = lastIndex;
            }
        }

        std::vector<double> globalCurrent(globalEclIndex.size());
        std::transform(globalEclIndex.begin(), globalEclIndex.end(), globalCurrent.begin(),
                       [](double v){ return 1+10.0*v;});

        std::vector<double> current(3);

        commAboveBelow.beginReset();
        for (std::size_t i = 0; i < current.size(); ++i)
        {
            auto gi = comm.rank() + comm.size() * i;

            if (gi==0)
            {
                commAboveBelow.pushBackEclIndex(-1, globalEclIndex[gi]);
            }
            else
            {
                commAboveBelow.pushBackEclIndex(globalEclIndex[gi-1], globalEclIndex[gi]);
            }
            current[i] = globalCurrent[gi];
        }
        commAboveBelow.endReset();
        auto above = commAboveBelow.communicateAbove(-10.0, current.data(), current.size());
        if (comm.rank() == 0)
            BOOST_CHECK(above[0]==-10.0);

        BOOST_CHECK(above.size() == current.size());

        for (std::size_t i = 0; i < current.size(); ++i)
        {
            auto gi = comm.rank() + comm.size() * i;
            if (gi > 0)
            {
                BOOST_CHECK(above[i]==globalCurrent[gi-1]);
            }
        }
        auto below = commAboveBelow.communicateBelow(-10.0, current.data(), current.size());
        if (comm.rank() == comm.size() - 1)
            BOOST_CHECK(below.back() == -10.0);

        BOOST_CHECK(below.size() == current.size());

        for (std::size_t i = 0; i < current.size(); ++i)
        {
            auto gi = comm.rank() + comm.size() * i;
            if (gi < globalCurrent.size() - 1)
            {
                BOOST_CHECK(below[i]==globalCurrent[gi+1]);
            }
        }
    }
}
