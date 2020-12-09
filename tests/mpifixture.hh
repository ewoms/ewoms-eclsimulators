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
#ifndef EWOMS_SIMULATOR_TEST_MPIFIXTURE_HH
#define EWOMS_SIMULATOR_TEST_MPIFIXTURE_HH

#include <dune/common/parallel/mpihelper.hh>
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
#endif
