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

#ifndef BDABRIDGE_HH
#define BDABRIDGE_HH

#include <config.h>
#include "dune/istl/solver.hh" // for struct InverseOperatorResult

#include "dune/istl/bcrsmatrix.hh"

#if HAVE_CUDA
#include <ewoms/eclsimulators/linalg/bda/cusparsesolverbackend.hh>
#endif

namespace Ewoms
{

typedef Dune::InverseOperatorResult InverseOperatorResult;

/// BdaBridge acts as interface between opm-simulators with the cusparseSolver
/// if CUDA was not found during CMake, function bodies of this class are empty
class BdaBridge
{
private:
#if HAVE_CUDA
    std::unique_ptr<cusparseSolverBackend> backend;
    bool use_gpu;
#endif

public:
    /// Construct a BdaBridge
    /// \param[in] use_gpu                    true iff the cusparseSolver is used, is passed via command-line: '--use-gpu=[true|false]'
    /// \param[in] linear_solver_verbosity    verbosity of cusparseSolver
    /// \param[in] maxit                      maximum number of iterations for cusparseSolver
    /// \param[in] tolerance                  required relative tolerance for cusparseSolver
    BdaBridge(bool use_gpu, int linear_solver_verbosity, int maxit, double tolerance);

    /// Solve linear system, A*x = b
    /// \param[in] mat     matrix A, should be of type Dune::BCRSMatrix
    /// \param[in] b       vector b, should be of type Dune::BlockVector
    /// \param[in] result  summary of solver result
    template <class BridgeMatrix, class BridgeVector>
    void solve_system(BridgeMatrix *mat, BridgeVector &b, InverseOperatorResult &result);

    /// Get the resulting x vector
    /// \param[inout] x    vector x, should be of type Dune::BlockVector
    template <class BridgeVector>
    void get_result(BridgeVector &x);

}; // end class BdaBridge

}

#endif
