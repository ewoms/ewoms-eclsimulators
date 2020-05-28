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

#ifndef EWOMS_CUSPARSESOLVER_BACKEND_HH
#define EWOMS_CUSPARSESOLVER_BACKEND_HH

#include "cublas_v2.h"
#include "cusparse_v2.h"

#include "ewoms/eclsimulators/linalg/bda/bdaresult.hh"
#include <ewoms/eclsimulators/linalg/bda/wellcontributions.hh>

namespace Ewoms
{

/// This class implements a cusparse-based ilu0-bicgstab solver on GPU
class cusparseSolverBackend{

private:

    int minit;
    int maxit;
    double tolerance;

    cublasHandle_t cublasHandle;
    cusparseHandle_t cusparseHandle;
    cudaStream_t stream;
    cusparseMatDescr_t descr_B, descr_M, descr_L, descr_U;
    bsrilu02Info_t info_M;
    bsrsv2Info_t info_L, info_U;
    // b: bsr matrix, m: preconditioner
    double *d_bVals, *d_mVals;
    int *d_bCols, *d_mCols;
    int *d_bRows, *d_mRows;
    double *d_x, *d_b, *d_r, *d_rw, *d_p;
    double *d_pw, *d_s, *d_t, *d_v;
    void *d_buffer;
    int N, Nb, nnz, nnzb;
    double *vals_contiguous;

    int block_size;

    bool initialized = false;
    bool analysis_done = false;

    // verbosity
    // 0: print nothing during solves, only when initializing
    // 1: print number of iterations and final norm
    // 2: also print norm each iteration
    // 3: also print timings of different backend functions

    int verbosity = 0;

    /// Solve linear system using ilu0-bicgstab
    /// \param[in] wellContribs   contains all WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void gpu_pbicgstab(WellContributions& wellContribs, BdaResult& res);

    /// Initialize GPU and allocate memory
    /// \param[in] N                number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] nnz              number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] dim              size of block
    void initialize(int N, int nnz, int dim);

    /// Clean memory
    void finalize();

    /// Copy linear system to GPU
    /// \param[in] vals        array of nonzeroes, each block is stored row-wise, contains nnz values
    /// \param[in] rows        array of rowPointers, contains N/dim+1 values
    /// \param[in] cols        array of columnIndices, contains nnz values
    /// \param[in] b           input vector, contains N values
    void copy_system_to_gpu(double *vals, int *rows, int *cols, double *b);

    // Update linear system on GPU, don't copy rowpointers and colindices, they stay the same
    /// \param[in] vals        array of nonzeroes, each block is stored row-wise, contains nnz values
    /// \param[in] rows        array of rowPointers, contains N/dim+1 values, only used if COPY_ROW_BY_ROW is true
    /// \param[in] b           input vector, contains N values
    void update_system_on_gpu(double *vals, int *rows, double *b);

    /// Reset preconditioner on GPU, ilu0-decomposition is done inplace by cusparse
    void reset_prec_on_gpu();

    /// Analyse sparsity pattern to extract parallelism
    /// \return true iff analysis was successful
    bool analyse_matrix();

    /// Perform ilu0-decomposition
    /// \return true iff decomposition was successful
    bool create_preconditioner();

    /// Solve linear system
    /// \param[in] wellContribs   contains all WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void solve_system(WellContributions& wellContribs, BdaResult &res);

public:

    enum class cusparseSolverStatus {
        CUSPARSE_SOLVER_SUCCESS,
        CUSPARSE_SOLVER_ANALYSIS_FAILED,
        CUSPARSE_SOLVER_CREATE_PRECONDITIONER_FAILED,
        CUSPARSE_SOLVER_UNKNOWN_ERROR
    };

    /// Construct a cusparseSolver
    /// \param[in] linear_solver_verbosity    verbosity of cusparseSolver
    /// \param[in] maxit                      maximum number of iterations for cusparseSolver
    /// \param[in] tolerance                  required relative tolerance for cusparseSolver
    cusparseSolverBackend(int linear_solver_verbosity, int maxit, double tolerance);

    /// Destroy a cusparseSolver, and free memory
    ~cusparseSolverBackend();

    /// Solve linear system, A*x = b, matrix A must be in blocked-CSR format
    /// \param[in] N              number of rows, divide by dim to get number of blockrows
    /// \param[in] nnz            number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] dim            size of block
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] rows           array of rowPointers, contains N/dim+1 values
    /// \param[in] cols           array of columnIndices, contains nnz values
    /// \param[in] b              input vector, contains N values
    /// \param[in] wellContribs   contains all WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    /// \return                   status code
    cusparseSolverStatus solve_system(int N, int nnz, int dim, double *vals, int *rows, int *cols, double *b, WellContributions& wellContribs, BdaResult &res);

    /// Post processing after linear solve, now only copies resulting x vector back
    /// \param[inout] x        resulting x vector, caller must guarantee that x points to a valid array
    void post_process(double *x);

}; // end class cusparseSolverBackend

}

#endif

