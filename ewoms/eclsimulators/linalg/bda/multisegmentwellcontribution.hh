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

#ifndef MULTISEGMENTWELLCONTRIBUTION_HH
#define MULTISEGMENTWELLCONTRIBUTION_HH

#if !HAVE_CUDA
#error "This header file can only be included if the CUDA libraries are available"
#endif

#include <config.h>

#include <vector>

#include <cuda_runtime.h>

namespace Ewoms
{

    /// This class serves to duplicate the functionality of the MultisegmentWell
    /// A MultisegmentWell uses C, D and B and performs y -= (C^T * (D^-1 * (B*x)))
    /// B and C are matrices, with M rows and N columns, where N is the size of the matrix. They contain blocks of MultisegmentWell::numEq by MultisegmentWell::numWellEq.
    /// D is a MxM matrix, the square blocks have size MultisegmentWell::numWellEq.
    /// B*x and D*B*x are a vector with M*numWellEq doubles
    /// C*D*B*x is a vector with N*numEq doubles.

    class MultisegmentWellContribution
    {

    private:
        unsigned int dim;                        // size of blockvectors in vectors x and y, equal to MultisegmentWell::numEq
        unsigned int dim_wells;                  // size of blocks in C, B and D, equal to MultisegmentWell::numWellEq
        unsigned int N;                          // number of rows in vectors x and y, N == dim*Nb
        unsigned int Nb;                         // number of blockrows in x and y
        unsigned int M;                          // number of rows, M == dim_wells*Mb
        unsigned int Mb;                         // number of blockrows in C, D and B

        cudaStream_t stream; // not actually used yet, will be when MultisegmentWellContribution are applied on GPU

        // C and B are stored in BCRS format, D is stored in CSC format (Dune::UMFPack)
        // Sparsity pattern for C is not stored, since it is the same as B
        unsigned int DnumBlocks;             // number of blocks in D
        unsigned int BnumBlocks;             // number of blocks in C and B
        std::vector<double> Cvals;
        std::vector<double> Dvals;
        std::vector<double> Bvals;
        std::vector<int> Dcols;              // Columnpointers, contains M+1 entries
        std::vector<unsigned int> Bcols;
        std::vector<int> Drows;              // Rowindicies, contains DnumBlocks*dim*dim_wells entries
        std::vector<unsigned int> Brows;
        std::vector<double> z1;          // z1 = B * x
        std::vector<double> z2;          // z2 = D^-1 * B * x
        void *UMFPACK_Symbolic, *UMFPACK_Numeric;

    public:

        /// Set a cudaStream to be used
        /// \param[in] stream           the cudaStream that is used
        void setCudaStream(cudaStream_t stream);

        /// Create a new MultisegmentWellContribution
        /// Matrices C and B are passed in Blocked CSR, matrix D in CSC
        /// The variables representing C, B and D will go out of scope when MultisegmentWell::addWellContribution() ends
        /// \param[in] dim              size of blocks in blockvectors x and y, equal to MultisegmentWell::numEq
        /// \param[in] dim_wells        size of blocks of C, B and D, equal to MultisegmentWell::numWellEq
        /// \param[in] Nb               number of blocks in vectors x and y
        /// \param[in] Mb               number of blockrows in C, B and D
        /// \param[in] BnumBlocks       number of blocks in C and B
        /// \param[in] Bvalues          nonzero values of matrix B
        /// \param[in] BcolIndices      columnindices of blocks of matrix B
        /// \param[in] BrowPointers     rowpointers of matrix B
        /// \param[in] DnumBlocks       number of blocks in D
        /// \param[in] Dvalues          nonzero values of matrix D
        /// \param[in] DcolPointers     columnpointers of matrix D
        /// \param[in] DrowIndices      rowindices of matrix D
        /// \param[in] Cvalues          nonzero values of matrix C
        MultisegmentWellContribution(unsigned int dim, unsigned int dim_wells,
            unsigned int Nb, unsigned int Mb,
            unsigned int BnumBlocks, std::vector<double> &Bvalues, std::vector<unsigned int> &BcolIndices, std::vector<unsigned int> &BrowPointers,
            unsigned int DnumBlocks, double *Dvalues, int *DcolPointers, int *DrowIndices,
            std::vector<double> &Cvalues);

        /// Destroy a MultisegmentWellContribution, and free memory
        ~MultisegmentWellContribution();

        /// Apply the MultisegmentWellContribution on CPU
        /// performs y -= (C^T * (D^-1 * (B*x))) for MultisegmentWell
        /// \param[in] h_x          vector x, must be on CPU
        /// \param[inout] h_y       vector y, must be on CPU
        void apply(double *h_x, double *h_y);

    };

} //namespace Ewoms

#endif
