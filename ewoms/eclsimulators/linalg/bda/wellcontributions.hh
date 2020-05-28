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

#ifndef WELLCONTRIBUTIONS_HH
#define WELLCONTRIBUTIONS_HH

#include <config.h>

#if ! HAVE_CUDA
  #error "This file should only be included if CUDA is found"
#endif

#include <vector>

#include <ewoms/eclsimulators/linalg/bda/multisegmentwellcontribution.hh>

#include <cuda_runtime.h>

namespace Ewoms
{

    /// This class serves to eliminate the need to include the WellContributions into the matrix (with --matrix-add-well-contributions=true) for the cusparseSolver
    /// If the --matrix-add-well-contributions commandline parameter is true, this class should not be used
    /// So far, StandardWell and MultisegmentWell are supported
    /// A single instance (or pointer) of this class is passed to the cusparseSolver.
    /// For StandardWell, this class contains all the data and handles the computation. For MultisegmentWell, the vector 'multisegments' contains all the data. For more information, check the MultisegmentWellContribution class.

    /// A StandardWell uses C, D and B and performs y -= (C^T * (D^-1 * (B*x)))
    /// B and C are vectors, disguised as matrices and contain blocks of StandardWell::numEq by StandardWell::numStaticWellEq
    /// D is a block, disguised as matrix, the square block has size StandardWell::numStaticWellEq. D is actually stored as D^-1
    /// B*x and D*B*x are a vector with numStaticWellEq doubles
    /// C*D*B*x is a blocked matrix with a symmetric sparsity pattern, contains square blocks with size numEq. For every columnindex i, j in StandardWell::duneB_, there is a block on (i, j) in C*D*B*x.
    ///
    /// This class is used in 3 phases:
    /// - get total size of all wellcontributions that must be stored here
    /// - allocate memory
    /// - copy data of wellcontributions
    class WellContributions
    {

    private:
        unsigned int num_blocks = 0;             // total number of blocks in all wells
        unsigned int dim;                        // number of columns in blocks in B and C, equal to StandardWell::numEq
        unsigned int dim_wells;                  // number of rows in blocks in B and C, equal to StandardWell::numStaticWellEq
        unsigned int num_std_wells = 0;          // number of StandardWells in this object
        unsigned int num_ms_wells = 0;           // number of MultisegmentWells in this object, must equal multisegments.size()
        unsigned int num_blocks_so_far = 0;      // keep track of where next data is written
        unsigned int num_std_wells_so_far = 0;   // keep track of where next data is written
        unsigned int *val_pointers = nullptr;    // val_pointers[wellID] == index of first block for this well in Ccols and Bcols
        unsigned int N;                          // number of rows (not blockrows) in vectors x and y
        bool allocated = false;
        std::vector<MultisegmentWellContribution*> multisegments;
        cudaStream_t stream;

        // data for StandardWells, could remain nullptrs if not used
        double *d_Cnnzs = nullptr;
        double *d_Dnnzs = nullptr;
        double *d_Bnnzs = nullptr;
        int *d_Ccols = nullptr;
        int *d_Bcols = nullptr;
        double *d_z1 = nullptr;
        double *d_z2 = nullptr;
        unsigned int *d_val_pointers = nullptr;

        double *h_x = nullptr, *h_y = nullptr;  // CUDA pinned memory for GPU memcpy

    public:

        /// StandardWell has C, D and B matrices that need to be copied
        enum class MatrixType {
            C,
            D,
            B
        };

        /// Set a cudaStream to be used
        /// \param[in] stream           the cudaStream that is used to launch the kernel in
        void setCudaStream(cudaStream_t stream);

        /// Create a new WellContributions, implementation is empty
        WellContributions(){};

        /// Destroy a WellContributions, and free memory
        ~WellContributions();

        /// Apply all Wells in this object
        /// performs y -= (C^T * (D^-1 * (B*x))) for all Wells
        /// \param[in] d_x        vector x, must be on GPU
        /// \param[inout] d_y     vector y, must be on GPU
        void apply(double *d_x, double *d_y);

        /// Allocate memory for the StandardWells
        void alloc();

        /// Indicate how large the blocks of the StandardWell (C and B) are
        /// \param[in] dim         number of columns
        /// \param[in] dim_wells   number of rows
        void setBlockSize(unsigned int dim, unsigned int dim_wells);

        /// Indicate how large the next StandardWell is, this function cannot be called after alloc() is called
        /// \param[in] numBlocks   number of blocks in C and B of next StandardWell
        void addNumBlocks(unsigned int numBlocks);

        /// Store a matrix in this object, in blocked csr format, can only be called after alloc() is called
        /// \param[in] type        indicate if C, D or B is sent
        /// \param[in] colIndices  columnindices of blocks in C or B, ignored for D
        /// \param[in] values      array of nonzeroes
        /// \param[in] val_size    number of blocks in C or B, ignored for D
        void addMatrix(MatrixType type, int *colIndices, double *values, unsigned int val_size);

        /// Add a MultisegmentWellContribution, actually creates an object on heap that is destroyed in the destructor
        /// Matrices C and B are passed in Blocked CSR, matrix D in CSC
        /// \param[in] dim              size of blocks in vectors x and y, equal to MultisegmentWell::numEq
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
        void addMultisegmentWellContribution(unsigned int dim, unsigned int dim_wells,
            unsigned int Nb, unsigned int Mb,
            unsigned int BnumBlocks, std::vector<double> &Bvalues, std::vector<unsigned int> &BcolIndices, std::vector<unsigned int> &BrowPointers,
            unsigned int DnumBlocks, double *Dvalues, int *DcolPointers, int *DrowIndices,
            std::vector<double> &Cvalues);

        /// Return the number of wells added to this object
        /// \return the number of wells added to this object
        unsigned int getNumWells(){
            return num_std_wells + num_ms_wells;
        }
    };

} //namespace Ewoms

#endif
