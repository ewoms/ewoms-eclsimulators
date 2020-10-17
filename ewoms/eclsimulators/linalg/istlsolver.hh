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
#ifndef EWOMS_ISTLSOLVER_EEBOS_HH
#define EWOMS_ISTLSOLVER_EEBOS_HH

#include <ewoms/common/parametersystem.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/eclsimulators/linalg/extractparallelgridinformationtoistl.hh>
#include <ewoms/eclsimulators/linalg/flexiblesolver.hh>
#include <ewoms/eclsimulators/linalg/parallelistlinformation.hh>
#include <ewoms/eclsimulators/linalg/welloperators.hh>
#include <ewoms/eclsimulators/linalg/writesystemmatrixhelper.hh>
#include <ewoms/eclsimulators/linalg/findoverlaprowsandcolumns.hh>
#include <ewoms/eclsimulators/linalg/getquasiimpesweights.hh>
#include <ewoms/eclsimulators/linalg/setuppropertytree.hh>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EFlowIstlSolver, INHERITS_FROM(EFlowIstlSolverParams));

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(GlobalEqVector);
NEW_PROP_TAG(SparseMatrixAdapter);
NEW_PROP_TAG(Indices);
NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(EclWellModel);
NEW_PROP_TAG(OwnerCellsFirst);

template <class TypeTag, class MyTypeTag>
struct EclWellModel;

//! Set the type of a global jacobian matrix for linear solvers that are based on
//! dune-istl.
SET_PROP(EFlowIstlSolver, SparseMatrixAdapter)
{
private:
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };
    typedef Dune::FieldMatrix<Scalar, numEq, numEq> Block;

public:
    typedef typename Ewoms::Linear::IstlSparseMatrixAdapter<Block> type;
};

END_PROPERTIES

namespace Ewoms
{

    /// This class solves the fully implicit black-oil system by
    /// solving the reduced system (after eliminating well variables)
    /// as a block-structured matrix (one block for all cell variables) for a fixed
    /// number of cell variables np .
    template <class TypeTag>
    class ISTLSolver
    {
    protected:
        using GridView = GET_PROP_TYPE(TypeTag, GridView);
        using Scalar = GET_PROP_TYPE(TypeTag, Scalar);
        using SparseMatrixAdapter = GET_PROP_TYPE(TypeTag, SparseMatrixAdapter);
        using Vector = GET_PROP_TYPE(TypeTag, GlobalEqVector);
        using Indices = GET_PROP_TYPE(TypeTag, Indices);
        using WellModel = GET_PROP_TYPE(TypeTag, EclWellModel);
        using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
        using Matrix = typename SparseMatrixAdapter::IstlMatrix;
        using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
        using FlexibleSolverType = Dune::FlexibleSolver<Matrix, Vector>;
        using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
        using WellModelOperator = WellModelAsLinearOperator<WellModel, Vector, Vector>;

#if HAVE_MPI
        using CommunicationType = Dune::OwnerOverlapCopyCommunication<int,int>;
#else
        using CommunicationType = Dune::CollectiveCommunication< int >;
#endif

    public:
        using AssembledLinearOperatorType = Dune::AssembledLinearOperator< Matrix, Vector, Vector >;

        static void registerParameters()
        {
            EFlowLinearSolverParameters::registerParameters<TypeTag>();
        }

        /// Construct a system solver.
        /// \param[in] parallelInformation In the case of a parallel run
        ///                                with dune-istl the information about the parallelization.
        explicit ISTLSolver(const Simulator& simulator)
            : simulator_(simulator),
              iterations_( 0 ),
              converged_(false),
              matrix_()
        {
            const bool on_io_rank = (simulator.gridView().comm().rank() == 0);
#if HAVE_MPI
            comm_.reset( new CommunicationType( simulator_.vanguard().grid().comm() ) );
#endif
            parameters_.template init<TypeTag>();
            prm_ = setupPropertyTree<TypeTag>(parameters_);
            const std::string gpu_mode = EWOMS_GET_PARAM(TypeTag, std::string, GpuMode);
            if (gpu_mode != "none") {
                EWOMS_THROW(std::logic_error,"So far, eflow does not support GPU based linear solvers!");
            }
            extractParallelGridInformationToISTL(simulator_.vanguard().grid(), parallelInformation_);

            // For some reason simulator_.model().elementMapper() is not initialized at this stage
            // Hence const auto& elemMapper = simulator_.model().elementMapper(); does not work.
            // Set it up manually
            using ElementMapper =
                Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
            ElementMapper elemMapper(simulator_.vanguard().gridView(), Dune::mcmgElementLayout());
            detail::findOverlapAndInterior(simulator_.vanguard().grid(), elemMapper, overlapRows_, interiorRows_);

            useWellConn_ = EWOMS_GET_PARAM(TypeTag, bool, MatrixAddWellContributions);
            const bool ownersFirst = EWOMS_GET_PARAM(TypeTag, bool, OwnerCellsFirst);
            if (!ownersFirst) {
                const std::string msg = "The linear solver no longer supports --owner-cells-first=false.";
                if (on_io_rank) {
                    OpmLog::error(msg);
                }
                EWOMS_THROW_NOLOG(std::runtime_error, msg);
            }

            interiorCellNum_ = detail::numMatrixRowsToUseInSolver(simulator_.vanguard().grid(), true);

            // Print parameters to PRT/DBG logs.
            if (on_io_rank) {
                std::ostringstream os;
                os << "Property tree for linear solver:\n";
                boost::property_tree::write_json(os, prm_, true);
                OpmLog::note(os.str());
            }
        }

        // nothing to clean here
        void eraseMatrix() {
        }

        void prepare(const SparseMatrixAdapter& M, Vector& b)
        {
            static bool firstcall = true;
#if HAVE_MPI
            if (firstcall && parallelInformation_.type() == typeid(ParallelISTLInformation)) {
                // Parallel case.
                const ParallelISTLInformation* parinfo = Ewoms::any_cast<ParallelISTLInformation>(&parallelInformation_);
                assert(parinfo);
                const size_t size = M.istlMatrix().N();
                parinfo->copyValuesTo(comm_->indexSet(), comm_->remoteIndices(), size, 1);
            }
#endif

            // update matrix entries for solvers.
            if (firstcall) {
                // eebos will not change the matrix object. Hence simply store a pointer
                // to the original one with a deleter that does nothing.
                // Outch! We need to be able to scale the linear system! Hence const_cast
                matrix_ = const_cast<Matrix*>(&M.istlMatrix());
            } else {
                // Pointers should not change
                if ( &(M.istlMatrix()) != matrix_ ) {
                        EWOMS_THROW(std::logic_error, "Matrix objects are expected to be reused when reassembling!"
                                  <<" old pointer was " << matrix_ << ", new one is " << (&M.istlMatrix()) );
                }
            }
            rhs_ = &b;

            if (isParallel() && prm_.get<std::string>("preconditioner.type") != "ParOverILU0") {
                makeOverlapRowsInvalid(getMatrix());
            }
            prepareFlexibleSolver();
            firstcall = false;
        }

        void setResidual(Vector& /* b */) {
            // rhs_ = &b; // Must be handled in prepare() instead.
        }

        void getResidual(Vector& b) const {
            b = *rhs_;
        }

        void setMatrix(const SparseMatrixAdapter& /* M */) {
            // matrix_ = &M.istlMatrix(); // Must be handled in prepare() instead.
        }

        bool solve(Vector& x) {
            // Write linear system if asked for.
            const int verbosity = prm_.get<int>("verbosity", 0);
            const bool write_matrix = verbosity > 10;
            if (write_matrix) {
                Ewoms::Helper::writeSystem(simulator_, //simulator is only used to get names
                                         getMatrix(),
                                         *rhs_,
                                         comm_.get());
            }

            // Solve system.
            Dune::InverseOperatorResult result;

            // use flexible istl solver.
            flexibleSolver_->apply(x, *rhs_, result);

            // Check convergence, iterations etc.
            checkConvergence(result);

            return converged_;
        }


        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x

        /// \copydoc NewtonIterationBlackoilInterface::iterations
        int iterations () const { return iterations_; }

        /// \copydoc NewtonIterationBlackoilInterface::parallelInformation
        const Ewoms::any& parallelInformation() const { return parallelInformation_; }

    protected:

        // 3x3 matrix block inversion was unstable at least 2.3 until and including
        // 2.5.0. There may still be some issue with the 4x4 matrix block inversion
        // we therefore still use the custom block inversion of eWoms
        typedef ParallelOverlappingILU0<Dune::BCRSMatrix<Dune::FieldMatrix<typename Matrix::field_type,
                                                                           Matrix::block_type::rows,
                                                                           Matrix::block_type::cols> >,
                                        Vector, Vector> SeqPreconditioner;

#if HAVE_MPI
        typedef Dune::OwnerOverlapCopyCommunication<int, int> Comm;
        // 3x3 matrix block inversion was unstable from at least 2.3 until and
        // including 2.5.0
        typedef ParallelOverlappingILU0<Matrix,Vector,Vector,Comm> ParPreconditioner;
#endif

        void checkConvergence( const Dune::InverseOperatorResult& result ) const
        {
            // store number of iterations
            iterations_ = result.iterations;
            converged_ = result.converged;

            // Check for failure of linear solver.
            if (!parameters_.ignoreConvergenceFailure_ && !result.converged) {
                const std::string msg("Convergence failure for linear solver.");
                EWOMS_THROW_NOLOG(NumericalIssue, msg);
            }
        }
    protected:

        bool isParallel() const {
#if HAVE_MPI
            return comm_->communicator().size() > 1;
#else
            return false;
#endif
        }

        void prepareFlexibleSolver()
        {

            std::function<Vector()> weightsCalculator = getWeightsCalculator();

            if (shouldCreateSolver()) {
                if (isParallel()) {
#if HAVE_MPI
                    if (useWellConn_) {
                        using ParOperatorType = Dune::OverlappingSchwarzOperator<Matrix, Vector, Vector, Comm>;
                        linearOperatorForFlexibleSolver_ = std::make_unique<ParOperatorType>(getMatrix(), *comm_);
                        flexibleSolver_ = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_, *comm_, prm_, weightsCalculator);
                    } else {
                        using ParOperatorType = WellModelGhostLastMatrixAdapter<Matrix, Vector, Vector, true>;
                        wellOperator_ = std::make_unique<WellModelOperator>(simulator_.problem().wellModel());
                        linearOperatorForFlexibleSolver_ = std::make_unique<ParOperatorType>(getMatrix(), *wellOperator_, interiorCellNum_);
                        flexibleSolver_ = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_, *comm_, prm_, weightsCalculator);
                    }
#endif
                } else {
                    if (useWellConn_) {
                        using SeqOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
                        linearOperatorForFlexibleSolver_ = std::make_unique<SeqOperatorType>(getMatrix());
                        flexibleSolver_ = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_, prm_, weightsCalculator);
                    } else {
                        using SeqOperatorType = WellModelMatrixAdapter<Matrix, Vector, Vector, false>;
                        wellOperator_ = std::make_unique<WellModelOperator>(simulator_.problem().wellModel());
                        linearOperatorForFlexibleSolver_ = std::make_unique<SeqOperatorType>(getMatrix(), *wellOperator_);
                        flexibleSolver_ = std::make_unique<FlexibleSolverType>(*linearOperatorForFlexibleSolver_, prm_, weightsCalculator);
                    }
                }
            }
            else
            {
                flexibleSolver_->preconditioner().update();
            }
        }

        /// Return true if we should (re)create the whole solver,
        /// instead of just calling update() on the preconditioner.
        bool shouldCreateSolver() const
        {
            // Decide if we should recreate the solver or just do
            // a minimal preconditioner update.
            if (!flexibleSolver_) {
                return true;
            }
            if (this->parameters_.cpr_reuse_setup_ == 0) {
                // Always recreate solver.
                return true;
            }
            if (this->parameters_.cpr_reuse_setup_ == 1) {
                // Recreate solver on the first iteration of every timestep.
                const int newton_iteration = this->simulator_.model().newtonMethod().numIterations();
                return newton_iteration == 0;
            }
            if (this->parameters_.cpr_reuse_setup_ == 2) {
                // Recreate solver if the last solve used more than 10 iterations.
                return this->iterations() > 10;
            }

            // Otherwise, do not recreate solver.
            assert(this->parameters_.cpr_reuse_setup_ == 3);

            return false;
        }

        /// Return an appropriate weight function if a cpr preconditioner is asked for.
        std::function<Vector()> getWeightsCalculator() const
        {
            std::function<Vector()> weightsCalculator;

            auto preconditionerType = prm_.get("preconditioner.type", "cpr");
            if (preconditionerType == "cpr" || preconditionerType == "cprt") {
                const bool transpose = preconditionerType == "cprt";
                const auto weightsType = prm_.get("preconditioner.weight_type", "quasiimpes");
                const auto pressureIndex = this->prm_.get("preconditioner.pressure_var_index", 1);
                if (weightsType == "quasiimpes") {
                    // weighs will be created as default in the solver
                    weightsCalculator = [this, transpose, pressureIndex]() {
                        return Ewoms::Amg::getQuasiImpesWeights<Matrix, Vector>(this->getMatrix(), pressureIndex, transpose);
                    };
                } else if (weightsType == "trueimpes") {
                    weightsCalculator = [this, pressureIndex]() {
                        return this->getTrueImpesWeights(pressureIndex);
                    };
                } else {
                    EWOMS_THROW(std::invalid_argument,
                              "Weights type " << weightsType << "not implemented for cpr."
                              << " Please use quasiimpes or trueimpes.");
                }
            }
            return weightsCalculator;
        }

        // Weights to make approximate pressure equations.
        // Calculated from the storage terms (only) of the
        // conservation equations, ignoring all other terms.
        Vector getTrueImpesWeights(int pressureVarIndex) const
        {
            Vector weights(rhs_->size());
            ElementContext elemCtx(simulator_);
            Ewoms::Amg::getTrueImpesWeights(pressureVarIndex, weights, simulator_.vanguard().gridView(),
                                          elemCtx, simulator_.model(),
                                          std::max(0, simulator_.taskletRunner().workerThreadIndex()));
            return weights;
        }

        /// Zero out off-diagonal blocks on rows corresponding to overlap cells
        /// Diagonal blocks on ovelap rows are set to diag(1.0).
        void makeOverlapRowsInvalid(Matrix& matrix) const
        {
            //value to set on diagonal
            const int numEq = Matrix::block_type::rows;
            typename Matrix::block_type diag_block(0.0);
            for (int eq = 0; eq < numEq; ++eq)
                diag_block[eq][eq] = 1.0;

            //loop over precalculated overlap rows and columns
            for (auto row = overlapRows_.begin(); row != overlapRows_.end(); row++ )
                {
                    int lcell = *row;
                    // Zero out row.
                    matrix[lcell] = 0.0;

                    //diagonal block set to diag(1.0).
                    matrix[lcell][lcell] = diag_block;
                }
        }

        Matrix& getMatrix()
        {
            return *matrix_;
        }

        const Matrix& getMatrix() const
        {
            return *matrix_;
        }

        const Simulator& simulator_;
        mutable int iterations_;
        mutable bool converged_;
        Ewoms::any parallelInformation_;

        // non-const to be able to scale the linear system
        Matrix* matrix_;
        Vector *rhs_;

        std::unique_ptr<FlexibleSolverType> flexibleSolver_;
        std::unique_ptr<AbstractOperatorType> linearOperatorForFlexibleSolver_;
        std::unique_ptr<WellModelAsLinearOperator<WellModel, Vector, Vector>> wellOperator_;
        std::vector<int> overlapRows_;
        std::vector<int> interiorRows_;
        std::vector<std::set<int>> wellConnectionsGraph_;

        bool useWellConn_;
        size_t interiorCellNum_;

        EFlowLinearSolverParameters parameters_;
        boost::property_tree::ptree prm_;
        bool scale_variables_;

        std::shared_ptr< CommunicationType > comm_;
    }; // end ISTLSolver

} // namespace Ewoms
#endif
