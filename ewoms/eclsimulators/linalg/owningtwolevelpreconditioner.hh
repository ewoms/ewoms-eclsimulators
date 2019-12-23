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

#ifndef EWOMS_OWNINGTWOLEVELPRECONDITIONER_HH
#define EWOMS_OWNINGTWOLEVELPRECONDITIONER_HH

#include <ewoms/eclsimulators/linalg/preconditionerwithupdate.hh>
#include <ewoms/eclsimulators/linalg/pressuresolverpolicy.hh>
#include <ewoms/eclsimulators/linalg/pressuretransferpolicy.hh>
#include <ewoms/eclsimulators/linalg/getquasiimpesweights.hh>
#include <ewoms/eclsimulators/linalg/twolevelmethodcpr.hh>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/paamg/amg.hh>

#include <boost/property_tree/ptree.hpp>

#include <fstream>
#include <type_traits>

namespace Dune
{

// Circular dependency between PreconditionerFactory [which can make an OwningTwoLevelPreconditioner]
// and OwningTwoLevelPreconditioner [which uses PreconditionerFactory to choose the fine-level smoother]
// must be broken, accomplished by forward-declaration here.
template <class Operator, class Comm = Dune::Amg::SequentialInformation>
class EclSimulatorsPreconditionerFactory;

// Must forward-declare FlexibleSolver as we want to use it as solver for the pressure system.
template <class MatrixTypeT, class VectorTypeT>
class FlexibleSolver;

/// A version of the two-level preconditioner that is:
/// - Self-contained, because it owns its policy components.
/// - Flexible, because it uses the runtime-flexible solver
///   and preconditioner factory.
template <class OperatorType,
          class VectorType,
          bool transpose = false,
          class Communication = Dune::Amg::SequentialInformation>
class OwningTwoLevelPreconditioner : public Dune::PreconditionerWithUpdate<VectorType, VectorType>
{
public:
    using pt = boost::property_tree::ptree;
    using MatrixType = typename OperatorType::matrix_type;
    using PrecFactory = EclSimulatorsPreconditionerFactory<OperatorType, Communication>;

    OwningTwoLevelPreconditioner(const OperatorType& linearoperator, const pt& prm)
        : linear_operator_(linearoperator)
        , finesmoother_(PrecFactory::create(linearoperator, prm.get_child("finesmoother")))
        , comm_(nullptr)
        , weights_(Ewoms::Amg::getQuasiImpesWeights<MatrixType, VectorType>(
              linearoperator.getmat(), prm.get<int>("pressure_var_index"), transpose))
        , levelTransferPolicy_(dummy_comm_, weights_, prm.get<int>("pressure_var_index"))
        , coarseSolverPolicy_(prm.get_child("coarsesolver"))
        , twolevel_method_(linearoperator,
                           finesmoother_,
                           levelTransferPolicy_,
                           coarseSolverPolicy_,
                           transpose ? 1 : 0,
                           transpose ? 0 : 1)
        , prm_(prm)
    {
        if (prm.get<int>("verbosity") > 10) {
            std::ofstream outfile(prm.get<std::string>("weights_filename"));
            if (!outfile) {
                throw std::runtime_error("Could not write weights");
            }
            Dune::writeMatrixMarket(weights_, outfile);
        }
    }

    OwningTwoLevelPreconditioner(const OperatorType& linearoperator, const pt& prm, const Communication& comm)
        : linear_operator_(linearoperator)
        , finesmoother_(PrecFactory::create(linearoperator, prm.get_child("finesmoother"), comm))
        , comm_(&comm)
        , weights_(Ewoms::Amg::getQuasiImpesWeights<MatrixType, VectorType>(
              linearoperator.getmat(), prm.get<int>("pressure_var_index"), transpose))
        , levelTransferPolicy_(*comm_, weights_, prm.get<int>("pressure_var_index"))
        , coarseSolverPolicy_(prm.get_child("coarsesolver"))
        , twolevel_method_(linearoperator,
                           finesmoother_,
                           levelTransferPolicy_,
                           coarseSolverPolicy_,
                           transpose ? 1 : 0,
                           transpose ? 0 : 1)
        , prm_(prm)
    {
        if (prm.get<int>("verbosity") > 10) {
            std::ofstream outfile(prm.get<std::string>("weights_filename"));
            if (!outfile) {
                throw std::runtime_error("Could not write weights");
            }
            Dune::writeMatrixMarket(weights_, outfile);
        }
    }

    virtual void pre(VectorType& x, VectorType& b) override
    {
        twolevel_method_.pre(x, b);
    }

    virtual void apply(VectorType& v, const VectorType& d) override
    {
        twolevel_method_.apply(v, d);
    }

    virtual void post(VectorType& x) override
    {
        twolevel_method_.post(x);
    }

    virtual void update() override
    {
        Ewoms::Amg::getQuasiImpesWeights<MatrixType, VectorType>(
            linear_operator_.getmat(), prm_.get<int>("pressure_var_index"), transpose, weights_);
        updateImpl(comm_);
    }

    virtual Dune::SolverCategory::Category category() const override
    {
        return linear_operator_.category();
    }

private:
    using PressureMatrixType = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
    using PressureVectorType = Dune::BlockVector<Dune::FieldVector<double, 1>>;
    using SeqCoarseOperatorType = Dune::MatrixAdapter<PressureMatrixType, PressureVectorType, PressureVectorType>;
    using ParCoarseOperatorType
        = Dune::OverlappingSchwarzOperator<PressureMatrixType, PressureVectorType, PressureVectorType, Communication>;
    using CoarseOperatorType = std::conditional_t<std::is_same<Communication, Dune::Amg::SequentialInformation>::value,
                                                  SeqCoarseOperatorType,
                                                  ParCoarseOperatorType>;
    using LevelTransferPolicy = Ewoms::PressureTransferPolicy<OperatorType, CoarseOperatorType, Communication, transpose>;
    using CoarseSolverPolicy = Dune::Amg::PressureSolverPolicy<CoarseOperatorType,
                                                               FlexibleSolver<PressureMatrixType, PressureVectorType>,
                                                               LevelTransferPolicy>;
    using TwoLevelMethod
        = Dune::Amg::TwoLevelMethodCpr<OperatorType, CoarseSolverPolicy, Dune::Preconditioner<VectorType, VectorType>>;

    // Handling parallel vs serial instantiation of preconditioner factory.
    template <class Comm>
    void updateImpl(const Comm*)
    {
        // Parallel case.
        finesmoother_ = PrecFactory::create(linear_operator_, prm_.get_child("finesmoother"), *comm_);
        twolevel_method_.updatePreconditioner(finesmoother_, coarseSolverPolicy_);
    }

    void updateImpl(const Dune::Amg::SequentialInformation*)
    {
        // Serial case.
        finesmoother_ = PrecFactory::create(linear_operator_, prm_.get_child("finesmoother"));
        twolevel_method_.updatePreconditioner(finesmoother_, coarseSolverPolicy_);
    }

    const OperatorType& linear_operator_;
    std::shared_ptr<Dune::Preconditioner<VectorType, VectorType>> finesmoother_;
    const Communication* comm_;
    VectorType weights_;
    LevelTransferPolicy levelTransferPolicy_;
    CoarseSolverPolicy coarseSolverPolicy_;
    TwoLevelMethod twolevel_method_;
    boost::property_tree::ptree prm_;
    Communication dummy_comm_;
};

} // namespace Dune

#endif // EWOMS_OWNINGTWOLEVELPRECONDITIONER_HH
