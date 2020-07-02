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

#ifndef EWOMS_PRESSURE_SOLVER_POLICY_HH
#define EWOMS_PRESSURE_SOLVER_POLICY_HH

#include <ewoms/eclsimulators/linalg/pressuretransferpolicy.hh>

#include <boost/property_tree/ptree.hpp>

#include <dune/istl/solver.hh>
#include <dune/istl/owneroverlapcopy.hh>

namespace Dune
{
namespace Amg
{
    namespace pt = boost::property_tree;

    template <class OperatorType, class Solver, class LevelTransferPolicy>
    class PressureSolverPolicy
    {
    public:
        /** @brief The type of the linear operator used. */
        using Operator = OperatorType;
        /**
         * @brief Constructs the coarse solver policy.
         * @param prm Parameter tree specifying the solver details.
         */
        explicit PressureSolverPolicy(const pt::ptree prm)
            : prm_(prm)
        {
        }

    private:
        using X = typename Operator::range_type;
        /**
         * @brief A wrapper that makes an inverse operator out of AMG.
         *
         * The operator will use one step of AMG to approximately solve
         * the coarse level system.
         */
        struct PressureInverseOperator : public Dune::InverseOperator<X, X>
        {
#if HAVE_MPI
            template <typename GlobalIndex, typename LocalIndex>
            PressureInverseOperator(Operator& op,
                                    const boost::property_tree::ptree& prm,
                                    const Dune::OwnerOverlapCopyCommunication<GlobalIndex, LocalIndex>& comm)
                : linsolver_()
            {
                assert(op.category() == Dune::SolverCategory::overlapping);
                linsolver_ = std::make_unique<Solver>(op, comm, prm, std::function<X()>());
            }
#endif // HAVE_MPI

            PressureInverseOperator(Operator& op,
                                    const boost::property_tree::ptree& prm,
                                    const SequentialInformation&)
                : linsolver_()
            {
                assert(op.category() != Dune::SolverCategory::overlapping);
                linsolver_ = std::make_unique<Solver>(op, prm, std::function<X()>());
            }

            Dune::SolverCategory::Category category() const override
            {
                return linsolver_->category();
            }

            void apply(X& x, X& b, double reduction, Dune::InverseOperatorResult& res) override
            {
                linsolver_->apply(x, b, reduction, res);
            }

            void apply(X& x, X& b, Dune::InverseOperatorResult& res) override
            {
                linsolver_->apply(x, b, res);
            }

            void updatePreconditioner()
            {
                linsolver_->preconditioner().update();
            }

        private:
            std::unique_ptr<Solver> linsolver_;
        };

    public:
        /** @brief The type of solver constructed for the coarse level. */
        using CoarseLevelSolver = PressureInverseOperator;

        /**
         * @brief Constructs a coarse level solver.
         *
         * @param transferPolicy The policy describing the transfer between levels.
         * @return A pointer to the constructed coarse level solver.
         */
        template <class LTP>
        void setCoarseOperator(LTP& transferPolicy)
        {
            coarseOperator_ = transferPolicy.getCoarseLevelOperator();
        }
        template <class LTP>
        CoarseLevelSolver* createCoarseLevelSolver(LTP& transferPolicy)
        {
            coarseOperator_ = transferPolicy.getCoarseLevelOperator();
            auto& tp = dynamic_cast<LevelTransferPolicy&>(transferPolicy); // TODO: make this unnecessary.
            PressureInverseOperator* inv
                = new PressureInverseOperator(*coarseOperator_, prm_, tp.getCoarseLevelCommunication());
            return inv;
        }

    private:
        /** @brief The coarse level operator. */
        std::shared_ptr<Operator> coarseOperator_;
        pt::ptree prm_;
    };
} // namespace Amg
} // namespace Dune

#endif // EWOMS_PRESSURE_SOLVER_POLICY_HH
