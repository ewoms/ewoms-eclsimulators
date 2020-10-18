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
#ifndef EWOMS_FLEXIBLE_SOLVER_HH
#define EWOMS_FLEXIBLE_SOLVER_HH

#include <ewoms/eclsimulators/linalg/preconditionerfactory.hh>

#include <dune/istl/solver.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <boost/property_tree/ptree.hpp>

namespace Dune
{

/// A solver class that encapsulates all needed objects for a linear solver
/// (operator, scalar product, iterative solver and preconditioner) and sets
/// them up based on runtime parameters, using the PreconditionerFactory for
/// setting up preconditioners.
template <class MatrixTypeT, class VectorTypeT>
class FlexibleSolver : public Dune::InverseOperator<VectorTypeT, VectorTypeT>
{
public:
    using MatrixType = MatrixTypeT;
    using VectorType = VectorTypeT;

    /// Base class type of the operator passed to the solver.
    using AbstractOperatorType = Dune::AssembledLinearOperator<MatrixType, VectorType, VectorType>;
    /// Base class type of the contained preconditioner.
    using AbstractPrecondType = Dune::PreconditionerWithUpdate<VectorType, VectorType>;

    /// Create a sequential solver.
    FlexibleSolver(AbstractOperatorType& op,
                   const boost::property_tree::ptree& prm,
                   const std::function<VectorType()>& weightsCalculator = std::function<VectorType()>());

    /// Create a parallel solver (if Comm is e.g. OwnerOverlapCommunication).
    template <class Comm>
    FlexibleSolver(AbstractOperatorType& op,
                   const Comm& comm,
                   const boost::property_tree::ptree& prm,
                   const std::function<VectorType()>& weightsCalculator = std::function<VectorType()>());

    virtual void apply(VectorType& x, VectorType& rhs, Dune::InverseOperatorResult& res) override;

    virtual void apply(VectorType& x, VectorType& rhs, double reduction, Dune::InverseOperatorResult& res) override;

    /// Access the contained preconditioner.
    AbstractPrecondType& preconditioner();

    virtual Dune::SolverCategory::Category category() const override;

private:
    using AbstractScalarProductType = Dune::ScalarProduct<VectorType>;
    using AbstractSolverType = Dune::InverseOperator<VectorType, VectorType>;

    // Machinery for making sequential or parallel operators/preconditioners/scalar products.
    template <class Comm>
    void initOpPrecSp(AbstractOperatorType& op, const boost::property_tree::ptree& prm,
                      const std::function<VectorType()> weightsCalculator, const Comm& comm);

    void initOpPrecSp(AbstractOperatorType& op, const boost::property_tree::ptree& prm,
                      const std::function<VectorType()> weightsCalculator, const Dune::Amg::SequentialInformation&);

    void initSolver(const boost::property_tree::ptree& prm, const bool is_iorank);

    // Main initialization routine.
    // Call with Comm == Dune::Amg::SequentialInformation to get a serial solver.
    template <class Comm>
    void init(AbstractOperatorType& op,
              const Comm& comm,
              const boost::property_tree::ptree& prm,
              const std::function<VectorType()> weightsCalculator);

    AbstractOperatorType* linearoperator_for_solver_;
    std::shared_ptr<AbstractOperatorType> linearoperator_for_precond_;
    std::shared_ptr<AbstractPrecondType> preconditioner_;
    std::shared_ptr<AbstractScalarProductType> scalarproduct_;
    std::shared_ptr<AbstractSolverType> linsolver_;
};

} // namespace Dune

#endif // EWOMS_FLEXIBLE_SOLVER_HH
