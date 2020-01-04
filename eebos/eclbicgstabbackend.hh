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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Ewoms::Linear::EclBiCGStabSolverBackend
 */
#ifndef EWOMS_ECL_BICGSTAB_BACKEND_HH
#define EWOMS_ECL_BICGSTAB_BACKEND_HH

#include "eclcriterion.hh"

#include <ewoms/numerics/linear/parallelbicgstabbackend.hh>

#include <memory>

namespace Ewoms {
namespace Linear {
template <class TypeTag>
class EclBiCGStabSolverBackend;
}} // namespace Linear, Ewoms

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclBiCGStabSolverBackend, INHERITS_FROM(ParallelBiCGStabLinearSolver));

NEW_PROP_TAG(LinearSolverMaxError);

SET_TYPE_PROP(EclBiCGStabSolverBackend,
              LinearSolverBackend,
              Ewoms::Linear::EclBiCGStabSolverBackend<TypeTag>);

SET_SCALAR_PROP(EclBiCGStabSolverBackend, LinearSolverTolerance, 5e-2);

END_PROPERTIES

namespace Ewoms {
namespace Linear {
/*!
 * \ingroup Linear
 *
 * \brief Implements a linear solver backend targeted at ECL problems.
 *
 * This class utilizes a parallel stabilized BiCG linar solver with a convergence
 * criterion that is tuned for the non-linear solver of eebos. The preconditioner be
 * specified like for the generic ParallelBiCGStabBackend (the default is ILU-0)
 */
template <class TypeTag>
class EclBiCGStabSolverBackend : public ParallelBiCGStabSolverBackend<TypeTag>
{
    typedef ParallelBiCGStabSolverBackend<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;

    typedef typename ParentType::ParallelOperator ParallelOperator;
    typedef typename ParentType::OverlappingVector OverlappingVector;
    typedef typename ParentType::ParallelPreconditioner ParallelPreconditioner;
    typedef typename ParentType::ParallelScalarProduct ParallelScalarProduct;

    typedef typename ParentType::RawLinearSolver RawLinearSolver;

public:
    EclBiCGStabSolverBackend(const Simulator& simulator)
        : ParentType(simulator)
    { }

    static void registerParameters()
    {
        ParentType::registerParameters();

        // this solver backend ignores the LinearSolverMaxError parameter. To avoid
        // confusion, we prevent it from showing up in the usage message.
        EWOMS_HIDE_PARAM(TypeTag, LinearSolverMaxError);
    }

protected:
    friend ParallelBaseBackend<TypeTag>;
    friend ParentType;

    std::shared_ptr<RawLinearSolver> prepareSolver_(ParallelOperator& parOperator,
                                                    ParallelScalarProduct& parScalarProduct,
                                                    ParallelPreconditioner& parPreCond)
    {
        typedef EclCriterion<TypeTag, OverlappingVector> EclCrit;

        this->convCrit_.reset(new EclCrit(this->simulator_));

        auto bicgstabSolver =
            std::make_shared<RawLinearSolver>(parPreCond, *this->convCrit_, parScalarProduct);

        int verbosity = 0;
        if (parOperator.overlap().myRank() == 0)
            verbosity = EWOMS_GET_PARAM(TypeTag, int, LinearSolverVerbosity);

        bicgstabSolver->setVerbosity(verbosity);
        bicgstabSolver->setMaxIterations(EWOMS_GET_PARAM(TypeTag, int, LinearSolverMaxIterations));
        bicgstabSolver->setLinearOperator(&parOperator);
        bicgstabSolver->setRhs(this->overlappingb_);

        return bicgstabSolver;
    }
};

}} // namespace Linear, Ewoms

#endif
