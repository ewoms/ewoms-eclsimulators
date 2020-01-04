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
 * \copydoc Ewoms::Linear::EclCriterion
 */
#ifndef EWOMS_ECL_CRITERION_HH
#define EWOMS_ECL_CRITERION_HH

#include <ewoms/numerics/linear/convergencecriterion.hh>

#include <iostream>

BEGIN_PROPERTIES

NEW_PROP_TAG(UseVolumetricResidual);
NEW_PROP_TAG(EnableConstraints);
NEW_PROP_TAG(Indices);

END_PROPERTIES

namespace Ewoms {
namespace Linear {

/*! \addtogroup Linear
 * \{
 */

/*!
 * \brief Convergence criterion which uses a convergence criterion for the linear solver
 *        that is consistent with the convergence criterion of eebos' non-linear solver.
 */
template <class TypeTag, class Vector>
class EclCriterion : public ConvergenceCriterion<Vector>
{
    typedef typename Vector::field_type Scalar;
    typedef typename Vector::block_type BlockType;

    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename Grid::CollectiveCommunication CollectiveCommunication;

    static const int numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    EclCriterion(const Simulator& simulator)
        : simulator_(simulator)
    {
        residualReductionTolerance_ = EWOMS_GET_PARAM(TypeTag, Scalar, LinearSolverTolerance);
    }

    /*!
     * \brief Sets the residual reduction tolerance.
     */
    void setResidualReductionTolerance(Scalar tol)
    { residualReductionTolerance_ = tol; }

    /*!
     * \brief Returns the tolerance of the residual reduction of the solution.
     */
    Scalar residualReductionTolerance() const
    { return residualReductionTolerance_; }

    /*!
     * \brief Returns the reduction of the maximum of the residual compared to the
     *        initial solution.
     */
    Scalar residualReduction() const
    { return residualError_/initialResidualError_; }

    /*!
     * \brief Sets the maximum absolute tolerated residual.
     */
    void setAbsResidualTolerance(Scalar tol)
    { absResidualTolerance_ = tol; }

    /*!
     * \brief Returns the tolerated maximum of the the infinity norm of the absolute
     *        residual.
     */
    Scalar absResidualTolerance() const
    { return absResidualTolerance_; }

    /*!
     * \brief Returns the infinity norm of the absolute residual.
     */
    Scalar absResidual() const
    { return residualError_; }

    /*!
     * \copydoc ConvergenceCriterion::setInitial(const Vector& , const Vector& )
     */
    void setInitial(const Vector& curSol, const Vector& curResid) override
    {
        updateErrors_(curSol, curSol, curResid);

        // to avoid divisions by zero, make sure that we don't use an initial error of 0
        lastResidualError_ = residualError_;
        initialResidualError_ =
            std::max<Scalar>(residualError_, std::numeric_limits<Scalar>::min()*1e3);
    }

    /*!
     * \copydoc ConvergenceCriterion::update(const Vector&, const Vector&, const Vector&)
     */
    void update(const Vector& curSol, const Vector& changeIndicator, const Vector& curResid) override
    { updateErrors_(curSol, changeIndicator, curResid);  }

    /*!
     * \copydoc ConvergenceCriterion::converged()
     */
    bool converged() const override
    {
        // we're converged if we achieved at least the specified reduction of the residual.
        return residualReduction() <= residualReductionTolerance();
    }

    /*!
     * \copydoc ConvergenceCriterion::failed()
     */
    bool failed() const override
    { return false; }

    /*!
     * \copydoc ConvergenceCriterion::accuracy()
     *
     * For the accuracy we only take the residual reduction into account,
     */
    Scalar accuracy() const override
    { return residualReduction(); }

    /*!
     * \copydoc ConvergenceCriterion::printInitial()
     */
    void printInitial(std::ostream& os = std::cout) const override
    {
        os << std::setw(20) << "iteration ";
        os << std::setw(20) << "residual ";
        os << std::setw(20) << "reduction ";
        os << std::setw(20) << "reductionRate ";
        os << std::endl;
    }

    /*!
     * \copydoc ConvergenceCriterion::print()
     */
    void print(Scalar iter, std::ostream& os = std::cout) const override
    {
        const Scalar eps = std::numeric_limits<Scalar>::min()*1e10;

        os << std::setw(20) << iter << " ";
        os << std::setw(20) << absResidual() << " ";
        os << std::setw(20) << residualReduction() << " ";
        os << std::setw(20) << lastResidualError_/std::max<Scalar>(residualError_, eps) << " ";
        os << std::endl << std::flush;
    }

private:
    // update the weighted absolute residual
    void updateErrors_(const Vector& curSol EWOMS_UNUSED, const Vector& changeIndicator,  const Vector& currentResidual)
    {
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

        const auto& constraintsMap = this->model_().linearizer().constraintsMap();
        this->lastResidualError_ = this->residualError_;

        // calculate the error as the maximum weighted tolerance of
        // the solution's residual
        this->residualError_ = 0.0;
        Dune::FieldVector<Scalar, numEq> componentSumError(0.0);
        size_t n = this->model_().numGridDof();
        for (unsigned dofIdx = 0; dofIdx < n; ++dofIdx) {
            // do not consider auxiliary DOFs for the error
            if (this->model_().dofTotalVolume(dofIdx) <= 0.0)
                continue;

            if (!this->model_().isLocalDof(dofIdx))
                continue;

            // also do not consider DOFs which are constraint
            if (GET_PROP_VALUE(TypeTag, EnableConstraints)) {
                if (constraintsMap.count(dofIdx) > 0)
                    continue;
            }

            const auto& r = currentResidual[dofIdx];
            for (unsigned eqIdx = 0; eqIdx < r.size(); ++eqIdx) {
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                    && eqIdx == Indices::conti0EqIdx + Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx))
                {
                    // ignore water
                    continue;
                }


                Scalar tmpError = r[eqIdx] * this->model_().eqWeight(dofIdx, eqIdx);

                // in the case of a volumetric formulation, the residual in the above is
                // per cubic meter
                if (GET_PROP_VALUE(TypeTag, UseVolumetricResidual))
                    tmpError *= this->model_().dofTotalVolume(dofIdx);

                componentSumError[eqIdx] += std::abs(tmpError);
            }
        }

        // take the other processes into account
        componentSumError = this->comm_().sum(componentSumError);

        residualError_ = 0.0;
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx)
            residualError_ += std::abs(componentSumError[eqIdx]);
    }

    const CollectiveCommunication& comm_() const
    { return simulator_.gridView().comm(); }

    const Model& model_() const
    { return simulator_.model(); }

    const Simulator& simulator_;

    // the total pore volume of the whole reservoir
    Scalar sumPv_;

    // the infinity norm of the residual of the last iteration
    Scalar lastResidualError_;

    // the infinity norm of the residual of the current iteration
    Scalar residualError_;

    // the infinity norm of the residual of the initial solution
    Scalar initialResidualError_;

    // the minimum reduction of the residual norm where the solution is to be considered
    // converged
    Scalar residualReductionTolerance_;

    // the maximum residual norm for the residual for the solution to be considered to be
    // converged
    Scalar absResidualTolerance_;
};

//! \} end documentation

}} // end namespace Linear, Ewoms

#endif
