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
/*!
 * \file
 *
 * \copydoc Ewoms::EclNewtonMethod
 */
#ifndef EWOMS_ECL_NEWTON_METHOD_HH
#define EWOMS_ECL_NEWTON_METHOD_HH

#include <ewoms/numerics/models/blackoil/blackoilnewtonmethod.hh>
#include <ewoms/common/signum.hh>
#include <ewoms/eclio/opmlog/opmlog.hh>

#include <ewoms/common/unused.hh>

BEGIN_PROPERTIES

NEW_PROP_TAG(EclNewtonSumTolerance);
NEW_PROP_TAG(EclNewtonStrictIterations);
NEW_PROP_TAG(EclNewtonRelaxedVolumeFraction);
NEW_PROP_TAG(EclNewtonSumToleranceExponent);
NEW_PROP_TAG(EclNewtonRelaxedTolerance);

END_PROPERTIES

namespace Ewoms {

/*!
 * \brief A newton solver which is eebos specific.
 */
template <class TypeTag>
class EclNewtonMethod : public BlackOilNewtonMethod<TypeTag>
{
    typedef BlackOilNewtonMethod<TypeTag> ParentType;
    using DiscNewtonMethod = GET_PROP_TYPE(TypeTag, DiscNewtonMethod);

    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using SolutionVector = GET_PROP_TYPE(TypeTag, SolutionVector);
    using GlobalEqVector = GET_PROP_TYPE(TypeTag, GlobalEqVector);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

    static const unsigned numEq = GET_PROP_VALUE(TypeTag, NumEq);


    friend NewtonMethod<TypeTag>;
    friend DiscNewtonMethod;
    friend ParentType;

public:
    EclNewtonMethod(Simulator& simulator) : ParentType(simulator)
    {
        errorPvFraction_ = 1.0;
        relaxedMaxPvFraction_ = EWOMS_GET_PARAM(TypeTag, Scalar, EclNewtonRelaxedVolumeFraction);

        sumTolerance_ = 0.0; // this gets determined in the error calculation proceedure
        relaxedTolerance_ = EWOMS_GET_PARAM(TypeTag, Scalar, EclNewtonRelaxedTolerance);

        numStrictIterations_ = EWOMS_GET_PARAM(TypeTag, int, EclNewtonStrictIterations);
    }

    /*!
     * \brief Register all run-time parameters for the Newton method.
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclNewtonSumTolerance,
                             "The maximum error tolerated by the Newton"
                             "method for considering a solution to be "
                             "converged");
        EWOMS_REGISTER_PARAM(TypeTag, int, EclNewtonStrictIterations,
                             "The number of Newton iterations where the"
                             " volumetric error is considered.");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclNewtonRelaxedVolumeFraction,
                             "The fraction of the pore volume of the reservoir "
                             "where the volumetric error may be voilated during "
                             "strict Newton iterations.");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclNewtonSumToleranceExponent,
                             "The the exponent used to scale the sum tolerance by "
                             "the total pore volume of the reservoir.");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, EclNewtonRelaxedTolerance,
                             "The maximum error which the volumetric residual "
                             "may exhibit if it is in a 'relaxed' "
                             "region during a strict iteration.");
    }

    /*!
     * \brief Returns true if the error of the solution is below the
     *        tolerance.
     */
    bool converged() const
    {
        if (errorPvFraction_ < relaxedMaxPvFraction_)
            return (this->error_ < relaxedTolerance_ && errorSum_ < sumTolerance_) ;
        else if (this->numIterations() > numStrictIterations_)
            return (this->error_ < relaxedTolerance_ && errorSum_ < sumTolerance_) ;

        return this->error_ <= this->tolerance() && errorSum_ <= sumTolerance_;
    }

    void preSolve_(const SolutionVector& currentSolution  EWOMS_UNUSED,
                   const GlobalEqVector& currentResidual)
    {
        const auto& constraintsMap = this->model().linearizer().constraintsMap();
        this->lastError_ = this->error_;
        Scalar newtonMaxError = EWOMS_GET_PARAM(TypeTag, Scalar, NewtonMaxError);

        // calculate the error as the maximum weighted tolerance of
        // the solution's residual
        this->error_ = 0.0;
        Dune::FieldVector<Scalar, numEq> componentSumError;
        std::fill(componentSumError.begin(), componentSumError.end(), 0.0);
        Scalar sumPv = 0.0;
        errorPvFraction_ = 0.0;
        const Scalar dt = this->simulator_.timeStepSize();
        for (unsigned dofIdx = 0; dofIdx < currentResidual.size(); ++dofIdx) {
            // do not consider auxiliary DOFs for the error
            if (dofIdx >= this->model().numGridDof()
                || this->model().dofTotalVolume(dofIdx) <= 0.0)
                continue;

            if (!this->model().isLocalDof(dofIdx))
                continue;

            // also do not consider DOFs which are constraint
            if (this->enableConstraints_()) {
                if (constraintsMap.count(dofIdx) > 0)
                    continue;
            }

            const auto& r = currentResidual[dofIdx];
            Scalar pvValue =
                this->simulator_.problem().referencePorosity(dofIdx, /*timeIdx=*/0)
                * this->model().dofTotalVolume(dofIdx);
            sumPv += pvValue;
            bool cnvViolated = false;

            Scalar dofVolume = this->model().dofTotalVolume(dofIdx);

            for (unsigned eqIdx = 0; eqIdx < r.size(); ++eqIdx) {
                Scalar tmpError = r[eqIdx] * dt * this->model().eqWeight(dofIdx, eqIdx) / pvValue;
                Scalar tmpError2 = r[eqIdx] * this->model().eqWeight(dofIdx, eqIdx);

                // in the case of a volumetric formulation, the residual in the above is
                // per cubic meter
                if (GET_PROP_VALUE(TypeTag, UseVolumetricResidual)) {
                    tmpError *= dofVolume;
                    tmpError2 *= dofVolume;
                }

                this->error_ = Ewoms::max(std::abs(tmpError), this->error_);

                if (std::abs(tmpError) > this->tolerance_)
                    cnvViolated = true;

                componentSumError[eqIdx] += std::abs(tmpError2);
            }
            if (cnvViolated)
                errorPvFraction_ += pvValue;
        }

        // take the other processes into account
        this->error_ = this->comm_.max(this->error_);
        componentSumError = this->comm_.sum(componentSumError);
        sumPv = this->comm_.sum(sumPv);
        errorPvFraction_ = this->comm_.sum(errorPvFraction_);

        errorPvFraction_ /= sumPv;

        errorSum_ = 0;
        for (unsigned eqIdx = 0; eqIdx < numEq; ++eqIdx)
            errorSum_ = std::max(std::abs(componentSumError[eqIdx]), errorSum_);

        // scale the tolerance for the total error with the pore volume. by default, the
        // exponent is 1/3, i.e., cubic root.
        Scalar x = EWOMS_GET_PARAM(TypeTag, Scalar, EclNewtonSumTolerance);
        Scalar y = EWOMS_GET_PARAM(TypeTag, Scalar, EclNewtonSumToleranceExponent);
        sumTolerance_ = x*std::pow(sumPv, y);

        this->endIterMsg() << " (max: " << this->tolerance_ << ", violated for " << errorPvFraction_*100 << "% of the pore volume), aggegate error: " << errorSum_ << " (max: " << sumTolerance_ << ")";

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (this->error_ > newtonMaxError)
            throw Ewoms::NumericalIssue("Newton: Error "+std::to_string(double(this->error_))
                                        +" is larger than maximum allowed error of "
                                        +std::to_string(double(newtonMaxError)));

        // make sure that the error never grows beyond the maximum
        // allowed one
        if (errorSum_ > newtonMaxError)
            throw Ewoms::NumericalIssue("Newton: Sum of the error "+std::to_string(double(errorSum_))
                                        +" is larger than maximum allowed error of "
                                        +std::to_string(double(newtonMaxError)));
    }

    void endIteration_(SolutionVector& nextSolution,
                       const SolutionVector& currentSolution)
    {
        ParentType::endIteration_(nextSolution, currentSolution);
        OpmLog::debug( "Newton iteration " + std::to_string(this->numIterations_) + ""
                  + " error: " + std::to_string(double(this->error_))
                  + this->endIterMsg().str());
        this->endIterMsg().str("");
    }

private:
    Scalar errorPvFraction_;
    Scalar errorSum_;

    Scalar relaxedTolerance_;
    Scalar relaxedMaxPvFraction_;

    Scalar sumTolerance_;

    int numStrictIterations_;
};
} // namespace Ewoms

#endif
