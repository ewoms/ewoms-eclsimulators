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

#ifndef EWOMS_MSWELLHELPERS_HH
#define EWOMS_MSWELLHELPERS_HH

#include <ewoms/eclsimulators/utils/deferredloggingerrorhelpers.hh>
#include <ewoms/eclsimulators/utils/deferredlogger.hh>
#include <ewoms/eclio/errormacros.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/msw/spiralicd.hh>
#include <dune/istl/solvers.hh>
#if HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif // HAVE_UMFPACK
#include <cmath>

namespace Ewoms {

namespace mswellhelpers
{
    // obtain y = D^-1 * x with a direct solver
    template <typename MatrixType, typename VectorType>
    VectorType
    invDXDirect(const MatrixType& D, VectorType x)
    {
#if HAVE_UMFPACK
        VectorType y(x.size());
        y = 0.;

        Dune::UMFPack<MatrixType> linsolver(D, 0);

        // Object storing some statistics about the solving process
        Dune::InverseOperatorResult res;

        // Solve
        linsolver.apply(y, x, res);

        // Checking if there is any inf or nan in y
        // it will be the solution before we find a way to catch the singularity of the matrix
        for (size_t i_block = 0; i_block < y.size(); ++i_block) {
            for (size_t i_elem = 0; i_elem < y[i_block].size(); ++i_elem) {
                if (std::isinf(y[i_block][i_elem]) || std::isnan(y[i_block][i_elem]) ) {
                    EWOMS_THROW(Ewoms::NumericalIssue, "nan or inf value found in invDXDirect due to singular matrix");
                }
            }
        }

        return y;
#else
        // this is not thread safe
        EWOMS_THROW(std::runtime_error, "Cannot use invDXDirect() without UMFPACK. "
                  "Reconfigure ewoms-eclsimulators with SuiteSparse/UMFPACK support and recompile.");
#endif // HAVE_UMFPACK
    }

    // obtain y = D^-1 * x with a BICSSTAB iterative solver
    template <typename MatrixType, typename VectorType>
    VectorType
    invDX(const MatrixType& D, VectorType x, Ewoms::DeferredLogger& deferred_logger)
    {
        // the function will change the value of x, so we should not use reference of x here.

        // TODO: store some of the following information to avoid to call it again and again for
        // efficiency improvement.
        // Bassically, only the solve / apply step is different.

        VectorType y(x.size());
        y = 0.;

        Dune::MatrixAdapter<MatrixType, VectorType, VectorType> linearOperator(D);

        // Sequential incomplete LU decomposition as the preconditioner
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2,7)
        Dune::SeqILU<MatrixType, VectorType, VectorType> preconditioner(D, 1.0);
#else
        Dune::SeqILU0<MatrixType, VectorType, VectorType> preconditioner(D, 1.0);
#endif
        // Dune::SeqILUn<MatrixType, VectorType, VectorType> preconditioner(D, 1, 0.92);
        // Dune::SeqGS<MatrixType, VectorType, VectorType> preconditioner(D, 1, 1);
        // Dune::SeqJac<MatrixType, VectorType, VectorType> preconditioner(D, 1, 1);

        // Preconditioned BICGSTAB solver
        Dune::BiCGSTABSolver<VectorType> linsolver(linearOperator,
                                                   preconditioner,
                                                   1.e-8, // desired residual reduction factor
                                                   250, // maximum number of iterations
                                                   0); // verbosity of the solver */

        // Object storing some statistics about the solving process
        Dune::InverseOperatorResult res;

        // Solve
        linsolver.apply(y, x, res);

        if ( !res.converged ) {
            EWOMS_DEFLOG_THROW(Ewoms::NumericalIssue, "the invDX does not get converged! ", deferred_logger);
        }

        return y;
    }

    template <typename ValueType>
    inline ValueType haalandFormular(const ValueType& re, const double diameter, const double roughness)
    {
        const ValueType value = -3.6 * Ewoms::log10(6.9 / re + std::pow(roughness / (3.7 * diameter), 10. / 9.) );

        // sqrt(1/f) should be non-positive
        assert(value >= 0.0);

        return 1. / (value * value);
    }

    template <typename ValueType>
    inline ValueType calculateFrictionFactor(const double area, const double diameter,
                                          const ValueType& w, const double roughness, const ValueType& mu)
    {

        ValueType f = 0.;
        // Reynolds number
        const ValueType re = Ewoms::abs( diameter * w / (area * mu));

        if ( re == 0.0 ) {
            // make sure it is because the mass rate is zero
            assert(w == 0.);
            return 0.0;
        }

        const ValueType re_value1 = 2000.;
        const ValueType re_value2 = 4000.;

        if (re < re_value1) {
            f = 16. / re;
        } else if (re > re_value2){
            f = haalandFormular(re, diameter, roughness);
        } else { // in between
            const ValueType f1 = 16. / re_value1;
            const ValueType f2 = haalandFormular(re_value2, diameter, roughness);

            f = (f2 - f1) / (re_value2 - re_value1) * (re - re_value1) + f1;
        }
        return f;
    }

    // calculating the friction pressure loss
    // l is the segment length
    // area is the segment cross area
    // diameter is the segment inner diameter
    // w is mass flow rate through the segment
    // density is density
    // roughness is the absolute roughness
    // mu is the average phase viscosity
    template <typename ValueType>
    ValueType frictionPressureLoss(const double l, const double diameter, const double area, const double roughness,
                                   const ValueType& density, const ValueType& w, const ValueType& mu)
    {
        const ValueType f = calculateFrictionFactor(area, diameter, w, roughness, mu);
        // \Note: a factor of 2 needs to be here based on the dimensional analysis
        return 2. * f * l * w * w / (area * area * diameter * density);
    }

    template <typename ValueType>
    ValueType velocityHead(const double area, const ValueType& mass_rate, const ValueType& density)
    {
        return (0.5 * mass_rate * mass_rate / (area * area * density));
    }

    // water in oil emulsion viscosity
    // TODO: maybe it should be two different ValueTypes. When we calculate the viscosity for transitional zone
    template <typename ValueType>
    ValueType WIOEmulsionViscosity(const ValueType& oil_viscosity, const ValueType& water_liquid_fraction,
                                   const double max_visco_ratio)
    {
        const ValueType temp_value = 1. / (1. - (0.8415 / 0.7480 * water_liquid_fraction) );
        const ValueType viscosity_ratio = Ewoms::pow(temp_value, 2.5);

        if (viscosity_ratio <= max_visco_ratio) {
            return oil_viscosity * viscosity_ratio;
        } else {
            return oil_viscosity * max_visco_ratio;
        }
    }

    // oil in water emulsion viscosity
    template <typename ValueType>
    ValueType OIWEmulsionViscosity(const ValueType& water_viscosity, const ValueType& water_liquid_fraction,
                                   const double max_visco_ratio)
    {
        const ValueType temp_value = 1. / (1. - (0.6019 / 0.6410) * (1. - water_liquid_fraction) );
        const ValueType viscosity_ratio = Ewoms::pow(temp_value, 2.5);

        if (viscosity_ratio <= max_visco_ratio) {
            return water_viscosity * viscosity_ratio;
        } else {
            return water_viscosity * max_visco_ratio;
        }
    }

    // calculating the viscosity of oil-water emulsion at local conditons
    template <typename ValueType>
    ValueType emulsionViscosity(const ValueType& water_fraction, const ValueType& water_viscosity,
                                const ValueType& oil_fraction, const ValueType& oil_viscosity,
                                const SpiralICD& sicd)
    {
        const double width_transition = sicd.widthTransitionRegion();

        // it is just for now, we should be able to treat it.
        if (width_transition <= 0.) {
            EWOMS_THROW(std::runtime_error, "Not handling non-positive transition width now");
        }

        const double critical_value = sicd.criticalValue();
        const ValueType transition_start_value = critical_value - width_transition / 2.0;
        const ValueType transition_end_value = critical_value + width_transition / 2.0;

        const ValueType liquid_fraction = water_fraction + oil_fraction;
        // if there is no liquid, we just return zero
        if (liquid_fraction == 0.) {
            return 0.;
        }

        const ValueType water_liquid_fraction = water_fraction / liquid_fraction;

        const double max_visco_ratio = sicd.maxViscosityRatio();
        if (water_liquid_fraction <= transition_start_value) {
            return WIOEmulsionViscosity(oil_viscosity, water_liquid_fraction, max_visco_ratio);
        } else if(water_liquid_fraction >= transition_end_value) {
            return OIWEmulsionViscosity(water_viscosity, water_liquid_fraction, max_visco_ratio);
        } else { // in the transition region
            const ValueType viscosity_start_transition = WIOEmulsionViscosity(oil_viscosity, transition_start_value, max_visco_ratio);
            const ValueType viscosity_end_transition = OIWEmulsionViscosity(water_viscosity, transition_end_value, max_visco_ratio);
            const ValueType emulsion_viscosity = (viscosity_start_transition * (transition_end_value - water_liquid_fraction)
                                               + viscosity_end_transition * (water_liquid_fraction - transition_start_value) ) / width_transition;
            return emulsion_viscosity;
        }
    }

} // namespace mswellhelpers

}

#endif
