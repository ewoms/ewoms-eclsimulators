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

#ifndef EWOMS_WELLHELPERS_HH
#define EWOMS_WELLHELPERS_HH

#include <ewoms/eclgrids/gridhelpers.hh>

#include <vector>

namespace Ewoms {

    namespace wellhelpers
    {

        inline
        double computeHydrostaticCorrection(const double well_ref_depth, const double vfp_ref_depth,
                                            const double rho, const double gravity) {
            const double dh = vfp_ref_depth - well_ref_depth;
            const double dp = rho * gravity * dh;

            return dp;
        }

        template <int dim, class C2F, class FC>
        std::array<double, dim>
        getCubeDim(const C2F& c2f,
                   FC         begin_face_centroids,
                   int        cell)
        {
            std::array< std::vector<double>, dim > X;
            {
                const std::vector<double>::size_type
                    nf = std::distance(c2f[cell].begin(),
                                       c2f[cell].end  ());

                for (int d = 0; d < dim; ++d) {
                    X[d].reserve(nf);
                }
            }

            typedef typename C2F::row_type::const_iterator FI;

            for (FI f = c2f[cell].begin(), e = c2f[cell].end(); f != e; ++f) {
                using Ewoms::UgGridHelpers::increment;
                using Ewoms::UgGridHelpers::getCoordinate;

                const FC& fc = increment(begin_face_centroids, *f, dim);

                for (int d = 0; d < dim; ++d) {
                    X[d].push_back(getCoordinate(fc, d));
                }
            }

            std::array<double, dim> cube;
            for (int d = 0; d < dim; ++d) {
                typedef std::vector<double>::iterator VI;
                typedef std::pair<VI,VI>              PVI;

                const PVI m = std::minmax_element(X[d].begin(), X[d].end());

                cube[d] = *m.second - *m.first;
            }

            return cube;
        }

    } // namespace wellhelpers

}

#endif
