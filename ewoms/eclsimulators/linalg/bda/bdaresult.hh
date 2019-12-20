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

#ifndef BDARESULT_HH
#define BDARESULT_HH

namespace Ewoms
{

/// This class is based on InverseOperatorResult struct from dune/istl/solver.hh
/// It is needed to prevent a compile error in basearray.hh, the nvcc compiler might not support all features in there
class BdaResult
{

public:
    int iterations = 0;         // number of iterations
    double reduction = 0.0;     // reduction of norm, norm_start / norm_final
    bool converged = false;     // true iff the linear solver reached the desired norm within maxit iterations
    double conv_rate = 0.0;     // average reduction of norm per iteration, usually calculated with 'static_cast<double>(pow(res.reduction,1.0/it));'
    double elapsed = 0.0;       // time in seconds to run the linear solver

    // Dune 2.6 has a member 'double condition_estimate = -1' in InverseOperatorResult

}; // end class BdaResult

}

#endif
