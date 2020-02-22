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
#ifndef EWOMS_EXTRACTPARALLELGRIDINFORMATIONTOISTL_HH
#define EWOMS_EXTRACTPARALLELGRIDINFORMATIONTOISTL_HH
#ifdef HAVE_EWOMS_ECLGRIDS

#include<ewoms/eclgrids/cpgrid.hh>
#include <any>

namespace Ewoms
{

/// \brief Extracts the information about the data decomposition from the grid for dune-istl
///
/// In the case that grid is a parallel grid this method will query it to get the information
/// about the data decompoisition and convert it to the format expected by the linear algebra
/// of dune-istl.
/// \warn if there is no support for dune-istl and MPI then this functio does not do anything.
/// \param[in] grid The grid to inspect.
/// \param[out] anyComm The handle to to store the information in. If grid is a parallel grid
/// then this will ecapsulate an instance of ParallelISTLInformation.

void extractParallelGridInformationToISTL(const Dune::CpGrid& grid, std::any& anyComm);

// Grid is not CpGrid --> do nothing.
template <class Grid>
void extractParallelGridInformationToISTL(const Grid&, std::any&)
{}

} // end namespace Ewoms
#endif //defined(HAVE_EWOMS_ECLGRIDS)
#endif // EWOMS_EXTRACTPARALLELGRIDINFORMATIONTOISTL_HH
