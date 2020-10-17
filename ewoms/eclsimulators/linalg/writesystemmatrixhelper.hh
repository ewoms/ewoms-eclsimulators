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
#ifndef EWOMS_WRITESYSTEMMATRIXHELPER_HH
#define EWOMS_WRITESYSTEMMATRIXHELPER_HH

#include <dune/istl/matrixmarket.hh>

namespace Ewoms
{
namespace Helper
{
    template <class SimulatorType, class MatrixType, class VectorType, class Communicator>
    void writeSystem(const SimulatorType& simulator,
                     const MatrixType& matrix,
                     const VectorType& rhs,
                     [[maybe_unused]] const Communicator* comm)
    {
        std::string dir = simulator.problem().outputDir();
        if (dir == ".") {
            dir = "";
        } else if (!dir.empty() && dir.back() != '/') {
            dir += "/";
        }
        namespace fs = Ewoms::filesystem;
        fs::path output_dir(dir);
        fs::path subdir("reports");
        output_dir = output_dir / subdir;
        if (!(fs::exists(output_dir))) {
            fs::create_directory(output_dir);
        }
        // Combine and return.
        std::ostringstream oss;
        oss << "prob_" << simulator.episodeIndex() << "_time_";
        oss << std::setprecision(15) << std::setw(12) << std::setfill('0') << simulator.time() << "_";
        int nit = simulator.model().newtonMethod().numIterations();
        oss << "_nit_" << nit << "_";
        std::string output_file(oss.str());
        fs::path full_path = output_dir / output_file;
        std::string prefix = full_path.string();
        {
            std::string filename = prefix + "matrix_istl";
#if HAVE_MPI
            if (comm != nullptr) { // comm is not set in serial runs
                Dune::storeMatrixMarket(matrix, filename, *comm, true);
            } else
#endif
            {
                Dune::storeMatrixMarket(matrix, filename + ".mm");
            }
        }
        {
            std::string filename = prefix + "rhs_istl";
#if HAVE_MPI
            if (comm != nullptr) { // comm is not set in serial runs
                Dune::storeMatrixMarket(rhs, filename, *comm, true);
            } else
#endif
            {
                Dune::storeMatrixMarket(rhs, filename + ".mm");
            }
        }
    }

} // namespace Helper
} // namespace Ewoms
#endif
