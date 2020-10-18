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
#include <config.h>

#include <ewoms/eclsimulators/linalg/setuppropertytree.hh>

#include <boost/version.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace Ewoms
{

/// Set up a property tree intended for FlexibleSolver by either reading
/// the tree from a JSON file or creating a tree giving the default solver
/// and preconditioner. If the latter, the parameters --linear-solver-reduction,
/// --linear-solver-maxiter and --linear-solver-verbosity are used, but if reading
/// from file the data in the JSON file will override any other options.
template<class TypeTag>
boost::property_tree::ptree
setupPropertyTree(EFlowLinearSolverParameters p) // Note: copying the parameters to potentially override.
{
    std::string conf = p.linsolver_;

    // Get configuration from file.
    if (conf.size() > 5 && conf.substr(conf.size() - 5, 5) == ".json") { // the ends_with() method is not available until C++20
#if BOOST_VERSION / 100 % 1000 > 48
        try {
            boost::property_tree::ptree prm;
            boost::property_tree::read_json(conf, prm);
            return prm;
        }
        catch (...) {
            EWOMS_THROW(std::invalid_argument, "Failed reading linear solver configuration from JSON file " << conf);
        }
#else
        EWOMS_THROW(std::invalid_argument,
                  "--linear-solver-configuration=file.json not supported with "
                      << "boost version. Needs version > 1.48.");
#endif
    }

    // Use CPR configuration.
    if ((conf == "cpr") || (conf == "cpr_trueimpes") || (conf == "cpr_quasiimpes")) {
        if (conf == "cpr") {
            // Treat "cpr" as short cut for the true IMPES variant.
            conf = "cpr_trueimpes";
        }
        if (!EWOMS_PARAM_IS_SET(TypeTag, int, LinearSolverMaxIter)) {
            // Use our own default unless it was explicitly overridden by user.
            p.linear_solver_maxiter_ = 20;
        }
        if (!EWOMS_PARAM_IS_SET(TypeTag, int, CprMaxEllIter)) {
            // Use our own default unless it was explicitly overridden by user.
            p.cpr_max_ell_iter_ = 1;
        }
        return setupCPR(conf, p);
    }

    if (conf == "amg") {
        return setupAMG(conf, p);
    }

    // Use ILU0 configuration.
    if (conf == "ilu0") {
        return setupILU(conf, p);
    }

    // No valid configuration option found.
    EWOMS_THROW(std::invalid_argument,
              conf << " is not a valid setting for --linear-solver-configuration."
              << " Please use ilu0, cpr, cpr_trueimpes, or cpr_quasiimpes");
}

} // namespace Ewoms
