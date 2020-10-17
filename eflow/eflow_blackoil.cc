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
#include "config.h"

#include <eflow/eflow_blackoil.hh>

#include <ewoms/common/resetlocale.hh>
#include <ewoms/eclgrids/cpgrid.hh>
#include <ewoms/eclsimulators/eflow/simulatorfullyimplicitblackoil.hh>
#include <ewoms/eclsimulators/eflow/eflowmain.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#else
#include <dune/common/parallel/mpihelper.hh>
#endif

namespace Ewoms {

void eflowBlackoilSetDeck(double setupTime, Deck& deck, EclipseState& eclState, Schedule& schedule, SummaryConfig& summaryConfig)
{
    using TypeTag = TTAG(EclEFlowProblem);
    using Vanguard = GET_PROP_TYPE(TypeTag, Vanguard);

    Vanguard::setExternalSetupTime(setupTime);
    Vanguard::setExternalDeck(&deck);
    Vanguard::setExternalEclState(&eclState);
    Vanguard::setExternalSchedule(&schedule);
    Vanguard::setExternalSummaryConfig(&summaryConfig);
}

std::unique_ptr<Ewoms::EFlowMain<TTAG(EclEFlowProblem)>>
eflowBlackoilMainInit(int argc, char** argv, bool outputCout, bool outputFiles)
{
    // we always want to use the default locale, and thus spare us the trouble
    // with incorrect locale settings.
    Ewoms::resetLocale();

#if HAVE_DUNE_FEM
    Dune::Fem::MPIManager::initialize(argc, argv);
#else
    Dune::MPIHelper::instance(argc, argv);
#endif

    return std::make_unique<Ewoms::EFlowMain<TTAG(EclEFlowProblem)>>(
        argc, argv, outputCout, outputFiles);
}

// ----------------- Main program -----------------
int eflowBlackoilMain(int argc, char** argv, bool outputCout, bool outputFiles)
{
    auto mainfunc = eflowBlackoilMainInit(argc, argv, outputCout, outputFiles);
    return mainfunc->execute();
}

}
