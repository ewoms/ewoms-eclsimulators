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

#include <ewoms/eclsimulators/wells/wellprodindexcalculator.hh>

#include <ewoms/eclio/parser/eclipsestate/schedule/well/connection.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/well.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/wellconnections.hh>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
    void checkSizeCompatibility(const Ewoms::WellProdIndexCalculator& wellPICalc,
                                const std::vector<double>&          connMobility)
    {
        if (connMobility.size() != wellPICalc.numConnections()) {
            throw std::logic_error {
                "Mobility vector size does not match expected number of connections"
            };
        }
    }

    double logRescale(const double r0, const double rw, const double rd, const double S)
    {
        const auto numerator = std::log(r0 / rw) + S;
        const auto denom     = std::log(rd / rw) + S;

        return numerator / denom;
    }

    void standardConnFactorsExplicitDrainRadius(const Ewoms::Well&     well,
                                                std::vector<double>& stdConnFact)
    {
        const auto& connections = well.getConnections();
        const auto rdrain = well.getDrainageRadius();

        std::transform(connections.begin(), connections.end(), stdConnFact.begin(),
            [rdrain](const Ewoms::Connection& conn)
        {
            return conn.CF() * logRescale(conn.r0(), conn.rw(), rdrain, conn.skinFactor());
        });
    }

    void standardConnFactorsDrainIsEquivalent(const Ewoms::Well&     well,
                                              std::vector<double>& stdConnFact)
    {
        const auto& connections = well.getConnections();

        std::transform(connections.begin(), connections.end(), stdConnFact.begin(),
            [](const Ewoms::Connection& conn)
        {
            return conn.CF();
        });
    }

    std::vector<double> calculateStandardConnFactors(const Ewoms::Well& well)
    {
        std::vector<double> stdConnFact(well.getConnections().size());

        if (well.getDrainageRadius() > 0.0) {
            // Well has an explicit drainage radius.  Apply logarithmic
            // scaling to the CTFs.
            standardConnFactorsExplicitDrainRadius(well, stdConnFact);
        }
        else {
            // Unspecified drainage radius.  Standard mobility connection
            // scaling factors are just the regular CTFs.
            standardConnFactorsDrainIsEquivalent(well, stdConnFact);
        }

        return stdConnFact;
    }
} // namespace Anonymous

Ewoms::WellProdIndexCalculator::WellProdIndexCalculator(const Well& well)
    : standardConnFactors_{ calculateStandardConnFactors(well) }
{}

void Ewoms::WellProdIndexCalculator::reInit(const Well& well)
{
    this->standardConnFactors_ = calculateStandardConnFactors(well);
}

double
Ewoms::WellProdIndexCalculator::
connectionProdIndStandard(const std::size_t connIdx,
                          const double      connMobility) const
{
    return this->standardConnFactors_[connIdx] * connMobility;
}

// ===========================================================================

std::vector<double>
Ewoms::connectionProdIndStandard(const WellProdIndexCalculator& wellPICalc,
                               const std::vector<double>&     connMobility)
{
    checkSizeCompatibility(wellPICalc, connMobility);

    const auto nConn = wellPICalc.numConnections();
    auto connPI = connMobility;
    for (auto connIx = 0*nConn; connIx < nConn; ++connIx) {
        connPI[connIx] = wellPICalc
            .connectionProdIndStandard(connIx, connMobility[connIx]);
    }

    return connPI;
}

double Ewoms::wellProdIndStandard(const WellProdIndexCalculator& wellPICalc,
                                const std::vector<double>&     connMobility)
{
    const auto connPI = connectionProdIndStandard(wellPICalc, connMobility);

    return std::accumulate(connPI.begin(), connPI.end(), 0.0);
}
