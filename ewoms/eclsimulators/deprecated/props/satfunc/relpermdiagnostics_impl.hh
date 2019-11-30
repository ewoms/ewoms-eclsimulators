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

#ifndef EWOMS_RELPERMDIAGNOSTICS_IMPL_HH
#define EWOMS_RELPERMDIAGNOSTICS_IMPL_HH

#include <vector>
#include <utility>

#include <ewoms/material/fluidmatrixinteractions/eclepsgridproperties.hh>
#include <ewoms/eclsimulators/deprecated/props/satfunc/relpermdiagnostics.hh>
#include <ewoms/eclgrids/utility/compressedtocartesian.hh>
#include <ewoms/eclgrids/gridhelpers.hh>

namespace Ewoms {

    template <class GridT>
    void RelpermDiagnostics::diagnosis(const Ewoms::EclipseState& eclState,
                                       const Ewoms::Deck& deck,
                                       const GridT& grid)
    {
        OpmLog::info("\n===============Saturation Functions Diagnostics===============\n");
        phaseCheck_(eclState);
        satFamilyCheck_(eclState);
        tableCheck_(eclState);
        unscaledEndPointsCheck_(deck, eclState);
        scaledEndPointsCheck_(deck, eclState, grid);
    }

    template <class GridT>
    void RelpermDiagnostics::scaledEndPointsCheck_(const Deck& deck,
                                                   const EclipseState& eclState,
                                                   const GridT& grid)
    {
        // All end points are subject to round-off errors, checks should account for it
        const float tolerance = 1e-6;
        const int nc = Ewoms::UgGridHelpers::numCells(grid);
        const auto& global_cell = Ewoms::UgGridHelpers::globalCell(grid);
        const auto dims = Ewoms::UgGridHelpers::cartDims(grid);
        const auto& compressedToCartesianIdx = Ewoms::compressedToCartesian(nc, global_cell);
        scaledEpsInfo_.resize(nc);
        EclEpsGridProperties epsGridProperties(eclState, false, compressedToCartesianIdx);
        const std::string tag = "Scaled endpoints";
        for (int c = 0; c < nc; ++c) {
            const std::string satnumIdx = std::to_string(epsGridProperties.satRegion(c));
            std::string cellIdx;
            {
                std::array<int, 3> ijk;
                const int cartIdx = compressedToCartesianIdx[c];
                ijk[0] = cartIdx % dims[0];
                ijk[1] = (cartIdx / dims[0]) % dims[1];
                ijk[2] = cartIdx / dims[0] / dims[1];
                cellIdx = "(" + std::to_string(ijk[0]) + ", " +
                    std::to_string(ijk[1]) + ", " +
                    std::to_string(ijk[2]) + ")";
            }
            scaledEpsInfo_[c].extractScaled(eclState, epsGridProperties, c);

            // SGU <= 1.0 - SWL
            if (scaledEpsInfo_[c].Sgu > (1.0 - scaledEpsInfo_[c].Swl + tolerance)) {
                const std::string msg = "For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SGU exceed 1.0 - SWL";
                OpmLog::warning(tag, msg);
            }

            // SGL <= 1.0 - SWU
            if (scaledEpsInfo_[c].Sgl > (1.0 - scaledEpsInfo_[c].Swu + tolerance)) {
                const std::string msg = "For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SGL exceed 1.0 - SWU";
                OpmLog::warning(tag, msg);
            }

            if (deck.hasKeyword("SCALECRS") && fluidSystem_ == FluidSystem::BlackOil) {
                // Mobilility check.
		    if ((scaledEpsInfo_[c].Sowcr + scaledEpsInfo_[c].Swcr) >= (1.0 + tolerance)) {
                    const std::string msg = "For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SOWCR + SWCR exceed 1.0";
                    OpmLog::warning(tag, msg);
                }

            if ((scaledEpsInfo_[c].Sogcr + scaledEpsInfo_[c].Sgcr + scaledEpsInfo_[c].Swl) >= (1.0 + tolerance)) {
                    const std::string msg = "For scaled endpoints input, cell" + cellIdx + " SATNUM = " + satnumIdx + ", SOGCR + SGCR + SWL exceed 1.0";
                    OpmLog::warning(tag, msg);
                }
            }
        }
    }

} //namespace Ewoms

#endif // EWOMS_RELPERMDIAGNOSTICS_IMPL_HH
