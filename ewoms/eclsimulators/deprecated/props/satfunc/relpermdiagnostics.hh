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

#ifndef EWOMS_RELPERMDIAGNOSTICS_HH
#define EWOMS_RELPERMDIAGNOSTICS_HH

#include <vector>
#include <utility>

#include <ewoms/eclio/utility/numeric/linearinterpolation.hh>
#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>
#include <ewoms/eclio/parser/deck/deck.hh>
#include <ewoms/eclio/opmlog/opmlog.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/ssfntable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/misctable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/msfntable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/swoftable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sgoftable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/slgoftable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/swfntable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sgfntable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sof3table.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sgcwmistable.hh>
#include <ewoms/eclio/parser/eclipsestate/tables/sorwmistable.hh>
#include <ewoms/material/fluidmatrixinteractions/eclepsscalingpoints.hh>

namespace Ewoms {

    class Sof2Table;
    class SgwfnTable;

    ///This class is intend to be a relpmer diganostics, to detect
    ///wrong input of relperm table and endpoints.
    class RelpermDiagnostics
    {
    public:
        ///This function is used to diagnosis relperm in
        ///eclipse data file. Errors and warings will be
        ///output if they're found.
        ///\param[in] eclState  eclipse state.
        ///\param[in] deck      ecliplise data file.
        ///\param[in] grid      unstructured grid.
        template <class GridT>
        void diagnosis(const EclipseState& eclState,
                       const Deck& deck,
                       const GridT& grid);

    private:
        enum FluidSystem {
            OilWater,
            OilGas,
            WaterGas,
            BlackOil,
            Solvent
        };

        FluidSystem fluidSystem_;

        enum SaturationFunctionFamily {
            FamilyI,
            FamilyII,
            NoFamily
        };

        SaturationFunctionFamily satFamily_;

        std::vector<Ewoms::EclEpsScalingPointsInfo<double> > unscaledEpsInfo_;
        std::vector<Ewoms::EclEpsScalingPointsInfo<double> > scaledEpsInfo_;

        ///Check the phase that used.
        void phaseCheck_(const EclipseState& es);

        ///Check saturation family I and II.
        void satFamilyCheck_(const EclipseState& eclState);

        ///Check saturation tables.
        void tableCheck_(const EclipseState& eclState);

        ///Check endpoints in the saturation tables.
        void unscaledEndPointsCheck_(const EclipseState& eclState);

        template <class GridT>
        void scaledEndPointsCheck_(const Deck& deck,
                                   const EclipseState& eclState,
                                   const GridT& grid);

        ///For every table, need to deal with case by case.
        void swofTableCheck_(const Ewoms::SwofTable& swofTables,
                             const int satnumIdx);
        void sgofTableCheck_(const Ewoms::SgofTable& sgofTables,
                             const int satnumIdx);
        void slgofTableCheck_(const Ewoms::SlgofTable& slgofTables,
                              const int satnumIdx);
        void swfnTableCheck_(const Ewoms::SwfnTable& swfnTables,
                             const int satnumIdx);
        void sgfnTableCheck_(const Ewoms::SgfnTable& sgfnTables,
                             const int satnumIdx);
        void sof3TableCheck_(const Ewoms::Sof3Table& sof3Tables,
                             const int satnumIdx);
        void sof2TableCheck_(const Ewoms::Sof2Table& sof2Tables,
                             const int satnumIdx);
        void sgwfnTableCheck_(const Ewoms::SgwfnTable& sgwfnTables,
                              const int satnumIdx);
        ///Tables for solvent model
        void sgcwmisTableCheck_(const Ewoms::SgcwmisTable& sgcwmisTables,
                                const int satnumIdx);
        void sorwmisTableCheck_(const Ewoms::SorwmisTable& sorwmisTables,
                                const int satnumIdx);
        void ssfnTableCheck_(const Ewoms::SsfnTable& ssfnTables,
                             const int satnumIdx);
        void miscTableCheck_(const Ewoms::MiscTable& miscTables,
                             const int miscnumIdx);
        void msfnTableCheck_(const Ewoms::MsfnTable& msfnTables,
                             const int satnumIdx);
    };

} //namespace Ewoms

#include <ewoms/eclsimulators/deprecated/props/satfunc/relpermdiagnostics_impl.hh>

#endif // EWOMS_RELPERMDIAGNOSTICS_HH
