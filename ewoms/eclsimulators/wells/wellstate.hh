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

#ifndef EWOMS_WELLSTATE_HH
#define EWOMS_WELLSTATE_HH

#include <ewoms/eclsimulators/deprecated/props/blackoilphases.hh>
#include <ewoms/eclio/output/data/wells.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/well/well.hh>
#include <ewoms/eclsimulators/wells/perforationdata.hh>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cassert>
#include <cstddef>

namespace Ewoms
{
    /// The state of a set of wells.
    class WellState
    {
    public:
        typedef std::array< int, 3 >  mapentry_t;
        typedef std::map< std::string, mapentry_t > WellMapType;

        /// Allocate and initialize if wells is non-null.
        /// Also tries to give useful initial values to the bhp() and
        /// wellRates() fields, depending on controls.  The
        /// perfRates() field is filled with zero, and perfPress()
        /// with -1e100.
        void init(const std::vector<double>& cellPressures,
                  const std::vector<Well>& wells_ecl,
                  const PhaseUsage& pu,
                  const std::vector<std::vector<PerforationData>>& well_perf_data,
                  const SummaryState& summary_state)
        {
            // clear old name mapping
            wellMap_.clear();

            well_perf_data_ = well_perf_data;

            {
                // const int nw = wells->number_of_wells;
                const int nw = wells_ecl.size();
                // const int np = wells->number_of_phases;
                const int np = pu.num_phases;
                open_for_output_.assign(nw, true);
                bhp_.resize(nw, 0.0);
                thp_.resize(nw, 0.0);
                temperature_.resize(nw, 273.15 + 20); // standard temperature for now
                wellrates_.resize(nw * np, 0.0);
                int connpos = 0;
                for (int w = 0; w < nw; ++w) {
                    const Well& well = wells_ecl[w];

                    // Initialize bhp(), thp(), wellRates().
                    initSingleWell(cellPressures, w, well, pu, summary_state);

                    // Setup wellname -> well index mapping.
                    const int num_perf_this_well = well_perf_data[w].size();
                    std::string name = well.name();
                    assert( name.size() > 0 );
                    mapentry_t& wellMapEntry = wellMap_[name];
                    wellMapEntry[ 0 ] = w;
                    wellMapEntry[ 1 ] = connpos;
                    wellMapEntry[ 2 ] = num_perf_this_well;
                    connpos += num_perf_this_well;
                }

                // The perforation rates and perforation pressures are
                // not expected to be consistent with bhp_ and wellrates_
                // after init().
                perfrates_.resize(connpos, 0.0);
                perfpress_.resize(connpos, -1e100);
            }
        }

        /// One bhp pressure per well.
        std::vector<double>& bhp() { return bhp_; }
        const std::vector<double>& bhp() const { return bhp_; }

        /// One thp pressure per well.
        std::vector<double>& thp() { return thp_; }
        const std::vector<double>& thp() const { return thp_; }

        /// One temperature per well.
        std::vector<double>& temperature() { return temperature_; }
        const std::vector<double>& temperature() const { return temperature_; }

        /// One rate per well and phase.
        std::vector<double>& wellRates() { return wellrates_; }
        const std::vector<double>& wellRates() const { return wellrates_; }

        /// One rate per well connection.
        std::vector<double>& perfRates() { return perfrates_; }
        const std::vector<double>& perfRates() const { return perfrates_; }

        /// One pressure per well connection.
        std::vector<double>& perfPress() { return perfpress_; }
        const std::vector<double>& perfPress() const { return perfpress_; }

        size_t getRestartBhpOffset() const {
            return 0;
        }

        size_t getRestartPerfPressOffset() const {
            return bhp_.size();
        }

        size_t getRestartPerfRatesOffset() const {
            return getRestartPerfPressOffset() + perfpress_.size();
        }

        size_t getRestartTemperatureOffset() const {
            return getRestartPerfRatesOffset() + perfrates_.size();
        }

        size_t getRestartWellRatesOffset() const {
            return getRestartTemperatureOffset() + temperature_.size();
        }

        const WellMapType& wellMap() const { return wellMap_; }
        WellMapType& wellMap() { return wellMap_; }

        /// The number of wells present.
        int numWells() const
        {
            return bhp().size();
        }

        /// The number of phases present.
        int numPhases() const
        {
            return wellRates().size() / numWells();
        }

        virtual void shutWell(int well_index) {
            this->open_for_output_[well_index] = false;
            this->thp_[well_index] = 0;
            this->bhp_[well_index] = 0;
            const int np = numPhases();
            for (int p = 0; p < np; ++p)
                this->wellrates_[np * well_index + p] = 0;
        }

        virtual data::Wells report(const PhaseUsage& pu, const int* globalCellIdxMap) const
        {
            using rt = data::Rates::opt;

            data::Wells dw;
            for( const auto& itr : this->wellMap_ ) {
                const auto well_index = itr.second[ 0 ];
                if (!this->open_for_output_[well_index])
                    continue;

                auto& well = dw[ itr.first ];
                well.bhp = this->bhp().at( well_index );
                well.thp = this->thp().at( well_index );
                well.temperature = this->temperature().at( well_index );

                const auto wellrate_index = well_index * pu.num_phases;
                const auto& wv = this->wellRates();
                if( pu.phase_used[BlackoilPhases::Aqua] ) {
                    well.rates.set( rt::wat, wv[ wellrate_index + pu.phase_pos[BlackoilPhases::Aqua] ] );
                }

                if( pu.phase_used[BlackoilPhases::Liquid] ) {
                    well.rates.set( rt::oil, wv[ wellrate_index + pu.phase_pos[BlackoilPhases::Liquid] ] );
                }

                if( pu.phase_used[BlackoilPhases::Vapour] ) {
                    well.rates.set( rt::gas, wv[ wellrate_index + pu.phase_pos[BlackoilPhases::Vapour] ] );
                }

                const int num_perf_well = this->well_perf_data_[well_index].size();
                well.connections.resize(num_perf_well);

                for( int i = 0; i < num_perf_well; ++i ) {
                    const auto active_index = this->well_perf_data_[well_index][i].cell_index;
                    auto& connection = well.connections[ i ];
                    connection.index = globalCellIdxMap[active_index];
                    connection.pressure = this->perfPress()[ itr.second[1] + i ];
                    connection.reservoir_rate = this->perfRates()[ itr.second[1] + i ];
                }
                assert(num_perf_well == int(well.connections.size()));
            }

            return dw;

        }

        virtual ~WellState() = default;
        WellState() = default;
        WellState(const WellState& rhs)  = default;
        WellState& operator=(const WellState& rhs) = default;

    private:
        std::vector<double> bhp_;
        std::vector<double> thp_;
        std::vector<double> temperature_;
        std::vector<double> wellrates_;
        std::vector<double> perfrates_;
        std::vector<double> perfpress_;
    protected:
        std::vector<bool>   open_for_output_;
    private:

        WellMapType wellMap_;

        void initSingleWell(const std::vector<double>& cellPressures,
                            const int w,
                            const Well& well,
                            const PhaseUsage& pu,
                            const SummaryState& summary_state)
        {
            assert(well.isInjector() || well.isProducer());

            // Set default zero initial well rates.
            // May be overwritten below.
            const int np = pu.num_phases;
            for (int p = 0; p < np; ++p) {
                wellrates_[np*w + p] = 0.0;
            }

            const int num_perf_this_well = well_perf_data_[w].size();
            if ( num_perf_this_well == 0 ) {
                // No perforations of the well. Initialize to zero.
                bhp_[w] = 0.;
                thp_[w] = 0.;
                return;
            }

            const auto inj_controls = well.isInjector() ? well.injectionControls(summary_state) : Well::InjectionControls(0);
            const auto prod_controls = well.isProducer() ? well.productionControls(summary_state) : Well::ProductionControls(0);

            const bool is_bhp = well.isInjector() ? (inj_controls.cmode == Well::InjectorCMode::BHP)
                : (prod_controls.cmode == Well::ProducerCMode::BHP);
            const double bhp_limit = well.isInjector() ? inj_controls.bhp_limit : prod_controls.bhp_limit;
            const bool is_grup = well.isInjector() ? (inj_controls.cmode == Well::InjectorCMode::GRUP)
                : (prod_controls.cmode == Well::ProducerCMode::GRUP);

            const double inj_surf_rate = well.isInjector() ? inj_controls.surface_rate : 0.0; // To avoid a "maybe-uninitialized" warning.

            if (well.getStatus() == Well::Status::STOP) {
                // Stopped well:
                // 1. Rates: zero well rates.
                // 2. Bhp: assign bhp equal to bhp control, if
                //    applicable, otherwise assign equal to
                //    first perforation cell pressure.
                if (is_bhp) {
                    bhp_[w] = bhp_limit;
                } else {
                    const int first_cell = well_perf_data_[w][0].cell_index;
                    bhp_[w] = cellPressures[first_cell];
                }
            } else if (is_grup) {
                // Well under group control.
                // 1. Rates: zero well rates.
                // 2. Bhp: initialize bhp to be a
                //    little above or below (depending on if
                //    the well is an injector or producer)
                //    pressure in first perforation cell.
                const int first_cell = well_perf_data_[w][0].cell_index;
                const double safety_factor = well.isInjector() ? 1.01 : 0.99;
                bhp_[w] = safety_factor*cellPressures[first_cell];
            } else {
                // Open well, under own control:
                // 1. Rates: initialize well rates to match
                //    controls if type is ORAT/GRAT/WRAT
                //    (producer) or RATE (injector).
                //    Otherwise, we cannot set the correct
                //    value here and initialize to zero rate.
                if (well.isInjector()) {
                    if (inj_controls.cmode == Well::InjectorCMode::RATE) {
                        switch (inj_controls.injector_type) {
                        case Well::InjectorType::WATER:
                            assert(pu.phase_used[BlackoilPhases::Aqua]);
                            wellrates_[np*w + pu.phase_pos[BlackoilPhases::Aqua]] = inj_surf_rate;
                            break;
                        case Well::InjectorType::GAS:
                            assert(pu.phase_used[BlackoilPhases::Vapour]);
                            wellrates_[np*w + pu.phase_pos[BlackoilPhases::Vapour]] = inj_surf_rate;
                            break;
                        case Well::InjectorType::OIL:
                            assert(pu.phase_used[BlackoilPhases::Liquid]);
                            wellrates_[np*w + pu.phase_pos[BlackoilPhases::Liquid]] = inj_surf_rate;
                            break;
                        case Well::InjectorType::MULTI:
                            // Not currently handled, keep zero init.
                            break;
                        }
                    } else {
                        // Keep zero init.
                    }
                } else {
                    assert(well.isProducer());
                    // Note negative rates for producing wells.
                    switch (prod_controls.cmode) {
                    case Well::ProducerCMode::ORAT:
                        assert(pu.phase_used[BlackoilPhases::Liquid]);
                        wellrates_[np*w + pu.phase_pos[BlackoilPhases::Liquid]] = -prod_controls.oil_rate;
                        break;
                    case Well::ProducerCMode::WRAT:
                        assert(pu.phase_used[BlackoilPhases::Aqua]);
                        wellrates_[np*w + pu.phase_pos[BlackoilPhases::Aqua]] = -prod_controls.water_rate;
                        break;
                    case Well::ProducerCMode::GRAT:
                        assert(pu.phase_used[BlackoilPhases::Vapour]);
                        wellrates_[np*w + pu.phase_pos[BlackoilPhases::Vapour]] = -prod_controls.gas_rate;
                        break;
                    default:
                        // Keep zero init.
                        break;
                    }
                }
                // 2. Bhp: initialize bhp to be target pressure if
                //    bhp-controlled well, otherwise set to a
                //    little above or below (depending on if
                //    the well is an injector or producer)
                //    pressure in first perforation cell.
                if (is_bhp) {
                    bhp_[w] = bhp_limit;
                } else {
                    const int first_cell = well_perf_data_[w][0].cell_index;
                    const double safety_factor = well.isInjector() ? 1.01 : 0.99;
                    bhp_[w] = safety_factor*cellPressures[first_cell];
                }
            }

            // 3. Thp: assign thp equal to thp target/limit, if such a limit exists,
            //    otherwise keep it zero.
            const bool has_thp = well.isInjector() ? inj_controls.hasControl(Well::InjectorCMode::THP)
                : prod_controls.hasControl(Well::ProducerCMode::THP);
            const double thp_limit = well.isInjector() ? inj_controls.thp_limit : prod_controls.thp_limit;
            if (has_thp) {
                thp_[w] = thp_limit;
            }

        }

    protected:
        std::vector<std::vector<PerforationData>> well_perf_data_;
    };

} // namespace Ewoms

#endif // EWOMS_WELLSTATE_HH
