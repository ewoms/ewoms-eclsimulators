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

#include <ewoms/eclsimulators/wells/vfpinjproperties.hh>
#include <ewoms/eclsimulators/deprecated/props/blackoilphases.hh>
#include <ewoms/eclio/errormacros.hh>

#include <ewoms/eclsimulators/wells/vfphelpers.hh>

namespace Ewoms {

VFPInjProperties::VFPInjProperties() {

}

VFPInjProperties::VFPInjProperties(const VFPInjTable* table){
    m_tables[table->getTableNum()] = table;
}

VFPInjProperties::VFPInjProperties(const VFPInjProperties::InjTable& tables) {
    for (const auto& table : tables) {
        m_tables[table.first] = table.second.get();
    }
}

double VFPInjProperties::bhp(int table_id,
                                 const double& aqua,
                                 const double& liquid,
                                 const double& vapour,
                                 const double& thp_arg) const {
    const VFPInjTable* table = detail::getTable(m_tables, table_id);

    detail::VFPEvaluation retval = detail::bhp(table, aqua, liquid, vapour, thp_arg);
    return retval.value;
}

double VFPInjProperties::thp(int table_id,
                             const double& aqua,
                             const double& liquid,
                             const double& vapour,
                             const double& bhp_arg) const {
    const VFPInjTable* table = detail::getTable(m_tables, table_id);

    //Find interpolation variables
    double flo = detail::getFlo(aqua, liquid, vapour, table->getFloType());

    const std::vector<double> thp_array = table->getTHPAxis();
    int nthp = thp_array.size();

    /**
     * Find the function bhp_array(thp) by creating a 1D view of the data
     * by interpolating for every value of thp. This might be somewhat
     * expensive, but let us assome that nthp is small
     */
    auto flo_i = detail::findInterpData(flo, table->getFloAxis());
    std::vector<double> bhp_array(nthp);
    for (int i=0; i<nthp; ++i) {
        auto thp_i = detail::findInterpData(thp_array[i], thp_array);
        bhp_array[i] = detail::interpolate(*table, flo_i, thp_i).value;
    }

    double retval = detail::findTHP(bhp_array, thp_array, bhp_arg);
    return retval;
}

const VFPInjTable* VFPInjProperties::getTable(const int table_id) const {
    return detail::getTable(m_tables, table_id);
}

bool VFPInjProperties::hasTable(const int table_id) const {
    return detail::hasTable(m_tables, table_id);
}

} // namespace Ewoms
