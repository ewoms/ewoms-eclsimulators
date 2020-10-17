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
#ifndef EWOMS_AUTODIFF_VFPPROPERTIES_HH_
#define EWOMS_AUTODIFF_VFPPROPERTIES_HH_

#include <ewoms/eclio/parser/eclipsestate/schedule/vfpinjtable.hh>
#include <ewoms/eclio/parser/eclipsestate/schedule/vfpprodtable.hh>

#include <map>

namespace Ewoms {

/**
 * A thin wrapper class that holds one VFPProdProperties and one
 * VFPInjProperties object.
 */
template<typename VFPInjProp, typename VFPProdProp>
class VFPProperties {
public:
    VFPProperties() {}

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_table  A *single* VFPINJ table or NULL (no table)
     * @param prod_table A *single* VFPPROD table or NULL (no table)
     */
    explicit VFPProperties(const VFPInjTable* inj_table, const VFPProdTable* prod_table)
    {
      if (inj_table != nullptr) {
        m_inj.reset(new VFPInjProp(inj_table));
      }
      if (prod_table != nullptr) {
        m_prod.reset(new VFPProdProp(prod_table));
      }
    }

    /**
     * Constructor
     * Takes *no* ownership of data.
     * @param inj_tables A map of different VFPINJ tables.
     * @param prod_tables A map of different VFPPROD tables.
     */
    VFPProperties(const std::map<int, std::shared_ptr<const VFPInjTable> >& inj_tables,
                  const std::map<int, std::shared_ptr<const VFPProdTable> >& prod_tables)
      : m_inj(new VFPInjProp(inj_tables)), m_prod(new VFPProdProp(prod_tables)) {}

    /**
     * Returns the VFP properties for injection wells
     */
    const VFPInjProp* getInj() const {
        return m_inj.get();
    }

    /**
     * Returns the VFP properties for production wells
     */
    const VFPProdProp* getProd() const {
        return m_prod.get();
    }

private:
    std::shared_ptr<VFPInjProp> m_inj;
    std::shared_ptr<VFPProdProp> m_prod;
};

} //Namespace

#endif /* EWOMS_AUTODIFF_VFPPROPERTIES_HH_ */
