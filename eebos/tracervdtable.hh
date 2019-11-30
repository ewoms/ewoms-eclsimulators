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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::TracerVdTable
 */
#ifndef EWOMS_TRACER_VD_TABLE_HH
#define EWOMS_TRACER_VD_TABLE_HH

#include <ewoms/eclio/parser/eclipsestate/tables/simpletable.hh>

namespace Ewoms {

/*!
 * \brief A class that contains tracer concentration vs depth table
 */
class TracerVdTable : public Ewoms::SimpleTable
{
public:
    TracerVdTable(const Ewoms::DeckItem& item)
    {
        this->m_schema.addColumn(Ewoms::ColumnSchema("DEPTH", Ewoms::Table::STRICTLY_INCREASING, Ewoms::Table::DEFAULT_NONE));
        this->m_schema.addColumn(Ewoms::ColumnSchema("TRACER_CONCENTRATION", Ewoms::Table::RANDOM, Ewoms::Table::DEFAULT_NONE));

        Ewoms::SimpleTable::init(item);
    }

    /*!
     * \brief Return the depth column
     */
    const Ewoms::TableColumn& getDepthColumn() const
    { return Ewoms::SimpleTable::getColumn(0); }

    /*!
     * \brief Return the tracer concentration column
     */
    const Ewoms::TableColumn& getTracerConcentration() const
    { return Ewoms::SimpleTable::getColumn(1); }
};
}

#endif // TRACERVDTABLE_HH

