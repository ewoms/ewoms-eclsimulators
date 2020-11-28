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
/*!
 * \file
 * \copydoc Ewoms::EclPolyhedralGridVanguard
 */
#ifndef EWOMS_ECL_POLYHEDRAL_GRID_VANGUARD_HH
#define EWOMS_ECL_POLYHEDRAL_GRID_VANGUARD_HH

#include "eclbasevanguard.hh"
#include "ecltransmissibility.hh"

#include <ewoms/eclgrids/polyhedralgrid.hh>

namespace Ewoms {
template <class TypeTag>
class EclPolyhedralGridVanguard;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclPolyhedralGridVanguard, INHERITS_FROM(EclBaseVanguard));

// declare the properties
SET_TYPE_PROP(EclPolyhedralGridVanguard, Vanguard, Ewoms::EclPolyhedralGridVanguard<TypeTag>);
SET_TYPE_PROP(EclPolyhedralGridVanguard, Grid, Dune::PolyhedralGrid<3, 3>);
SET_TYPE_PROP(EclPolyhedralGridVanguard, EquilGrid, GET_PROP_TYPE(TypeTag, Grid));

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 *
 * This class uses Dune::PolyhedralGrid as the simulation grid.
 */
template <class TypeTag>
class EclPolyhedralGridVanguard : public EclBaseVanguard<TypeTag>
{
    friend class EclBaseVanguard<TypeTag>;

    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);

public:
    using Grid = GET_PROP_TYPE(TypeTag, Grid);
    using EquilGrid = GET_PROP_TYPE(TypeTag, EquilGrid);
    using GridView = GET_PROP_TYPE(TypeTag, GridView);

private:
    typedef Grid* GridPointer;
    typedef Dune::CartesianIndexMapper<Grid> CartesianIndexMapper;
    typedef CartesianIndexMapper* CartesianIndexMapperPointer;

public:
    EclPolyhedralGridVanguard(Simulator& simulator)
        : EclBaseVanguard<TypeTag>(simulator),
          simulator_( simulator )
    {
        this->callImplementationInit();
    }

    ~EclPolyhedralGridVanguard()
    {
        delete cartesianIndexMapper_;
        delete grid_;
    }

    /*!
     * \brief Return a reference to the simulation grid.
     */
    Grid& grid()
    { return *grid_; }

    /*!
     * \brief Return a reference to the simulation grid.
     */
    const Grid& grid() const
    { return *grid_; }

    /*!
     * \brief Returns a refefence to the grid which should be used by the EQUIL
     *        initialization code.
     *
     * The EQUIL keyword is used to specify the initial condition of the reservoir in
     * hydrostatic equilibrium. Since the code which does this accepts only Ewoms::CpGrid,
     * this is not necessarily the
     * same as the grid which is used for the actual simulation.
     */
    const EquilGrid& equilGrid() const
    { return *grid_; }

    /*!
     * \brief Indicates that the initial condition has been computed and the memory used
     *        by the EQUIL grid can be released.
     *
     * Depending on the implementation, subsequent accesses to the EQUIL grid lead to
     * crashes.
     */
    void releaseEquilGrid()
    { /* do nothing: The EQUIL grid is the simulation grid! */ }

    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
    void loadBalance()
    { /* do nothing: PolyhedralGrid is not parallel! */ }

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the logically Cartesian index.
     */
    const CartesianIndexMapper& cartesianIndexMapper() const
    { return *cartesianIndexMapper_; }

    /*!
     * \brief Returns mapper from compressed to cartesian indices for the EQUIL grid
     *
     * Since PolyhedralGrid is not parallel, that's always the same as
     * cartesianIndexMapper().
     */
    const CartesianIndexMapper& equilCartesianIndexMapper() const
    { return *cartesianIndexMapper_; }

    /*!
     * \brief Free the memory occupied by the global transmissibility object.
     *
     * After writing the initial solution, this array should not be necessary anymore.
     */
    void releaseGlobalTransmissibilities()
    {
    }

    std::unordered_set<std::string> defunctWellNames() const
    { return defunctWellNames_; }

    const EclTransmissibility<TypeTag>& globalTransmissibility() const
    {
        return simulator_.problem().eclTransmissibilities();
    }

protected:
    void createGrids_()
    {
        grid_ = new Grid(this->eclState().getInputGrid(), this->eclState().fieldProps().porv(true));
        cartesianIndexMapper_ = new CartesianIndexMapper(*grid_);
    }

    void filterConnections_()
    {
        // not handling the removal of completions for this type of grid yet.
    }

    Simulator& simulator_;

    GridPointer grid_;
    CartesianIndexMapperPointer cartesianIndexMapper_;

    std::unordered_set<std::string> defunctWellNames_;
};

} // namespace Ewoms

#endif
