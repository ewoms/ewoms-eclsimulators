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
 * \copydoc Ewoms::EclCpGridVanguard
 */
#ifndef EWOMS_ECL_CP_GRID_VANGUARD_HH
#define EWOMS_ECL_CP_GRID_VANGUARD_HH

#include "eclbasevanguard.hh"
#include "ecltransmissibility.hh"
#include "femcpgridcompat.hh"

#include <ewoms/eclgrids/cpgrid.hh>
#include <ewoms/eclgrids/cpgrid/gridhelpers.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/common/version.hh>

namespace Ewoms {
template <class TypeTag>
class EclCpGridVanguard;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclCpGridVanguard, INHERITS_FROM(EclBaseVanguard));

// declare the properties
SET_TYPE_PROP(EclCpGridVanguard, Vanguard, Ewoms::EclCpGridVanguard<TypeTag>);
SET_TYPE_PROP(EclCpGridVanguard, Grid, Dune::CpGrid);
SET_TYPE_PROP(EclCpGridVanguard, EquilGrid, typename GET_PROP_TYPE(TypeTag, Grid));

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 *
 * This class uses Dune::CpGrid as the simulation grid.
 */
template <class TypeTag>
class EclCpGridVanguard : public EclBaseVanguard<TypeTag>
{
    friend class EclBaseVanguard<TypeTag>;
    typedef EclBaseVanguard<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;

public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, EquilGrid) EquilGrid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

private:
    typedef Dune::CartesianIndexMapper<Grid> CartesianIndexMapper;

public:
    EclCpGridVanguard(Simulator& simulator)
        : EclBaseVanguard<TypeTag>(simulator), mpiRank()
    {
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif
        this->callImplementationInit();
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
    {
        assert(mpiRank == 0);
        return *equilGrid_;
    }

    /*!
     * \brief Indicates that the initial condition has been computed and the memory used
     *        by the EQUIL grid can be released.
     *
     * Depending on the implementation, subsequent accesses to the EQUIL grid lead to
     * crashes.
     */
    void releaseEquilGrid()
    {
        equilGrid_.reset();
        equilCartesianIndexMapper_.reset();
    }

    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
    void loadBalance()
    {
#if HAVE_MPI
        int mpiSize = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

        if (mpiSize > 1) {
            // the CpGrid's loadBalance() method likes to have the transmissibilities as
            // its edge weights. since this is (kind of) a layering violation and
            // transmissibilities are relatively expensive to compute, we only do it if
            // more than a single process is involved in the simulation.
            cartesianIndexMapper_.reset(new CartesianIndexMapper(*grid_));
            if (grid_->size(0))
            {
                globalTrans_.reset(new EclTransmissibility<TypeTag>(*this));
                globalTrans_->update();
            }

            Dune::EdgeWeightMethod edgeWeightsMethod = this->edgeWeightsMethod();

            // convert to transmissibility for faces
            // TODO: grid_->numFaces() is not generic. use grid_->size(1) instead? (might
            // not work)
            const auto& gridView = grid_->leafGridView();
            unsigned numFaces = grid_->numFaces();
            std::vector<double> faceTrans(numFaces, 0.0);
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
            ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());
#else
            ElementMapper elemMapper(this->gridView());
#endif
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++ elemIt) {
                const auto& elem = *elemIt;
                auto isIt = gridView.ibegin(elem);
                const auto& isEndIt = gridView.iend(elem);
                for (; isIt != isEndIt; ++ isIt) {
                    const auto& is = *isIt;
                    if (!is.neighbor())
                        continue;

                    unsigned I = elemMapper.index(is.inside());
                    unsigned J = elemMapper.index(is.outside());

                    // FIXME (?): this is not portable!
                    unsigned faceIdx = is.id();

                    faceTrans[faceIdx] = globalTrans_->transmissibility(I, J);
                }
            }

            //distribute the grid and switch to the distributed view.
            {
                const auto wells = this->schedule().getWellsatEnd();
                defunctWellNames_ = std::get<1>(grid_->loadBalance(edgeWeightsMethod, &wells, faceTrans.data()));
            }
            grid_->switchToDistributedView();

            cartesianIndexMapper_.reset();

            if ( ! equilGrid_ )
            {
                // for processes that do not hold the global grid we filter here using the local grid.
                // If we would filter in filterConnection_ our partition would be empty and the connections of all
                // wells would be removed.
                const auto eclipseGrid = Ewoms::UgGridHelpers::createEclipseGrid(grid(), this->eclState().getInputGrid());
                this->schedule().filterConnections(eclipseGrid);
            }
        }
#endif

        cartesianIndexMapper_.reset(new CartesianIndexMapper(*grid_));

        this->updateGridView_();
    }

    /*!
     * \brief Free the memory occupied by the global transmissibility object.
     *
     * After writing the initial solution, this array should not be necessary anymore.
     */
    void releaseGlobalTransmissibilities()
    {
        globalTrans_.reset();
    }

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the logically Cartesian index.
     */
    const CartesianIndexMapper& cartesianIndexMapper() const
    { return *cartesianIndexMapper_; }

    /*!
     * \brief Returns mapper from compressed to cartesian indices for the EQUIL grid
     */
    const CartesianIndexMapper& equilCartesianIndexMapper() const
    {
        assert(mpiRank == 0);
        assert(equilCartesianIndexMapper_);
        return *equilCartesianIndexMapper_;
    }

    std::unordered_set<std::string> defunctWellNames() const
    { return defunctWellNames_; }

    const EclTransmissibility<TypeTag>& globalTransmissibility() const
    {
        assert( globalTrans_ != nullptr );
        return *globalTrans_;
    }

    void releaseGlobalTransmissibility()
    {
        globalTrans_.reset();
    }

protected:
    void createGrids_()
    {
        const auto& gridProps = this->eclState().get3DProperties();
        const std::vector<double>& porv = gridProps.getDoubleGridProperty("PORV").getData();

        grid_.reset(new Dune::CpGrid());
        grid_->processEclipseFormat(this->eclState().getInputGrid(),
                                    /*isPeriodic=*/false,
                                    /*flipNormals=*/false,
                                    /*clipZ=*/false,
                                    porv,
                                    this->eclState().getInputNNC());

        // we use separate grid objects: one for the calculation of the initial condition
        // via EQUIL and one for the actual simulation. The reason is that the EQUIL code
        // is allergic to distributed grids and the simulation grid is distributed before
        // the initial condition is calculated.
        // After loadbalance grid_ will contain a global and distribute view.
        // equilGrid_being a shallow copy only the global view.
        if (mpiRank == 0)
        {
            equilGrid_.reset(new Dune::CpGrid(*grid_));
            equilCartesianIndexMapper_.reset(new CartesianIndexMapper(*equilGrid_));
        }
    }

    // removing some connection located in inactive grid cells
    void filterConnections_()
    {
        // We only filter if we hold the global grid. Otherwise the filtering
        // is done after load balancing as in the future the other processes
        // will hold an empty partition for the global grid and hence filtering
        // here would remove all well connections.
        if (equilGrid_)
        {
            const auto eclipseGrid = Ewoms::UgGridHelpers::createEclipseGrid(equilGrid(), this->eclState().getInputGrid());
            this->schedule().filterConnections(eclipseGrid);
        }
    }

    std::unique_ptr<Grid> grid_;
    std::unique_ptr<EquilGrid> equilGrid_;
    std::unique_ptr<CartesianIndexMapper> cartesianIndexMapper_;
    std::unique_ptr<CartesianIndexMapper> equilCartesianIndexMapper_;

    std::unique_ptr<EclTransmissibility<TypeTag> > globalTrans_;
    std::unordered_set<std::string> defunctWellNames_;
    int mpiRank;
};

} // namespace Ewoms

#endif
