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

#ifndef EWOMS_WELLCONNECTIONAUXILIARYMODULE_HH
#define EWOMS_WELLCONNECTIONAUXILIARYMODULE_HH

#include <ewoms/numerics/discretizations/common/baseauxiliarymodule.hh>

#include <ewoms/eclgrids/cpgrid.hh>

#include <ewoms/eclio/parser/eclipsestate/eclipsestate.hh>

namespace Ewoms
{
template<class TypeTag>
class WellConnectionAuxiliaryModule
    : public Ewoms::BaseAuxiliaryModule<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;

public:

    using NeighborSet = typename
        Ewoms::BaseAuxiliaryModule<TypeTag>::NeighborSet;

    WellConnectionAuxiliaryModule(const Schedule& schedule,
                                  const Dune::CpGrid& grid)
    {
        // Create cartesian to compressed mapping
        const auto& globalCell = grid.globalCell();
        const auto& cartesianSize = grid.logicalCartesianSize();

        auto size = cartesianSize[0]*cartesianSize[1]*cartesianSize[2];

        std::vector<int> cartesianToCompressed(size, -1);
        auto begin = globalCell.begin();

        for ( auto cell = begin, end= globalCell.end(); cell != end; ++cell )
        {
          cartesianToCompressed[ *cell ] = cell - begin;
        }

        const auto& schedule_wells = schedule.getWellsatEnd();
        wells_.reserve(schedule_wells.size());

        // initialize the additional cell connections introduced by wells.
        for ( const auto well : schedule_wells )
        {
            std::vector<int> compressed_well_perforations;
            // All possible completions of the well
            const auto& completionSet = well.getConnections();
            compressed_well_perforations.reserve(completionSet.size());

            for ( size_t c=0; c < completionSet.size(); c++ )
            {
                const auto& completion = completionSet.get(c);
                int i = completion.getI();
                int j = completion.getJ();
                int k = completion.getK();
                int cart_grid_idx = i + cartesianSize[0]*(j + cartesianSize[1]*k);
                int compressed_idx = cartesianToCompressed[cart_grid_idx];

                if ( compressed_idx >= 0 ) // Ignore completions in inactive/remote cells.
                {
                    compressed_well_perforations.push_back(compressed_idx);
                }
            }

            if ( ! compressed_well_perforations.empty() )
            {
                std::sort(compressed_well_perforations.begin(),
                          compressed_well_perforations.end());

                wells_.push_back(compressed_well_perforations);
            }
        }
    }

    unsigned numDofs() const
    {
        // No extra dofs are inserted for wells.
        return 0;
    }

    void addNeighbors(std::vector<NeighborSet>& neighbors) const
    {
        for(const auto well_perforations : wells_)
        {
            for(const auto& perforation : well_perforations)
                neighbors[perforation].insert(well_perforations.begin(),
                                              well_perforations.end());
        }
    }

    void applyInitial()
    {}

    void linearize(SparseMatrixAdapter& , GlobalEqVector&)
    {
        // Linearization is done in StandardDenseWells
    }

    private:

    std::vector<std::vector<int> > wells_;
};

} // end namespace Ewoms
#endif
