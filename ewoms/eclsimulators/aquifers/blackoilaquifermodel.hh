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
#ifndef EWOMS_BLACKOILAQUIFERMODEL_HH
#define EWOMS_BLACKOILAQUIFERMODEL_HH

#include <eebos/eclbaseaquifermodel.hh>

#include <ewoms/eclio/parser/eclipsestate/aquancon.hh>
#include <ewoms/eclio/parser/eclipsestate/aquiferct.hh>
#include <ewoms/eclio/parser/eclipsestate/aquifetp.hh>

#include <ewoms/eclio/output/data/aquifer.hh>

#include <ewoms/eclsimulators/aquifers/aquifercartertracy.hh>
#include <ewoms/eclsimulators/aquifers/aquiferfetkovich.hh>

#include <ewoms/eclgrids/cpgrid.hh>
#include <ewoms/eclgrids/polyhedralgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <ewoms/common/densead/math.hh>

#include <vector>
#include <type_traits>

namespace Ewoms
{

template<class Grid>
class SupportsFaceTag
    : public std::integral_constant<bool, false>
{};

template<>
class SupportsFaceTag<Dune::CpGrid>
    : public std::integral_constant<bool, true>
{};

template<>
class SupportsFaceTag<Dune::PolyhedralGrid<3, 3>>
    : public std::integral_constant<bool, true>
{};

#if HAVE_DUNE_ALUGRID
template<>
class SupportsFaceTag<Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>>
    : public std::integral_constant<bool, true>
{};
#endif

/// Class for handling the blackoil well model.
template <typename TypeTag>
class BlackoilAquiferModel
{
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);

public:
    explicit BlackoilAquiferModel(Simulator& simulator);

    void initialSolutionApplied();
    void initFromRestart(const std::vector<data::AquiferData>& aquiferSoln);

    void beginEpisode();
    void beginTimeStep();
    void beginIteration();
    // add the water rate due to aquifers to the source term.
    template <class Context>
    void addToSource(RateVector& rates, const Context& context, unsigned spaceIdx, unsigned timeIdx) const;
    void endIteration();
    void endTimeStep();
    void endEpisode();

    Ewoms::data::Aquifers aquiferData() const;

    template <class Restarter>
    void serialize(Restarter& res);

    template <class Restarter>
    void deserialize(Restarter& res);

protected:
    // ---------      Types      ---------
    using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);
    using Scalar = GET_PROP_TYPE(TypeTag, Scalar);

    typedef AquiferCarterTracy<TypeTag> AquiferCarterTracy_object;
    typedef AquiferFetkovich<TypeTag> AquiferFetkovich_object;

    Simulator& simulator_;

    // TODO: probably we can use one variable to store both types of aquifers, because
    // they share the base class
    mutable std::vector<AquiferCarterTracy_object> aquifers_CarterTracy;
    mutable std::vector<AquiferFetkovich_object> aquifers_Fetkovich;

    // This initialization function is used to connect the parser objects with the ones needed by AquiferCarterTracy
    void init();

    bool aquiferActive() const;
    bool aquiferCarterTracyActive() const;
    bool aquiferFetkovichActive() const;
};

} // namespace Ewoms

#include "blackoilaquifermodel_impl.hh"
#endif
