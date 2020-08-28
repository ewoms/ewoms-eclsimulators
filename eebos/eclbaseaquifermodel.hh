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
 * \copydoc Ewoms::EclBaseAquiferModel
 */
#ifndef EWOMS_ECL_BASE_AQUIFER_MODEL_HH
#define EWOMS_ECL_BASE_AQUIFER_MODEL_HH

#include <ewoms/eclio/output/data/aquifer.hh>

#include <ewoms/common/propertysystem.hh>

#include <exception>
#include <stdexcept>
#include <vector>

BEGIN_PROPERTIES

NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(RateVector);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup EclBaseAquiferModel
 *
 * \brief The base class which specifies the API of aquifer models.
 *
 * This class only provides the API for the actual aquifer model, it does not do
 * anything on its own.
 */
template <class TypeTag>
class EclBaseAquiferModel
{
    using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
    using RateVector = GET_PROP_TYPE(TypeTag, RateVector);

public:
    EclBaseAquiferModel(Simulator& simulator)
        : simulator_(simulator)
    {}

    /*!
     * \brief Called once the problem has been fully initialized and the initial
     *        condition has been applied.
     */
    void initialSolutionApplied()
    { }

    /*!
     * \brief Called if aquifers are being initialized from values retrieved
     *        from a restart file.
     *
     * \param[in] aquiferSoln Set of aquifer-related initial values, mostly
     *        pertaining to analytic aquifers.  Contains at minimum the
     *        aquifer pressure and the base run's total produced liquid
     *        volume from the model's aquifers.
     */
    void initFromRestart(const std::vector<data::AquiferData>& aquiferSoln EWOMS_UNUSED)
    {
        throw std::logic_error {
            "Initialization from restart data not supported "
            "for base aquifer model"
        };
    }

    /*!
     * \brief This method is called when a new episode (report step) starts.
     */
    void beginEpisode()
    { }

    /*!
     * \brief This method is called when a new time step (substep) starts.
     */
    void beginTimeStep()
    { }

    /*!
     * \brief This method is called before each Newton-Raphson iteration.
     */
    void beginIteration()
    { }

    /*!
     * \brief Add the water which enters or leaves the reservoir due to aquifiers.
     */
    template <class Context>
    void addToSource(RateVector& rate EWOMS_UNUSED,
                     const Context& context EWOMS_UNUSED,
                     unsigned spaceIdx EWOMS_UNUSED,
                     unsigned timeIdx EWOMS_UNUSED) const
    { }

    /*!
     * \brief This method is called after each Newton-Raphson successful iteration.
     *
     * I.e., no exceptions were thrown during the linearization and linear solution
     * procedures.
     */
    void endIteration()
    { }

    /*!
     * \brief This method is called after each successful time step (substep).
     *
     * I.e., all iterations of the time step were successful and the Newton-Raphson
     * algorithm converged.
     */
    void endTimeStep()
    { }

    /*!
     * \brief This method is called once an episode (report step) has been finished
     *        successfully.
     */
    void endEpisode()
    { }

    /*!
     * \brief Write the internal state of the aquifer model to disk using an ad-hoc file
     *        format.
     */
    template <class Restarter>
    void serialize(Restarter& res EWOMS_UNUSED)
    { }

    /*!
     * \brief Load the internal state of the aquifer model to disk using an ad-hoc file
     *        format.
     */
    template <class Restarter>
    void deserialize(Restarter& res EWOMS_UNUSED)
    { }

protected:
    Simulator& simulator_;
};

} // namespace Ewoms

#endif
