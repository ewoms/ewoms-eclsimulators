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

#ifndef EWOMS_AQUIFERINTERFACE_HH
#define EWOMS_AQUIFERINTERFACE_HH

#include <ewoms/eclio/utility/numeric/linearinterpolation.hh>
#include <ewoms/eclio/parser/eclipsestate/aquancon.hh>
#include <ewoms/eclio/parser/eclipsestate/aquiferct.hh>
#include <ewoms/eclio/parser/eclipsestate/aquifetp.hh>

#include <ewoms/eclio/output/data/aquifer.hh>

#include <ewoms/common/mathtoolbox.hh>
#include <ewoms/common/densead/evaluation.hh>
#include <ewoms/common/densead/math.hh>
#include <ewoms/material/fluidstates/blackoilfluidstate.hh>

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace Ewoms
{
template <typename TypeTag>
class AquiferInterface
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) BlackoilIndices;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, IntensiveQuantities) IntensiveQuantities;

    enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableBrine = GET_PROP_VALUE(TypeTag, EnableBrine) };

    static const int numEq = BlackoilIndices::numEq;
    typedef double Scalar;

    typedef DenseAd::Evaluation<double, /*size=*/numEq> Eval;

    typedef Ewoms::BlackOilFluidState<Eval,
                                    FluidSystem,
                                    enableTemperature,
                                    enableEnergy,
                                    BlackoilIndices::gasEnabled,
                                    enableBrine,
                                    BlackoilIndices::numPhases>
        FluidState;

    static const auto waterCompIdx = FluidSystem::waterCompIdx;
    static const auto waterPhaseIdx = FluidSystem::waterPhaseIdx;

    // Constructor
    AquiferInterface(const Aquancon::AquanconOutput& connection,
                     const std::unordered_map<int, int>& cartesian_to_compressed,
                     const Simulator& eebosSimulator)
        : connection_(connection)
        , eebos_simulator_(eebosSimulator)
        , cartesian_to_compressed_(cartesian_to_compressed)
    {
        assert(this->connection_.influx_coeff.size() == this->connection_.global_index.size());
        assert(this->connection_.influx_coeff.size() == this->connection_.influx_multiplier.size());
        assert(this->connection_.influx_multiplier.size() == this->connection_.reservoir_face_dir.size());
    }

    // Deconstructor
    virtual ~AquiferInterface()
    {
    }

    void initFromRestart(const std::vector<data::AquiferData>& aquiferSoln)
    {
        auto xaqPos
            = std::find_if(aquiferSoln.begin(), aquiferSoln.end(), [this](const data::AquiferData& xaq) -> bool {
                  return xaq.aquiferID == this->connection_.aquiferID;
              });

        if (xaqPos == aquiferSoln.end()) {
            // No restart value applies to this aquifer.  Nothing to do.
            return;
        }

        this->assignRestartData(*xaqPos);

        this->W_flux_ = xaqPos->volume;
        this->pa0_ = xaqPos->initPressure;
        this->solution_set_from_restart_ = true;
    }

    void initialSolutionApplied()
    {
        initQuantities();
    }

    void beginTimeStep()
    {
        ElementContext elemCtx(eebos_simulator_);
        auto elemIt = eebos_simulator_.gridView().template begin<0>();
        const auto& elemEndIt = eebos_simulator_.gridView().template end<0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;

            elemCtx.updatePrimaryStencil(elem);

            int cellIdx = elemCtx.globalSpaceIndex(0, 0);
            int idx = cellToConnectionIdx_[cellIdx];
            if (idx < 0)
                continue;

            elemCtx.updateIntensiveQuantities(0);
            const auto& iq = elemCtx.intensiveQuantities(0, 0);
            pressure_previous_[idx] = Ewoms::getValue(iq.fluidState().pressure(waterPhaseIdx));
        }
    }

    template <class Context>
    void addToSource(RateVector& rates, const Context& context, unsigned spaceIdx, unsigned timeIdx)
    {
        unsigned cellIdx = context.globalSpaceIndex(spaceIdx, timeIdx);

        int idx = cellToConnectionIdx_[cellIdx];
        if (idx < 0)
            return;

        // We are dereferencing the value of IntensiveQuantities because cachedIntensiveQuantities return a const
        // pointer to IntensiveQuantities of that particular cell_id
        const IntensiveQuantities intQuants = context.intensiveQuantities(spaceIdx, timeIdx);
        // This is the pressure at td + dt
        updateCellPressure(pressure_current_, idx, intQuants);
        updateCellDensity(idx, intQuants);
        calculateInflowRate(idx, context.simulator());
        rates[BlackoilIndices::conti0EqIdx + FluidSystem::waterCompIdx]
            += Qai_[idx] / context.dofVolume(spaceIdx, timeIdx);
    }

protected:
    inline Scalar gravity_() const
    {
        return eebos_simulator_.problem().gravity()[2];
    }

    inline void initQuantities()
    {
        // We reset the cumulative flux at the start of any simulation, so, W_flux = 0
        if (!this->solution_set_from_restart_) {
            W_flux_ = 0.;
        }

        // We next get our connections to the aquifer and initialize these quantities using the initialize_connections
        // function
        initializeConnections();
        calculateAquiferCondition();
        calculateAquiferConstants();

        pressure_previous_.resize(cell_idx_.size(), 0.);
        pressure_current_.resize(cell_idx_.size(), 0.);
        Qai_.resize(cell_idx_.size(), 0.0);
    }

    inline void
    updateCellPressure(std::vector<Eval>& pressure_water, const int idx, const IntensiveQuantities& intQuants)
    {
        const auto& fs = intQuants.fluidState();
        pressure_water.at(idx) = fs.pressure(waterPhaseIdx);
    }

    inline void
    updateCellPressure(std::vector<Scalar>& pressure_water, const int idx, const IntensiveQuantities& intQuants)
    {
        const auto& fs = intQuants.fluidState();
        pressure_water.at(idx) = fs.pressure(waterPhaseIdx).value();
    }

    inline void updateCellDensity(const int idx, const IntensiveQuantities& intQuants)
    {
        const auto& fs = intQuants.fluidState();
        rhow_.at(idx) = fs.density(waterPhaseIdx);
    }

    template <class faceCellType, class ugridType>
    inline double getFaceArea(const faceCellType& faceCells,
                              const ugridType& ugrid,
                              const int faceIdx,
                              const int idx) const
    {
        // Check now if the face is outside of the reservoir, or if it adjoins an inactive cell
        // Do not make the connection if the product of the two cellIdx > 0. This is because the
        // face is within the reservoir/not connected to boundary. (We still have yet to check for inactive cell
        // adjoining)
        double faceArea = 0.;
        const auto cellNeighbour0 = faceCells(faceIdx, 0);
        const auto cellNeighbour1 = faceCells(faceIdx, 1);
        const auto defaultFaceArea = Ewoms::UgGridHelpers::faceArea(ugrid, faceIdx);
        const auto calculatedFaceArea
            = (!this->connection_.influx_coeff.at(idx)) ? defaultFaceArea : *(this->connection_.influx_coeff.at(idx));
        faceArea = (cellNeighbour0 * cellNeighbour1 > 0) ? 0. : calculatedFaceArea;
        if (cellNeighbour1 == 0) {
            faceArea = (cellNeighbour0 < 0) ? faceArea : 0.;
        } else if (cellNeighbour0 == 0) {
            faceArea = (cellNeighbour1 < 0) ? faceArea : 0.;
        }
        return faceArea;
    }

    virtual void endTimeStep() = 0;

    const Aquancon::AquanconOutput connection_;
    const Simulator& eebos_simulator_;
    const std::unordered_map<int, int> cartesian_to_compressed_;

    // Grid variables
    std::vector<size_t> cell_idx_;
    std::vector<Scalar> faceArea_connected_;
    std::vector<int> cellToConnectionIdx_;
    // Quantities at each grid id
    std::vector<Scalar> cell_depth_;
    std::vector<Scalar> pressure_previous_;
    std::vector<Eval> pressure_current_;
    std::vector<Eval> Qai_;
    std::vector<Eval> rhow_;
    std::vector<Scalar> alphai_;

    Scalar Tc_; // Time constant
    Scalar pa0_; // initial aquifer pressure

    Eval W_flux_;

    bool solution_set_from_restart_ {false};

    virtual void initializeConnections() = 0;

    virtual void assignRestartData(const data::AquiferData& xaq) = 0;

    virtual void calculateInflowRate(int idx, const Simulator& simulator) = 0;

    virtual void calculateAquiferCondition() = 0;

    virtual void calculateAquiferConstants() = 0;

    virtual Scalar calculateReservoirEquilibrium() = 0;
    // This function is used to initialize and calculate the alpha_i for each grid connection to the aquifer
};
} // namespace Ewoms
#endif
