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

#ifndef EWOMS_ECLSIM_AQUIFERCT_HH
#define EWOMS_ECLSIM_AQUIFERCT_HH

#include <ewoms/eclsimulators/aquifers/aquiferinterface.hh>

#include <ewoms/eclio/output/data/aquifer.hh>

#include <exception>
#include <stdexcept>

namespace Ewoms
{

template <typename TypeTag>
class AquiferCarterTracy : public AquiferInterface<TypeTag>
{
public:
    typedef AquiferInterface<TypeTag> Base;

    using typename Base::BlackoilIndices;
    using typename Base::ElementContext;
    using typename Base::Eval;
    using typename Base::FluidState;
    using typename Base::FluidSystem;
    using typename Base::IntensiveQuantities;
    using typename Base::RateVector;
    using typename Base::Scalar;
    using typename Base::Simulator;

    using Base::waterCompIdx;
    using Base::waterPhaseIdx;
    AquiferCarterTracy(const Aquancon::AquanconOutput& connection,
                       const std::unordered_map<int, int>& cartesian_to_compressed,
                       const Simulator& eebosSimulator,
                       const AquiferCT::AQUCT_data& aquct_data)
        : Base(connection, cartesian_to_compressed, eebosSimulator)
        , aquct_data_(aquct_data)
    {
    }

    void endTimeStep() override
    {
        for (const auto& Qai : Base::Qai_) {
            Base::W_flux_ += Qai * Base::eebos_simulator_.timeStepSize();
        }
    }

protected:
    // Variables constants
    const AquiferCT::AQUCT_data aquct_data_;
    Scalar beta_; // Influx constant
    // TODO: it is possible it should be a AD variable
    Scalar mu_w_; // water viscosity

    // This function is used to initialize and calculate the alpha_i for each grid connection to the aquifer
    inline void initializeConnections() override
    {
        const auto& eclState = Base::eebos_simulator_.vanguard().eclState();
        const auto& ugrid = Base::eebos_simulator_.vanguard().grid();
        const auto& grid = eclState.getInputGrid();

        Base::cell_idx_ = this->connection_.global_index;
        auto globalCellIdx = ugrid.globalCell();

        // We hack the cell depth values for now. We can actually get it from elementcontext pos
        Base::cell_depth_.resize(Base::cell_idx_.size(), aquct_data_.d0);
        Base::alphai_.resize(Base::cell_idx_.size(), 1.0);
        Base::faceArea_connected_.resize(Base::cell_idx_.size(), 0.0);

        auto cell2Faces = Ewoms::UgGridHelpers::cell2Faces(ugrid);
        auto faceCells = Ewoms::UgGridHelpers::faceCells(ugrid);

        // Translate the C face tag into the enum used by ewoms-eclio's TransMult class
        Ewoms::FaceDir::DirEnum faceDirection;

        // denom_face_areas is the sum of the areas connected to an aquifer
        Scalar denom_face_areas = 0.;
        Base::cellToConnectionIdx_.resize(Base::eebos_simulator_.gridView().size(/*codim=*/0), -1);
        for (size_t idx = 0; idx < Base::cell_idx_.size(); ++idx) {
            const int cell_index = Base::cartesian_to_compressed_.at(Base::cell_idx_[idx]);
            Base::cellToConnectionIdx_[cell_index] = idx;

            const auto cellFacesRange = cell2Faces[cell_index];
            for (auto cellFaceIter = cellFacesRange.begin(); cellFaceIter != cellFacesRange.end(); ++cellFaceIter) {
                // The index of the face in the compressed grid
                const int faceIdx = *cellFaceIter;

                // the logically-Cartesian direction of the face
                const int faceTag = Ewoms::UgGridHelpers::faceTag(ugrid, cellFaceIter);

                switch (faceTag) {
                case 0:
                    faceDirection = Ewoms::FaceDir::XMinus;
                    break;
                case 1:
                    faceDirection = Ewoms::FaceDir::XPlus;
                    break;
                case 2:
                    faceDirection = Ewoms::FaceDir::YMinus;
                    break;
                case 3:
                    faceDirection = Ewoms::FaceDir::YPlus;
                    break;
                case 4:
                    faceDirection = Ewoms::FaceDir::ZMinus;
                    break;
                case 5:
                    faceDirection = Ewoms::FaceDir::ZPlus;
                    break;
                default:
                    EWOMS_THROW(Ewoms::NumericalIssue,
                              "Initialization of Aquifer Carter Tracy problem. Make sure faceTag is correctly defined");
                }

                if (faceDirection == this->connection_.reservoir_face_dir.at(idx)) {
                    Base::faceArea_connected_.at(idx) = Base::getFaceArea(faceCells, ugrid, faceIdx, idx);
                    denom_face_areas += (this->connection_.influx_multiplier.at(idx) * Base::faceArea_connected_.at(idx));
                }
            }
            auto cellCenter = grid.getCellCenter(Base::cell_idx_.at(idx));
            Base::cell_depth_.at(idx) = cellCenter[2];
        }

        const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
        for (size_t idx = 0; idx < Base::cell_idx_.size(); ++idx) {
            Base::alphai_.at(idx) = (denom_face_areas < eps_sqrt)
                ? // Prevent no connection NaNs due to division by zero
                0.
                : (this->connection_.influx_multiplier.at(idx) * Base::faceArea_connected_.at(idx)) / denom_face_areas;
        }
    }

    void assignRestartData(const data::AquiferData& /* xaq */) override
    {
        throw std::runtime_error {"Restart-based initialization not currently supported "
                                  "for Carter-Tracey analytic aquifers"};
    }

    inline void getInfluenceTableValues(Scalar& pitd, Scalar& pitd_prime, const Scalar& td)
    {
        // We use the ewoms-common numeric linear interpolator
        pitd = Ewoms::linearInterpolation(aquct_data_.td, aquct_data_.pi, td);
        pitd_prime = Ewoms::linearInterpolationDerivative(aquct_data_.td, aquct_data_.pi, td);
    }

    inline Scalar dpai(int idx)
    {
        Scalar dp = Base::pa0_
            + Base::rhow_.at(idx).value() * Base::gravity_() * (Base::cell_depth_.at(idx) - aquct_data_.d0)
            - Base::pressure_previous_.at(idx);
        return dp;
    }

    // This function implements Eqs 5.8 and 5.9 of the EclipseTechnicalDescription
    inline void calculateEqnConstants(Scalar& a, Scalar& b, const int idx, const Simulator& simulator)
    {
        const Scalar td_plus_dt = (simulator.timeStepSize() + simulator.time()) / Base::Tc_;
        const Scalar td = simulator.time() / Base::Tc_;
        Scalar PItdprime = 0.;
        Scalar PItd = 0.;
        getInfluenceTableValues(PItd, PItdprime, td_plus_dt);
        a = 1.0 / Base::Tc_ * ((beta_ * dpai(idx)) - (Base::W_flux_.value() * PItdprime)) / (PItd - td * PItdprime);
        b = beta_ / (Base::Tc_ * (PItd - td * PItdprime));
    }

    // This function implements Eq 5.7 of the EclipseTechnicalDescription
    inline void calculateInflowRate(int idx, const Simulator& simulator) override
    {
        Scalar a, b;
        calculateEqnConstants(a, b, idx, simulator);
        Base::Qai_.at(idx)
            = Base::alphai_.at(idx) * (a - b * (Base::pressure_current_.at(idx) - Base::pressure_previous_.at(idx)));
    }

    inline void calculateAquiferConstants() override
    {
        // We calculate the influx constant
        beta_ = aquct_data_.c2 * aquct_data_.h * aquct_data_.theta * aquct_data_.phi_aq * aquct_data_.C_t
            * aquct_data_.r_o * aquct_data_.r_o;
        // We calculate the time constant
        Base::Tc_ = mu_w_ * aquct_data_.phi_aq * aquct_data_.C_t * aquct_data_.r_o * aquct_data_.r_o
            / (aquct_data_.k_a * aquct_data_.c1);
    }

    inline void calculateAquiferCondition() override
    {

        int pvttableIdx = aquct_data_.pvttableID - 1;
        Base::rhow_.resize(Base::cell_idx_.size(), 0.);
        if (!aquct_data_.p0.first) {
            Base::pa0_ = calculateReservoirEquilibrium();
        } else {
            Base::pa0_ = aquct_data_.p0.second;
        }

        // use the thermodynamic state of the first active cell as a
        // reference. there might be better ways to do this...
        ElementContext elemCtx(Base::eebos_simulator_);
        auto elemIt = Base::eebos_simulator_.gridView().template begin</*codim=*/0>();
        elemCtx.updatePrimaryStencil(*elemIt);
        elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
        const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
        // Initialize a FluidState object first
        FluidState fs_aquifer;
        // We use the temperature of the first cell connected to the aquifer
        // Here we copy the fluidstate of the first cell, so we do not accidentally mess up the reservoir fs
        fs_aquifer.assign(iq0.fluidState());
        Eval temperature_aq, pa0_mean;
        temperature_aq = fs_aquifer.temperature(0);
        pa0_mean = Base::pa0_;
        Eval mu_w_aquifer = FluidSystem::waterPvt().viscosity(pvttableIdx, temperature_aq, pa0_mean);
        mu_w_ = mu_w_aquifer.value();
    }

    // This function is for calculating the aquifer properties from equilibrium state with the reservoir
    // TODO: this function can be moved to the Inteface class, since it is the same for both Aquifer models
    inline Scalar calculateReservoirEquilibrium() override
    {
        // Since the global_indices are the reservoir index, we just need to extract the fluidstate at those indices
        std::vector<Scalar> pw_aquifer;
        Scalar water_pressure_reservoir;

        ElementContext elemCtx(Base::eebos_simulator_);
        const auto& gridView = Base::eebos_simulator_.gridView();
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            elemCtx.updatePrimaryStencil(elem);

            size_t cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            int idx = Base::cellToConnectionIdx_[cellIdx];
            if (idx < 0)
                continue;

            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq0.fluidState();

            water_pressure_reservoir = fs.pressure(waterPhaseIdx).value();
            Base::rhow_[idx] = fs.density(waterPhaseIdx);
            pw_aquifer.push_back(
                (water_pressure_reservoir
                 - Base::rhow_[idx].value() * Base::gravity_() * (Base::cell_depth_[idx] - aquct_data_.d0))
                * Base::alphai_[idx]);
        }

        // We take the average of the calculated equilibrium pressures.
        Scalar aquifer_pres_avg = std::accumulate(pw_aquifer.begin(), pw_aquifer.end(), 0.) / pw_aquifer.size();
        return aquifer_pres_avg;
    }
}; // class AquiferCarterTracy
} // namespace Ewoms

#endif
