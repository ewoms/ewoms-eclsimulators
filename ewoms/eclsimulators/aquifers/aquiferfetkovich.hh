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

#ifndef EWOMS_AQUIFETP_HH
#define EWOMS_AQUIFETP_HH

#include <ewoms/eclsimulators/aquifers/aquiferinterface.hh>

#include <ewoms/eclio/output/data/aquifer.hh>

#include <exception>
#include <stdexcept>

namespace Ewoms
{

template <typename TypeTag>
class AquiferFetkovich : public AquiferInterface<TypeTag>
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

    AquiferFetkovich(const std::vector<Aquancon::AquancCell>& connections,
                     const std::unordered_map<int, int>& cartesian_to_compressed,
                     const Simulator& eebosSimulator,
                     const Aquifetp::AQUFETP_data& aqufetp_data)
        : Base(aqufetp_data.aquiferID, connections, cartesian_to_compressed, eebosSimulator)
        , aqufetp_data_(aqufetp_data)
    {
    }

    void endTimeStep() override
    {
        for (const auto& q : this->Qai_) {
            this->W_flux_ += q * this->eebos_simulator_.timeStepSize();
            aquifer_pressure_ = aquiferPressure();
        }
    }

protected:
    // Aquifer Fetkovich Specific Variables
    // TODO: using const reference here will cause segmentation fault, which is very strange
    const Aquifetp::AQUFETP_data aqufetp_data_;
    Scalar aquifer_pressure_; // aquifer

    inline void initializeConnections() override
    {
        const auto& eclState = this->eebos_simulator_.vanguard().eclState();
        const auto& ugrid = this->eebos_simulator_.vanguard().grid();
        const auto& grid = eclState.getInputGrid();

        auto globalCellIdx = ugrid.globalCell();

        // We hack the cell depth values for now. We can actually get it from elementcontext pos
        this->cell_depth_.resize(this->size(), aqufetp_data_.d0);
        this->alphai_.resize(this->size(), 1.0);
        this->faceArea_connected_.resize(this->size(), 0.0);

        auto cell2Faces = Ewoms::UgGridHelpers::cell2Faces(ugrid);
        auto faceCells = Ewoms::UgGridHelpers::faceCells(ugrid);

        // Translate the C face tag into the enum used by ewoms-eclio's TransMult class
        Ewoms::FaceDir::DirEnum faceDirection;

        // denom_face_areas is the sum of the areas connected to an aquifer
        Scalar denom_face_areas = 0.;
        this->cellToConnectionIdx_.resize(this->eebos_simulator_.gridView().size(/*codim=*/0), -1);
        for (size_t idx = 0; idx < this->size(); ++idx) {
            const auto global_index = this->connections_[idx].global_index;
            const int cell_index = this->cartesian_to_compressed_.at(global_index);

            this->cellToConnectionIdx_[cell_index] = idx;
            const auto cellCenter = grid.getCellCenter(global_index);
            this->cell_depth_.at(idx) = cellCenter[2];

            if (!this->connections_[idx].influx_coeff.first) { // influx_coeff is defaulted
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
                                  "Initialization of Aquifer problem. Make sure faceTag is correctly defined");
                    }

                    if (faceDirection == this->connections_[idx].face_dir) {
                        this->faceArea_connected_[idx] = this->getFaceArea(faceCells, ugrid, faceIdx, idx);
                        break;
                    }
                }
            } else {
                this->faceArea_connected_.at(idx) = this->connections_[idx].influx_coeff.second;
            }
            denom_face_areas += (this->connections_[idx].influx_mult * this->faceArea_connected_.at(idx));
        }

        const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
        for (size_t idx = 0; idx < this->size(); ++idx) {
            this->alphai_.at(idx) = (denom_face_areas < eps_sqrt)
                ? // Prevent no connection NaNs due to division by zero
                0.
                : (this->connections_[idx].influx_mult * this->faceArea_connected_.at(idx)) / denom_face_areas;
        }
    }

    void assignRestartData(const data::AquiferData& xaq) override
    {
        if (xaq.type != data::AquiferType::Fetkovich) {
            throw std::invalid_argument {"Analytic aquifer data for unexpected aquifer type "
                                         "passed to Fetkovich aquifer"};
        }

        this->aquifer_pressure_ = xaq.pressure;
    }

    inline Eval dpai(int idx)
    {
        const Eval dp = aquifer_pressure_ - this->pressure_current_.at(idx)
            + this->rhow_[idx] * this->gravity_() * (this->cell_depth_[idx] - aqufetp_data_.d0);
        return dp;
    }

    // This function implements Eq 5.12 of the EclipseTechnicalDescription
    inline Scalar aquiferPressure()
    {
        Scalar Flux = this->W_flux_.value();
        Scalar pa_ = this->pa0_ - Flux / (aqufetp_data_.C_t * aqufetp_data_.V0);
        return pa_;
    }

    inline void calculateAquiferConstants() override
    {
        this->Tc_ = (aqufetp_data_.C_t * aqufetp_data_.V0) / aqufetp_data_.J;
    }
    // This function implements Eq 5.14 of the EclipseTechnicalDescription
    inline void calculateInflowRate(int idx, const Simulator& simulator) override
    {
        const Scalar td_Tc_ = simulator.timeStepSize() / this->Tc_;
        const Scalar coef = (1 - exp(-td_Tc_)) / td_Tc_;
        this->Qai_.at(idx) = this->alphai_[idx] * aqufetp_data_.J * dpai(idx) * coef;
    }

    inline void calculateAquiferCondition() override
    {
        this->rhow_.resize(this->size(), 0.);

        if (this->solution_set_from_restart_) {
            return;
        }

        if (!aqufetp_data_.p0.first) {
            this->pa0_ = calculateReservoirEquilibrium();
        } else {
            this->pa0_ = aqufetp_data_.p0.second;
        }
        aquifer_pressure_ = this->pa0_;
    }

    inline Scalar calculateReservoirEquilibrium() override
    {
        // Since the global_indices are the reservoir index, we just need to extract the fluidstate at those indices
        std::vector<Scalar> pw_aquifer;
        Scalar water_pressure_reservoir;

        ElementContext elemCtx(this->eebos_simulator_);
        const auto& gridView = this->eebos_simulator_.gridView();
        auto elemIt = gridView.template begin</*codim=*/0>();
        const auto& elemEndIt = gridView.template end</*codim=*/0>();
        for (; elemIt != elemEndIt; ++elemIt) {
            const auto& elem = *elemIt;
            elemCtx.updatePrimaryStencil(elem);
            size_t cellIdx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
            int idx = this->cellToConnectionIdx_[cellIdx];
            if (idx < 0)
                continue;

            elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            const auto& iq0 = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
            const auto& fs = iq0.fluidState();

            water_pressure_reservoir = fs.pressure(waterPhaseIdx).value();
            this->rhow_[idx] = fs.density(waterPhaseIdx);
            pw_aquifer.push_back(
                (water_pressure_reservoir
                 - this->rhow_[idx].value() * this->gravity_() * (this->cell_depth_[idx] - aqufetp_data_.d0))
                * this->alphai_[idx]);
        }

        // We take the average of the calculated equilibrium pressures.
        const Scalar sum_alpha = std::accumulate(this->alphai_.begin(), this->alphai_.end(), 0.);
        const Scalar aquifer_pres_avg = std::accumulate(pw_aquifer.begin(), pw_aquifer.end(), 0.) / sum_alpha;
        return aquifer_pres_avg;
    }
}; // Class AquiferFetkovich
} // namespace Ewoms
#endif
