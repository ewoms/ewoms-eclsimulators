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

namespace Ewoms
{

  template<typename TypeTag>
  class AquiferFetkovich: public AquiferInterface<TypeTag>
  {

  public:
    typedef AquiferInterface<TypeTag> Base;

    using typename Base::Simulator;
    using typename Base::ElementContext;
    using typename Base::FluidSystem;
    using typename Base::BlackoilIndices;
    using typename Base::RateVector;
    using typename Base::IntensiveQuantities;
    using typename Base::Eval;
    using typename Base::Scalar;
    using typename Base::FluidState;

    using Base::waterCompIdx;
    using Base::waterPhaseIdx;

    AquiferFetkovich( const Aquancon::AquanconOutput& connection,
                      const std::unordered_map<int, int>& cartesian_to_compressed,
                      const Simulator& ebosSimulator,
                      const Aquifetp::AQUFETP_data& aqufetp_data)
    : Base(connection, cartesian_to_compressed, ebosSimulator)
    , aqufetp_data_(aqufetp_data)
    {}

    void endTimeStep()
    {
      for (const auto& Qai: Base::Qai_) {
        Base::W_flux_ += Qai*Base::ebos_simulator_.timeStepSize();
        aquifer_pressure_ = aquiferPressure();
      }
    }

  protected:
    // Aquifer Fetkovich Specific Variables
    // TODO: using const reference here will cause segmentation fault, which is very strange
    const Aquifetp::AQUFETP_data aqufetp_data_;
    Scalar aquifer_pressure_; // aquifer

    inline void initializeConnections(const Aquancon::AquanconOutput& connection)
    {
      const auto& eclState = Base::ebos_simulator_.vanguard().eclState();
      const auto& ugrid = Base::ebos_simulator_.vanguard().grid();
      const auto& grid = eclState.getInputGrid();

      Base::cell_idx_ = connection.global_index;
      auto globalCellIdx = ugrid.globalCell();

      assert( Base::cell_idx_ == connection.global_index);
      assert( (Base::cell_idx_.size() == connection.influx_coeff.size()) );
      assert( (connection.influx_coeff.size() == connection.influx_multiplier.size()) );
      assert( (connection.influx_multiplier.size() == connection.reservoir_face_dir.size()) );

      // We hack the cell depth values for now. We can actually get it from elementcontext pos
      Base::cell_depth_.resize(Base::cell_idx_.size(), aqufetp_data_.d0);
      Base::alphai_.resize(Base::cell_idx_.size(), 1.0);
      Base::faceArea_connected_.resize(Base::cell_idx_.size(),0.0);

      auto cell2Faces = Ewoms::UgGridHelpers::cell2Faces(ugrid);
      auto faceCells  = Ewoms::UgGridHelpers::faceCells(ugrid);

      // Translate the C face tag into the enum used by ewoms-eclio's TransMult class
      Ewoms::FaceDir::DirEnum faceDirection;

      // denom_face_areas is the sum of the areas connected to an aquifer
      Scalar denom_face_areas = 0.;
      Base::cellToConnectionIdx_.resize(Base::ebos_simulator_.gridView().size(/*codim=*/0), -1);
      for (size_t idx = 0; idx < Base::cell_idx_.size(); ++idx)
      {
        const int cell_index = Base::cartesian_to_compressed_.at(Base::cell_idx_[idx]);
        Base::cellToConnectionIdx_[cell_index] = idx;

        const auto cellFacesRange = cell2Faces[cell_index];
        for(auto cellFaceIter = cellFacesRange.begin(); cellFaceIter != cellFacesRange.end(); ++cellFaceIter)
        {
          // The index of the face in the compressed grid
          const int faceIdx = *cellFaceIter;

          // the logically-Cartesian direction of the face
          const int faceTag = Ewoms::UgGridHelpers::faceTag(ugrid, cellFaceIter);

          switch(faceTag)
          {
            case 0: faceDirection = Ewoms::FaceDir::XMinus;
            break;
            case 1: faceDirection = Ewoms::FaceDir::XPlus;
            break;
            case 2: faceDirection = Ewoms::FaceDir::YMinus;
            break;
            case 3: faceDirection = Ewoms::FaceDir::YPlus;
            break;
            case 4: faceDirection = Ewoms::FaceDir::ZMinus;
            break;
            case 5: faceDirection = Ewoms::FaceDir::ZPlus;
            break;
            default: EWOMS_THROW(Ewoms::NumericalIssue,"Initialization of Aquifer problem. Make sure faceTag is correctly defined");
          }

          if (faceDirection == connection.reservoir_face_dir.at(idx))
          {
            Base::faceArea_connected_.at(idx) =  Base::getFaceArea(faceCells, ugrid, faceIdx, idx, connection);
            denom_face_areas += ( connection.influx_multiplier.at(idx) * Base::faceArea_connected_.at(idx) );
          }
        }
        auto cellCenter = grid.getCellCenter(Base::cell_idx_.at(idx));
        Base::cell_depth_.at(idx) = cellCenter[2];
      }

      const double eps_sqrt = std::sqrt(std::numeric_limits<double>::epsilon());
      for (size_t idx = 0; idx < Base::cell_idx_.size(); ++idx)
      {
        Base::alphai_.at(idx) = (denom_face_areas < eps_sqrt)? // Prevent no connection NaNs due to division by zero
        0.
        : ( connection.influx_multiplier.at(idx) * Base::faceArea_connected_.at(idx) )/denom_face_areas;
      }
    }

    inline Eval dpai(int idx)
    {
      const Eval dp = aquifer_pressure_ - Base::pressure_current_.at(idx)
                      + Base::rhow_[idx] * Base::gravity_()*(Base::cell_depth_[idx] - aqufetp_data_.d0);
      return dp;
    }

    // This function implements Eq 5.12 of the EclipseTechnicalDescription
    inline Scalar aquiferPressure()
    {
      Scalar Flux = Base::W_flux_.value();
      Scalar pa_ = Base::pa0_ - Flux / ( aqufetp_data_.C_t * aqufetp_data_.V0 );
      return pa_;
    }

    inline void calculateAquiferConstants()
    {
      Base::Tc_ = ( aqufetp_data_.C_t * aqufetp_data_.V0 ) / aqufetp_data_.J ;
    }
    // This function implements Eq 5.14 of the EclipseTechnicalDescription
    inline void calculateInflowRate(int idx, const Simulator& simulator)
    {
      const Scalar td_Tc_ = simulator.timeStepSize() / Base::Tc_ ;
      const Scalar coef = (1 - exp(-td_Tc_)) / td_Tc_;
      Base::Qai_.at(idx) = Base::alphai_[idx] * aqufetp_data_.J * dpai(idx) * coef;
    }

    inline void calculateAquiferCondition()
    {
      Base::rhow_.resize(Base::cell_idx_.size(),0.);
      if (!aqufetp_data_.p0)
      {
        Base::pa0_ = calculateReservoirEquilibrium();
      }
      else
      {
        Base::pa0_ = *(aqufetp_data_.p0);
      }
      aquifer_pressure_ = Base::pa0_ ;
    }

    inline Scalar calculateReservoirEquilibrium()
    {
      // Since the global_indices are the reservoir index, we just need to extract the fluidstate at those indices
      std::vector<Scalar> pw_aquifer;
      Scalar water_pressure_reservoir;

      ElementContext elemCtx(Base::ebos_simulator_);
      const auto& gridView = Base::ebos_simulator_.gridView();
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
        pw_aquifer.push_back( (water_pressure_reservoir - Base::rhow_[idx].value()*Base::gravity_()*(Base::cell_depth_[idx] - aqufetp_data_.d0))*Base::alphai_[idx] );
      }

      // We take the average of the calculated equilibrium pressures.
      const Scalar sum_alpha = std::accumulate(this->alphai_.begin(), this->alphai_.end(), 0.);
      const Scalar aquifer_pres_avg = std::accumulate(pw_aquifer.begin(), pw_aquifer.end(), 0.) / sum_alpha;
      return aquifer_pres_avg;
    }
  }; //Class AquiferFetkovich
} // namespace Ewoms
#endif
