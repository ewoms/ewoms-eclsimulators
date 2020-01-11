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

#ifndef EWOMS_BLACKOILPHASES_HH
#define EWOMS_BLACKOILPHASES_HH

namespace Ewoms
{

    class BlackoilPhases
    {
    public:
        static const int MaxNumPhases = 3;

        // "Crypto phases" are "phases" (or rather "conservation quantities") in the
        // sense that they can be active or not and canonical indices can be translated
        // to and from active ones. That said, they are not considered by num_phases or
        // MaxNumPhases. The crypto phases which are currently implemented are solvent,
        // polymer, energy, polymer molecular weight, foam and brine.
        static const int NumCryptoPhases = 6;

        // enum ComponentIndex { Water = 0, Oil = 1, Gas = 2 };
        enum PhaseIndex { Aqua = 0, Liquid = 1, Vapour = 2, Solvent = 3, Polymer = 4, Energy = 5, PolymerMW = 6, Foam = 7, Brine = 8 };
    };

    struct PhaseUsage : public BlackoilPhases
    {
        int num_phases;
        int phase_used[MaxNumPhases + NumCryptoPhases];
        int phase_pos[MaxNumPhases + NumCryptoPhases];
        bool has_solvent;
        bool has_polymer;
        bool has_energy;
        // polymer molecular weight
        bool has_polymermw;
        bool has_foam;
	bool has_brine;
    };

    /// Check or assign presence of a formed, free phase.  Limited to
    /// the 'BlackoilPhases'.
    ///
    /// Use a std::vector<PhasePresence> to represent the conditions
    /// in an entire model.
    class PhasePresence
    {
    public:
        PhasePresence()
            : present_(0)
        {}

        bool hasFreeWater() const { return present(BlackoilPhases::Aqua  ); }
        bool hasFreeOil  () const { return present(BlackoilPhases::Liquid); }
        bool hasFreeGas  () const { return present(BlackoilPhases::Vapour); }

        void setFreeWater() { insert(BlackoilPhases::Aqua  ); }
        void setFreeOil  () { insert(BlackoilPhases::Liquid); }
        void setFreeGas  () { insert(BlackoilPhases::Vapour); }

        bool operator==(const PhasePresence& other) const { return present_ == other.present_; }
        bool operator!=(const PhasePresence& other) const { return !this->operator==(other); }

    private:
        unsigned char present_;

        bool present(const BlackoilPhases::PhaseIndex i) const
        {
            return present_ & (1 << i);
        }

        void insert(const BlackoilPhases::PhaseIndex i)
        {
            present_ |= (1 << i);
        }
    };

} // namespace Ewoms

#endif // EWOMS_BLACKOILPHASES_HH
