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
#include "config.h"
#include <ewoms/eclsimulators/eflow/main.hh>
#include <ewoms/numerics/models/blackoil/blackoilonephaseindices.hh>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclEFlowProblemSimple, INHERITS_FROM(EclEFlowProblem));
//! The indices required by the model
SET_PROP(EclEFlowProblemSimple, Indices)
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    typedef TTAG(EclEFlowProblem) BaseTypeTag;
    typedef GET_PROP_TYPE(BaseTypeTag, FluidSystem) FluidSystem;

public:
    typedef Ewoms::BlackOilOnePhaseIndices<GET_PROP_VALUE(TypeTag, EnableSolvent),
                                         GET_PROP_VALUE(TypeTag, EnablePolymer),
                                         GET_PROP_VALUE(TypeTag, EnableEnergy),
                                         GET_PROP_VALUE(TypeTag, EnableFoam),
                                         GET_PROP_VALUE(TypeTag, EnableBrine),
                                         /*PVOffset=*/0,
                                         /*enebledCompIdx=*/FluidSystem::waterCompIdx>
        type;
};

END_PROPERTIES

int main(int argc, char** argv)
{
    using TypeTag = TTAG(EclEFlowProblemSimple);
    auto mainObject = Ewoms::Main(argc, argv);
    return mainObject.runStatic<TypeTag>();
}
