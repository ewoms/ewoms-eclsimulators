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
#include "eflow/eflow_tag.hh"
#include <ewoms/numerics/models/blackoil/blackoilonephaseindices.hh>

BEGIN_PROPERTIES
NEW_TYPE_TAG(EclEFlowProblemSimple, INHERITS_FROM(EclEFlowProblem));
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(FluidSystem);
//! The indices required by the model
SET_PROP(EclEFlowProblemSimple, Indices)
{
private:
    // it is unfortunately not possible to simply use 'TypeTag' here because this leads
    // to cyclic definitions of some properties. if this happens the compiler error
    // messages unfortunately are *really* confusing and not really helpful.
    typedef TTAG(EclEFlowProblem) BaseTypeTag;
    typedef typename GET_PROP_TYPE(BaseTypeTag, FluidSystem) FluidSystem;

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
SET_PROP(EclEFlowProblemSimple, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum { enableTemperature = GET_PROP_VALUE(TypeTag, EnableTemperature) };
    enum { enableSolvent = GET_PROP_VALUE(TypeTag, EnableSolvent) };
    enum { enableEnergy = GET_PROP_VALUE(TypeTag, EnableEnergy) };
    enum { enableBrine = GET_PROP_VALUE(TypeTag, EnableBrine) };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
    static const bool compositionSwitchEnabled = Indices::gasEnabled;

public:
    // typedef Ewoms::BlackOilFluidSystemSimple<Scalar> type;
    typedef Ewoms::BlackOilFluidState<Evaluation,
                                    FluidSystem,
                                    enableTemperature,
                                    enableEnergy,
                                    compositionSwitchEnabled,
                                    enableBrine,
                                    Indices::numPhases>
        type;
};

// SET_PROP(EclEFlowProblemSimple, FluidSystem)
// {
// private:
//   //typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//   typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//   typedef typename GET_PROP_TYPE(TypeTag, Evaluation) Evaluation;
//   typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

// public:
//   typedef Ewoms::BlackOilFluidSystem<Scalar,Indices> type;
// };
END_PROPERTIES

int
main(int argc, char** argv)
{
    typedef TTAG(EclEFlowProblemSimple) TypeTag;
    return mainEFlow<TypeTag>(argc, argv);
}
