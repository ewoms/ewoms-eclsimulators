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
/*!
 * \file
 * \copydoc Ewoms::VtkEclTracerModule
 */
#ifndef EWOMS_VTK_ECL_TRACER_MODULE_HH
#define EWOMS_VTK_ECL_TRACER_MODULE_HH

#include <ewoms/numerics/io/vtkmultiwriter.hh>
#include <ewoms/numerics/io/baseoutputmodule.hh>

#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>
#include <ewoms/numerics/models/blackoil/blackoilproperties.hh>

#include <dune/common/fvector.hh>

#include <cstdio>

BEGIN_PROPERTIES

// create new type tag for the VTK tracer output
NEW_TYPE_TAG(VtkEclTracer);

// create the property tags needed for the tracer model
NEW_PROP_TAG(EnableVtkOutput);
NEW_PROP_TAG(VtkOutputFormat);
NEW_PROP_TAG(VtkWriteEclTracerConcentration);

// set default values for what quantities to output
SET_BOOL_PROP(VtkEclTracer, VtkWriteEclTracerConcentration, false);

END_PROPERTIES

namespace Ewoms {
    /*!
 * \ingroup Vtk
 *
 * \brief VTK output module for the tracer model's parameters.
 */
    template <class TypeTag>
    class VtkEclTracerModule : public BaseOutputModule<TypeTag>
    {
        typedef BaseOutputModule<TypeTag> ParentType;

        using Simulator = GET_PROP_TYPE(TypeTag, Simulator);
        using ElementContext = GET_PROP_TYPE(TypeTag, ElementContext);

        using GridView = GET_PROP_TYPE(TypeTag, GridView);

        static const int vtkFormat = GET_PROP_VALUE(TypeTag, VtkOutputFormat);
        typedef Ewoms::VtkMultiWriter<GridView, vtkFormat> VtkMultiWriter;

        typedef typename ParentType::ScalarBuffer ScalarBuffer;

    public:
        VtkEclTracerModule(const Simulator& simulator)
            : ParentType(simulator)
        { }

        /*!
     * \brief Register all run-time parameters for the tracer VTK output
     * module.
     */
        static void registerParameters()
        {
            EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteEclTracerConcentration,
                                 "Include the tracer concentration "
                                 "in the VTK output files");
        }

        /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
        void allocBuffers()
        {
            if (eclTracerConcentrationOutput_()){
                const auto& tracerModel = this->simulator_.problem().tracerModel();
                eclTracerConcentration_.resize(tracerModel.numTracers());
                for(size_t tracerIdx=0; tracerIdx<eclTracerConcentration_.size();++tracerIdx){

                    this->resizeScalarBuffer_(eclTracerConcentration_[tracerIdx]);
                }
            }

        }

        /*!
     * \brief Modify the internal buffers according to the intensive quantities relevant for
     *        an element
     */
        void processElement(const ElementContext& elemCtx)
        {
            if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
                return;

            const auto& tracerModel = elemCtx.problem().tracerModel();

            for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
                unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

                if (eclTracerConcentrationOutput_()){
                    for(size_t tracerIdx=0; tracerIdx<eclTracerConcentration_.size();++tracerIdx){
                        eclTracerConcentration_[tracerIdx][globalDofIdx] = tracerModel.tracerConcentration(tracerIdx, globalDofIdx);
                    }
                }
            }
        }

        /*!
     * \brief Add all buffers to the VTK output writer.
     */
        void commitBuffers(BaseOutputWriter& baseWriter)
        {
            VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
            if (!vtkWriter)
                return;

            if (eclTracerConcentrationOutput_()){
                const auto& tracerModel = this->simulator_.problem().tracerModel();
                for(size_t tracerIdx=0; tracerIdx<eclTracerConcentration_.size();++tracerIdx){
                    const std::string tmp = "tracerConcentration_" + tracerModel.tracerName(tracerIdx);
                    this->commitScalarBuffer_(baseWriter,tmp.c_str(), eclTracerConcentration_[tracerIdx]);
                }
            }

        }

    private:
        static bool eclTracerConcentrationOutput_()
        {
            static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteEclTracerConcentration);
            return val;
        }

        std::vector<ScalarBuffer> eclTracerConcentration_;
    };
} // namespace Ewoms

#endif
