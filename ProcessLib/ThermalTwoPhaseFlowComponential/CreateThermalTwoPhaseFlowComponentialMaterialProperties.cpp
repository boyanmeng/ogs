/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermalTwoPhaseFlowComponentialMaterialProperties.h"
#include "BaseLib/Algorithm.h"
#include "BaseLib/Logging.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"
#include "ThermalTwoPhaseFlowComponentialMaterialProperties.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowComponential
{
std::unique_ptr<ThermalTwoPhaseFlowComponentialMaterialProperties>
createThermalTwoPhaseFlowComponentialMaterialProperties()
{
    DBUG("Reading material properties of two-phase flow process.");

    std::unique_ptr<MaterialLib::Fluid::WaterVaporProperties> vapor_property =
        std::make_unique<MaterialLib::Fluid::WaterVaporProperties>();

    return std::make_unique<ThermalTwoPhaseFlowComponentialMaterialProperties>(
        std::move(vapor_property));
}

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
