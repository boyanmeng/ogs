/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include "ProcessLib/ThermalTwoPhaseFlowComponential/ThermalTwoPhaseFlowComponentialMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowComponential
{
std::unique_ptr<ThermalTwoPhaseFlowComponentialMaterialProperties>
createThermalTwoPhaseFlowComponentialMaterialProperties();

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
