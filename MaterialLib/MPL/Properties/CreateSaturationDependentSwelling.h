/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}
namespace ParameterLib
{
struct CoordinateSystem;
}

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createSaturationDependentSwelling(
    BaseLib::ConfigTree const& config,
    ParameterLib::CoordinateSystem const* const local_coordinate_system);
}  // namespace MaterialPropertyLib
