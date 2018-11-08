/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/math/constants/constants.hpp>

namespace BaseLib
{
class ConfigTree;
}
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct Pipe
{
    double const diameter;
    double const wall_thickness;
    double const wall_thermal_conductivity;

    double area() const
    {
        constexpr double pi = boost::math::constants::pi<double>();
        return pi * diameter * diameter / 4;
    }
};

Pipe createPipe(BaseLib::ConfigTree const& config);

}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
