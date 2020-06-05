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

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ThermalTwoPhaseFlowComponentialMaterialProperties.h"

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace ThermalTwoPhaseFlowComponential
{
struct ThermalTwoPhaseFlowComponentialProcessData
{
    Eigen::VectorXd const specific_body_force;

    bool const has_gravity;
    bool const has_mass_lumping;
    ParameterLib::Parameter<double> const& temperature;
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    std::unique_ptr<ThermalTwoPhaseFlowComponentialMaterialProperties> material;
};

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
