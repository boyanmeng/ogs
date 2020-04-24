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

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace ThermalTwoPhaseFlowComponential
{
struct ThermalTwoPhaseFlowComponentialProcessData
{
    Eigen::VectorXd const _specific_body_force;

    bool const _has_gravity;
    bool const _has_mass_lumping;
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
};

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
