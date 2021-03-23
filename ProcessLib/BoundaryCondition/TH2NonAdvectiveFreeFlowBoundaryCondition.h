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

#include "GenericNaturalBoundaryCondition.h"
#include "TH2NonAdvectiveFreeFlowBoundaryConditionLocalAssembler.h"
#include "MeshLib/PropertyVector.h"

namespace ProcessLib
{
using TH2NonAdvectiveFreeFlowBoundaryCondition =
    GenericNaturalBoundaryCondition<
        TH2NonAdvectiveFreeFlowBoundaryConditionData,
        TH2NonAdvectiveFreeFlowBoundaryConditionLocalAssembler>;

std::unique_ptr<TH2NonAdvectiveFreeFlowBoundaryCondition>
createTH2NonAdvectiveFreeFlowBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const global_dim, Process const& process,
    unsigned const shapefunction_order);

}  // namespace ProcessLib
