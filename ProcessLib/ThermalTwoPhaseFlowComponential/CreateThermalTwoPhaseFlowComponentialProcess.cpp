/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "CreateThermalTwoPhaseFlowComponentialProcess.h"
#include <cassert>

#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "ThermalTwoPhaseFlowComponentialProcess.h"
#include "ThermalTwoPhaseFlowComponentialProcessData.h"
#include "CreateThermalTwoPhaseFlowComponentialMaterialProperties.h"
#include "ThermalTwoPhaseFlowComponentialMaterialProperties.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowComponential
{
std::unique_ptr<Process> createThermalTwoPhaseFlowComponentialProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves,
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "THERMAL_TWOPHASE_COMPONENTIAL");

    DBUG("Create ThermalTwoPhaseFlowComponentialProcess model.");
    //! \ogs_file_param{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__process_variables__liquid_pressure}
         "liquid_pressure",
         //! \ogs_file_param_special{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__process_variables__overall_molar_fraction_air}
         "overall_molar_fraction_air",
         //! \ogs_file_param_special{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__process_variables__overall_molar_fraction_contaminant}
         "overall_molar_fraction_contaminant",
         //! \ogs_file_param_special{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__process_variables__temperature}
         "temperature"});
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    process_variables.push_back(std::move(per_process_variables));

    SecondaryVariableCollection secondary_variables;

    ProcessLib::createSecondaryVariables(config, secondary_variables);
    // Specific body force
    std::vector<double> const b =
        //! \ogs_file_param{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__specific_body_force}
        config.getConfigParameter<std::vector<double>>("specific_body_force");
    assert(!b.empty() && b.size() < 4);
    Eigen::VectorXd specific_body_force(b.size());
    bool const has_gravity = MathLib::toVector(b).norm() > 0;
    if (has_gravity)
    {
        std::copy_n(b.data(), b.size(), specific_body_force.data());
    }

    //! \ogs_file_param{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__mass_lumping}
    auto const mass_lumping = config.getConfigParameter<bool>("mass_lumping");

    std::unique_ptr<ThermalTwoPhaseFlowComponentialMaterialProperties>
        material = nullptr;

    auto const material_ids = materialIDs(mesh);
    if (material_ids)
    {
        INFO("The twophase flow is in heterogeneous porous media.");
    }
    else
    {
        INFO("The twophase flow is in homogeneous porous media.");
    }

    material = createThermalTwoPhaseFlowComponentialMaterialProperties();

    auto media_map =
        MaterialPropertyLib::createMaterialSpatialDistributionMap(media, mesh);
    double const error_tolerance =
        config.getConfigParameter<double>("error_tolerance");

    ThermalTwoPhaseFlowComponentialProcessData process_data{
        specific_body_force,
        has_gravity,
        mass_lumping,
        error_tolerance, std::move(media_map),
        std::move(material)};

    return std::make_unique<ThermalTwoPhaseFlowComponentialProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables), curves);
}

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
