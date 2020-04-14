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
#include "ParameterLib/Utils.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/ThermalTwoPhaseFlowComponential/CreateThermalTwoPhaseFlowComponentialMaterialProperties.h"
#include "ProcessLib/ThermalTwoPhaseFlowComponential/ThermalTwoPhaseFlowComponentialMaterialProperties.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "ThermalTwoPhaseFlowComponentialProcess.h"
#include "ThermalTwoPhaseFlowComponentialProcessData.h"

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
    config.checkConfigParameter("type", "THERMALTWOPHASEFLOW_COMPONENTIAL");

    DBUG("Create ThermalTwoPhaseFlowComponentialProcess model.");
    //! \ogs_file_param{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    auto per_process_variables = findProcessVariables(
        variables, pv_config,
        {//! \ogs_file_param_special{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__process_variables__liquid_pressure}
         "liquid_pressure",
         //! \ogs_file_param_special{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__process_variables__overall_mass_density}
         "overall_mass_density"});
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
    // diffusion coeff
    auto& diff_coeff_b = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__diffusion_coeff_component_b}
        "diffusion_coeff_component_b", parameters, 1, &mesh);
    auto& diff_coeff_a = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__diffusion_coeff_component_a}
        "diffusion_coeff_component_a", parameters, 1, &mesh);
    auto& temperature = ParameterLib::findParameter<double>(
        config,
        //! \ogs_file_param_special{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__temperature}
        "temperature", parameters, 1, &mesh);

    //! \ogs_file_param{prj__processes__process__THERMALTWOPHASEFLOW_COMPONENTIAL__material_property}
    auto const& mat_config = config.getConfigSubtree("material_property");

    boost::optional<MeshLib::PropertyVector<int> const&> material_ids;
    if (mesh.getProperties().existsPropertyVector<int>("MaterialIDs"))
    {
        INFO("The twophase flow is in heterogeneous porous media.");
        auto const& mat_ids =
            mesh.getProperties().getPropertyVector<int>("MaterialIDs");
        material_ids = *mat_ids;
    }
    else
    {
        INFO("The twophase flow is in homogeneous porous media.");
    }

    std::unique_ptr<ThermalTwoPhaseFlowComponentialMaterialProperties> material =
        createThermalTwoPhaseFlowComponentialMaterialProperties(mat_config, material_ids,
                                                 parameters);

    ThermalTwoPhaseFlowComponentialProcessData process_data{
        specific_body_force, has_gravity, mass_lumping,       diff_coeff_b,
        diff_coeff_a,        temperature, std::move(material)};

    return std::make_unique<ThermalTwoPhaseFlowComponentialProcess>(
        std::move(name), mesh, std::move(jacobian_assembler), parameters,
        integration_order, std::move(process_variables),
        std::move(process_data), std::move(secondary_variables), mat_config,
        curves);
}

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
