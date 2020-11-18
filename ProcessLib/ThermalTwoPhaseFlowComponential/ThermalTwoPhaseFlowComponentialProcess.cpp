/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermalTwoPhaseFlowComponentialProcess.h"

#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ThermalTwoPhaseFlowComponentialLocalAssembler.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowComponential
{
ThermalTwoPhaseFlowComponentialProcess::ThermalTwoPhaseFlowComponentialProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    ThermalTwoPhaseFlowComponentialProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
    /*curves*/)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables)),
      _process_data(std::move(process_data))
{
    DBUG("Create ThermalTwoPhaseFlowComponentialProcess model.");
}

void ThermalTwoPhaseFlowComponentialProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    ProcessLib::createLocalAssemblers<ThermalTwoPhaseFlowComponentialLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    _secondary_variables.addSecondaryVariable(
        "saturation",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermalTwoPhaseFlowComponentialLocalAssemblerInterface::getIntPtSaturation));

    _secondary_variables.addSecondaryVariable(
        "pressure_nonwetting",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &ThermalTwoPhaseFlowComponentialLocalAssemblerInterface::
                             getIntPtNonwettingPressure));
    _secondary_variables.addSecondaryVariable(
        "capillary_pressure",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermalTwoPhaseFlowComponentialLocalAssemblerInterface::
                getIntPtCapillaryPressure));
    _secondary_variables.addSecondaryVariable(
        "liquid_molar_fraction_air",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermalTwoPhaseFlowComponentialLocalAssemblerInterface::
                getIntPtLiquidMolFracAir));
    _secondary_variables.addSecondaryVariable(
        "liquid_molar_fraction_contaminant",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermalTwoPhaseFlowComponentialLocalAssemblerInterface::
                getIntPtLiquidMolFracContaminant));
    _secondary_variables.addSecondaryVariable(
        "gas_molar_fraction_water",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermalTwoPhaseFlowComponentialLocalAssemblerInterface::
                getIntPtGasMolFracWater));
    _secondary_variables.addSecondaryVariable(
        "gas_molar_fraction_contaminant",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &ThermalTwoPhaseFlowComponentialLocalAssemblerInterface::
                getIntPtGasMolFracContaminant));
    
}

void ThermalTwoPhaseFlowComponentialProcess::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble ThermalTwoPhaseFlowComponentialProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_table, t, dt, x, xdot, process_id, M, K,
        b, _coupled_solutions);
}

void ThermalTwoPhaseFlowComponentialProcess::assembleWithJacobianConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    GlobalVector const& xdot, const double dxdot_dx, const double dx_dx,
    int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
    GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian ThermalTwoPhaseFlowComponentialProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
       dof_table = {std::ref(*_local_to_global_index_map)};
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, pv.getActiveElementIDs(), dof_table, t, dt, x, xdot,
        dxdot_dx, dx_dx, process_id, M, K, b, Jac, _coupled_solutions);
}

void ThermalTwoPhaseFlowComponentialProcess::preTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x, double const t, double const dt,
    const int process_id)
{
    DBUG("PreTimestep ThermalTwoPhaseFlowComponentialProcess.");

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];

    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        pv.getActiveElementIDs(), *_local_to_global_index_map, *x[process_id],
        t, dt);
}

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
