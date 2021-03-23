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

#include "ProcessLib/Process.h"
#include "ThermalTwoPhaseFlowComponentialLocalAssembler.h"

namespace MathLib
{
class PiecewiseLinearInterpolation;
}
namespace MeshLib
{
class Mesh;
}  // namespace MeshLib
namespace NumLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowComponential
{
/**
 * \brief A class to simulate the non-isothermal componential two-phase flow process in
 * porous media
 */
class ThermalTwoPhaseFlowComponentialProcess final : public Process
{
public:
    ThermalTwoPhaseFlowComponentialProcess(
        std::string name,
        MeshLib::Mesh& mesh,
        std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        ThermalTwoPhaseFlowComponentialProcessData&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves);

    bool isLinear() const override { return false; }

    Eigen::Vector3d getFluxTH2(
        std::size_t const element_id, int pv_index,
                            MathLib::Point3d const& p, double const t,
                            std::vector<GlobalVector*> const& x) const override;

private:
    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, double const dt,
                                 std::vector<GlobalVector*> const& x,
                                 std::vector<GlobalVector*> const& xdot,
                                 int const process_id, GlobalMatrix& M,
                                 GlobalMatrix& K, GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, double const dt, std::vector<GlobalVector*> const& x,
        std::vector<GlobalVector*> const& xdot, const double dxdot_dx,
        const double dx_dx,
        int const process_id, GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b,
        GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(std::vector<GlobalVector*> const& x,
                                    const double t, const double dt,
                                    const int process_id) override;

    ThermalTwoPhaseFlowComponentialProcessData _process_data;

    std::vector<std::unique_ptr<ThermalTwoPhaseFlowComponentialLocalAssemblerInterface>>
        _local_assemblers;
};

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
