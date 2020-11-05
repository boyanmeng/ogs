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

#include <vector>
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "ThermalTwoPhaseFlowComponentialProcessData.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowComponential
{
template <typename NodalMatrixType>
struct IntegrationPointData final
{
    explicit IntegrationPointData()
        : sw(0.9),
          x_w_L(0.99),
          x_a_L(0.01),
          dsw_dpl(0.0),
          dxwG_dpl(0.0),
          dxaG_dpl(0.0),
          dsw_dXa(0.0),
          dxwG_dXa(0.0),
          dxaG_dXa(0.0),
          dsw_dT(0.0),
          dxwG_dT(0.0),
          dxaG_dT(0.0),
          dxwL_dpl(0.0),
          dxaL_dpl(0.0),
          dxwL_dXa(0.0),
          dxaL_dXa(0.0),
          dxwL_dT(0.0),
          dxaL_dT(0.0)
    {
    }
    // ThermalTwoPhaseFlowComponentialMaterialProperties& mat_property;
    double sw;
    double x_w_L;
    double x_a_L;
    double dsw_dpl;
    double dxwG_dpl;
    double dxaG_dpl;
    double dsw_dXa;
    double dxwG_dXa;
    double dxaG_dXa;
    double dsw_dT;
    double dxwG_dT;
    double dxaG_dT;
    double dxwL_dpl;
    double dxaL_dpl;
    double dxwL_dXa;
    double dxaL_dXa;
    double dxwL_dT;
    double dxaL_dT;

    double integration_weight;
    NodalMatrixType mass_operator;
    NodalMatrixType diffusion_operator;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
const unsigned NUM_NODAL_DOF = 3;

class ThermalTwoPhaseFlowComponentialLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtSaturation(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtNonwettingPressure(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> const& getIntPtCapillaryPressure(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> const& getIntPtLiquidMolFracAir(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> const& getIntPtGasMolFracWater(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class ThermalTwoPhaseFlowComponentialLocalAssembler
    : public ThermalTwoPhaseFlowComponentialLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using LocalMatrixType = typename LocalAssemblerTraits::LocalMatrix;
    using LocalVectorType = typename LocalAssemblerTraits::LocalVector;

public:
    ThermalTwoPhaseFlowComponentialLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermalTwoPhaseFlowComponentialProcessData const& process_data)
        : _element(element),
          _integration_method(integration_order),
          _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                            IntegrationMethod, GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _process_data(process_data),
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_nonwetting(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _capillary_pressure(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _liquid_molar_fraction_air(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _gas_molar_fraction_water(
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);  // reserve memory space, size is still 0
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back();   // create object at end, ip_data.size() = 1->2->3->4
            auto const& sm = _shape_matrices[ip];  // size = n_integration_points
            _ip_data[ip].integration_weight =
                sm.integralMeasure * sm.detJ *
                _integration_method.getWeightedPoint(ip).getWeight();
            _ip_data[ip].mass_operator.setZero(ShapeFunction::NPOINTS,
                                              ShapeFunction::NPOINTS);
            _ip_data[ip].diffusion_operator.setZero(ShapeFunction::NPOINTS,
                                                   ShapeFunction::NPOINTS);
            _ip_data[ip].mass_operator.noalias() =
                sm.N.transpose() * sm.N * _ip_data[ip].integration_weight;
            _ip_data[ip].diffusion_operator.noalias() =
                sm.dNdx.transpose() * sm.dNdx * _ip_data[ip].integration_weight;
        }
    }

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& local_xdot,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _shape_matrices[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSaturation(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_saturation.empty());
        return _saturation;
    }

    std::vector<double> const& getIntPtNonwettingPressure(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_pressure_nonwetting.empty());
        return _pressure_nonwetting;
    }

    std::vector<double> const& getIntPtCapillaryPressure(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_capillary_pressure.empty());
        return _capillary_pressure;
    }

    std::vector<double> const& getIntPtLiquidMolFracAir(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_liquid_molar_fraction_air.empty());
        return _liquid_molar_fraction_air;
    }

    std::vector<double> const& getIntPtGasMolFracWater(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_gas_molar_fraction_water.empty());
        return _gas_molar_fraction_water;
    }

private:
    MeshLib::Element const& _element;

    IntegrationMethod const _integration_method;
    std::vector<ShapeMatrices, Eigen::aligned_allocator<ShapeMatrices>>
        _shape_matrices;

    ThermalTwoPhaseFlowComponentialProcessData const& _process_data;
    std::vector<IntegrationPointData<NodalMatrixType>,
                Eigen::aligned_allocator<IntegrationPointData<NodalMatrixType>>>
        _ip_data;

    std::vector<double> _saturation;  /// used for secondary variable output
    std::vector<double> _pressure_nonwetting;
    std::vector<double> _capillary_pressure;
    std::vector<double> _liquid_molar_fraction_air;
    std::vector<double> _gas_molar_fraction_water;

    static const int liquid_pressure_matrix_index = 0;
    static const int overall_mol_frac_air_matrix_index = ShapeFunction::NPOINTS;
    static const int temperature_matrix_index = 2 * ShapeFunction::NPOINTS;

    static const int liquid_pressure_size = ShapeFunction::NPOINTS;
    static const int overall_mol_frac_air_size = ShapeFunction::NPOINTS;
    static const int temperature_size = ShapeFunction::NPOINTS;
};

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib

#include "ThermalTwoPhaseFlowComponentialLocalAssembler-impl.h"
