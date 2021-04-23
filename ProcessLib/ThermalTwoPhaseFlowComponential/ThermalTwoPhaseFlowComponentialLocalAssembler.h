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
#include "NumLib/Fem/InitShapeMatrices.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"

#include "ThermalTwoPhaseFlowComponentialProcessData.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowComponential
{
template <typename NodalMatrixType>
struct IntegrationPointData final
{
    double integration_weight;
    NodalMatrixType mass_operator;
    NodalMatrixType diffusion_operator;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};
const unsigned NUM_NODAL_DOF = 4;

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

    virtual std::vector<double> const& getIntPtWettingPressure(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> const& getIntPtLiquidMolFracContaminant(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> const& getIntPtGasMolFracWater(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> const& getIntPtGasMolFracContaminant(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> const& getIntPtRelPermWet(
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
          _shape_matrices(
              NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                        GlobalDim>(
              element, is_axially_symmetric, _integration_method)),
          _process_data(process_data),
          _saturation(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _pressure_wetting(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _liquid_molar_fraction_contaminant(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _gas_molar_fraction_water(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _gas_molar_fraction_contaminant(
              std::vector<double>(_integration_method.getNumberOfPoints())),
          _krel_wet(
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

    std::vector<double> const& getIntPtWettingPressure(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_pressure_wetting.empty());
        return _pressure_wetting;
    }

    std::vector<double> const& getIntPtLiquidMolFracContaminant(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_liquid_molar_fraction_contaminant.empty());
        return _liquid_molar_fraction_contaminant;
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

    std::vector<double> const& getIntPtGasMolFracContaminant(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_gas_molar_fraction_contaminant.empty());
        return _gas_molar_fraction_contaminant;
    }

    std::vector<double> const& getIntPtRelPermWet(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override
    {
        assert(!_krel_wet.empty());
        return _krel_wet;
    }

    Eigen::Vector3d getFluxTH2(
        MathLib::Point3d const& pnt_local_coords, int const pv_index,
        double const t, std::vector<double> const& local_x) const override
    {
        using MaterialLib::PhysicalConstant::IdealGasConstant;
        auto const& water_mol_mass =
            MaterialLib::PhysicalConstant::MolarMass::Water;
        auto const& air_mol_mass =
            MaterialLib::PhysicalConstant::MolarMass::Air;

        auto const local_pc = Eigen::Map<const NodalVectorType>(
            &local_x[capillary_pressure_matrix_index], capillary_pressure_size);
        auto const local_pg = Eigen::Map<const NodalVectorType>(
            &local_x[gas_pressure_matrix_index], gas_pressure_size);
        auto const local_Xc = Eigen::Map<const NodalVectorType>(
            &local_x[overall_mol_frac_contaminant_matrix_index],
            overall_mol_frac_contaminant_size);
        auto const local_T = Eigen::Map<const NodalVectorType>(
            &local_x[temperature_matrix_index], temperature_size);

        // Eval shape matrices at given point
        // Note: Axial symmetry is set to false here, because we only need dNdx
        // here, which is not affected by axial symmetry.
        auto const shape_matrices =
            NumLib::computeShapeMatrices<ShapeFunction, ShapeMatricesType,
                                         GlobalDim>(
                _element, false /*is_axially_symmetric*/,
                std::array{pnt_local_coords})[0];

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        MPL::VariableArray vars;

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");
        auto const& gas_phase = medium.phase("Gas");
        auto const& solid_phase = medium.phase("Solid");

        // components
        auto const& liquid_water = liquid_phase.component("wasser");
        auto const& dissolved_contaminant =
            liquid_phase.component("contaminant");
        auto const& water_vapor = gas_phase.component("wasser");
        auto const& gaseous_air = gas_phase.component("air");
        auto const& contaminant_vapor = gas_phase.component("contaminant");

        double pc_int_pt;
        NumLib::shapeFunctionInterpolate(local_pc, shape_matrices.N, pc_int_pt);
        double pg_int_pt;
        NumLib::shapeFunctionInterpolate(local_pg, shape_matrices.N, pg_int_pt);
        double Xc_int_pt;
        NumLib::shapeFunctionInterpolate(local_Xc, shape_matrices.N, Xc_int_pt);
        double T_int_pt;
        NumLib::shapeFunctionInterpolate(local_T, shape_matrices.N, T_int_pt);
        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T_int_pt;
        vars[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pc_int_pt;

        double const dt = std::numeric_limits<double>::quiet_NaN();

        double const density_water =
            liquid_water.property(MPL::PropertyType::density)
                .template value<double>(vars, pos, t, dt);

        double const pl = pg_int_pt - pc_int_pt;

        auto& saturation_model = medium.property(MPL::PropertyType::saturation);
        double const Sw =
            saturation_model.template value<double>(vars, pos, t, dt);
        vars[static_cast<int>(
            MaterialPropertyLib::Variable::liquid_saturation)] = Sw;

        double const ideal_gas_constant_times_T_int_pt =
            IdealGasConstant * T_int_pt;
        double const mol_density_water = density_water / water_mol_mass;
        double const mol_density_wet = mol_density_water;
        double const mol_density_nonwet =
            pg_int_pt / ideal_gas_constant_times_T_int_pt;
        double const mol_density_tot =
            Sw * mol_density_wet + (1 - Sw) * mol_density_nonwet;

        double const henry_contaminant =
            dissolved_contaminant.property(MPL::PropertyType::henry_constant)
                .template value<double>(vars, pos, t, dt);

        double const p_vapor_nonwet =
            _process_data.material->calculateVaporPressureNonwet(
                pc_int_pt, T_int_pt, density_water);

        double const k_c = pg_int_pt * henry_contaminant / mol_density_wet;
        double const k_w = pg_int_pt / p_vapor_nonwet;

        double const ntot_c =
            Sw * mol_density_wet * k_c + (1 - Sw) * mol_density_nonwet;
        double const x_contaminant_nonwet =
            Xc_int_pt * mol_density_tot / ntot_c;
        double const x_contaminant_wet = k_c * x_contaminant_nonwet;
        double const x_water_wet = 1 - x_contaminant_wet;
        double const x_water_nonwet = x_water_wet / k_w;
        double const x_air_nonwet = 1 - x_water_nonwet - x_contaminant_nonwet;

        double const contaminant_mol_mass =
            contaminant_vapor.property(MPL::PropertyType::molar_mass)
                .template value<double>(vars, pos, t, dt);
        double const density_contaminant_wet =
            mol_density_wet * contaminant_mol_mass * x_contaminant_wet;
        double const density_wet = density_water + density_contaminant_wet;

        double const density_water_nonwet =
            mol_density_nonwet * water_mol_mass * x_water_nonwet;
        double const density_air_nonwet =
            mol_density_nonwet * air_mol_mass * x_air_nonwet;
        double const density_contaminant_nonwet =
            mol_density_nonwet * contaminant_mol_mass * x_contaminant_nonwet;
        double const density_nonwet = density_air_nonwet +
                                      density_water_nonwet +
                                      density_contaminant_nonwet;

        auto const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<Eigen::Vector2d>(vars, pos, t, dt);
        auto const k_rel_wet = k_rel[0];
        auto const k_rel_nonwet = k_rel[1];

        auto const mu_nonwet = gas_phase.property(MPL::PropertyType::viscosity)
                                   .template value<double>(vars, pos, t, dt);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;
        auto const mu_wet = liquid_phase.property(MPL::PropertyType::viscosity)
                                .template value<double>(vars, pos, t, dt);
        double const lambda_wet = k_rel_wet / mu_wet;

        auto const K_int = MPL::formEigenTensor<GlobalDim>(
            medium.property(MPL::PropertyType::permeability)
                .value(vars, pos, t, dt));

        GlobalDimVectorType velocity_nonwet =
            -lambda_nonwet * K_int * shape_matrices.dNdx * local_pg;
        GlobalDimVectorType velocity_wet =
            -lambda_wet * K_int * shape_matrices.dNdx * (local_pg - local_pc);
        if (_process_data.has_gravity)
        {
            auto const& b = _process_data.specific_body_force;
            velocity_nonwet += lambda_nonwet * K_int * density_nonwet * b;
            velocity_wet += lambda_wet * K_int * density_wet * b;
        }

        double const mol_mass_nonwet =
            x_water_nonwet * water_mol_mass + x_air_nonwet * air_mol_mass +
            x_contaminant_nonwet * contaminant_mol_mass;
        double const X_water_nonwet =
            x_water_nonwet * water_mol_mass / mol_mass_nonwet;
        double const X_air_nonwet =
            x_air_nonwet * air_mol_mass / mol_mass_nonwet;
        double const X_contaminant_nonwet = 1 - X_water_nonwet - X_air_nonwet;

        // heat capacity
        double const heat_capacity_water =
            water_vapor.property(MPL::PropertyType::heat_capacity)
                .template value<double>(vars, pos, t, dt);
        double const heat_capacity_air =
            gaseous_air.property(MPL::PropertyType::heat_capacity)
                .template value<double>(vars, pos, t, dt);
        double const heat_capacity_contaminant =
            contaminant_vapor.property(MPL::PropertyType::heat_capacity)
                .template value<double>(vars, pos, t, dt);

        double const latent_heat_evaporation = 2258000.;

        // enthalpy
        double const enthalpy_water_nonwet =
            _process_data.material->getWaterVaporEnthalpySimple(
                T_int_pt, heat_capacity_water, pg_int_pt,
                latent_heat_evaporation);
        double const enthalpy_air_nonwet =
            _process_data.material->getAirEnthalpySimple(
                T_int_pt, heat_capacity_air, pg_int_pt);
        double const enthalpy_contaminant_nonwet =
            _process_data.material->getContaminantEnthalpySimple(
                T_int_pt, heat_capacity_contaminant, contaminant_mol_mass,
                pg_int_pt);

        double const enthalpy_nonwet =
            X_water_nonwet * enthalpy_water_nonwet +
            X_air_nonwet * enthalpy_air_nonwet +
            X_contaminant_nonwet * enthalpy_contaminant_nonwet;

        double const enthalpy_wet =
            _process_data.material->getLiquidWaterEnthalpySimple(
                T_int_pt, heat_capacity_water, pl);

        Eigen::Vector3d fluxTH2(0.0, 0.0, 0.0);
        switch (pv_index)
        {
            case 0:
                fluxTH2.head<GlobalDim>() =
                    mol_density_wet * x_water_wet * velocity_wet +
                    mol_density_nonwet * x_water_nonwet * velocity_nonwet;
                break;
            case 1:
                fluxTH2.head<GlobalDim>() =
                    mol_density_nonwet * x_air_nonwet * velocity_nonwet;
                break;
            case 2:
                fluxTH2.head<GlobalDim>() =
                    mol_density_wet * x_contaminant_wet * velocity_wet +
                    mol_density_nonwet * x_contaminant_nonwet * velocity_nonwet;
                break;
            case 3:
                fluxTH2.head<GlobalDim>() =
                    density_wet * enthalpy_wet * velocity_wet +
                    density_nonwet * enthalpy_nonwet * velocity_nonwet;
                break;
            default:
            {
            }
        }

        return fluxTH2;
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
    std::vector<double> _pressure_wetting;
    std::vector<double> _liquid_molar_fraction_contaminant;
    std::vector<double> _gas_molar_fraction_water;
    std::vector<double> _gas_molar_fraction_contaminant;
    std::vector<double> _krel_wet;

    static const int capillary_pressure_matrix_index = 0;
    static const int gas_pressure_matrix_index = ShapeFunction::NPOINTS;
    static const int overall_mol_frac_contaminant_matrix_index =
        2 * ShapeFunction::NPOINTS;
    static const int temperature_matrix_index = 3 * ShapeFunction::NPOINTS;

    static const int capillary_pressure_size = ShapeFunction::NPOINTS;
    static const int gas_pressure_size = ShapeFunction::NPOINTS;
    static const int overall_mol_frac_contaminant_size = ShapeFunction::NPOINTS;
    static const int temperature_size = ShapeFunction::NPOINTS;
};

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib

#include "ThermalTwoPhaseFlowComponentialLocalAssembler-impl.h"
