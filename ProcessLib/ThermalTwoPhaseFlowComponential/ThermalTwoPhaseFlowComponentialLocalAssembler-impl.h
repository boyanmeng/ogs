/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// TODO: nomenclature

#pragma once

#include "ThermalTwoPhaseFlowComponentialLocalAssembler.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/Function/Interpolation.h"
#include "ThermalTwoPhaseFlowComponentialProcessData.h"

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowComponential
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
void ThermalTwoPhaseFlowComponentialLocalAssembler<
    ShapeFunction, IntegrationMethod,
    GlobalDim>::assemble(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& /*local_xdot*/,
                         std::vector<double>& local_M_data,
                         std::vector<double>& local_K_data,
                         std::vector<double>& local_b_data)
{
    using MaterialLib::PhysicalConstant::IdealGasConstant;
    auto const& water_mol_mass =
        MaterialLib::PhysicalConstant::MolarMass::Water;
    auto const& air_mol_mass = MaterialLib::PhysicalConstant::MolarMass::Air;

    auto const local_matrix_size = local_x.size();

    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    auto Mwp =
        local_M.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Mwa =
        local_M.template block<nonwet_pressure_size, overall_mol_frac_air_size>(
            nonwet_pressure_matrix_index, overall_mol_frac_air_matrix_index);
    auto Mwt = local_M.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Map =
        local_M.template block<overall_mol_frac_air_size, nonwet_pressure_size>(
            overall_mol_frac_air_matrix_index, nonwet_pressure_matrix_index);
    auto Maa = local_M.template block<overall_mol_frac_air_size,
                                       overall_mol_frac_air_size>(
        overall_mol_frac_air_matrix_index, overall_mol_frac_air_matrix_index);
    auto Mat =
        local_M.template block<overall_mol_frac_air_size, temperature_size>(
            overall_mol_frac_air_matrix_index, temperature_matrix_index);

    auto Mep = local_M.template block<temperature_size,
                                      nonwet_pressure_size>(
        temperature_matrix_index,
        nonwet_pressure_matrix_index);
    auto Mea =
        local_M.template block<temperature_size,
                                       overall_mol_frac_air_size>(
            temperature_matrix_index,
        overall_mol_frac_air_matrix_index);
    auto Met = local_M.template block<temperature_size,
                                      temperature_size>(
        temperature_matrix_index, temperature_matrix_index);
    
    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kwp =
        local_K.template block<nonwet_pressure_size, nonwet_pressure_size>(
            nonwet_pressure_matrix_index, nonwet_pressure_matrix_index);
    auto Kwa =
        local_K.template block<nonwet_pressure_size, overall_mol_frac_air_size>(
            nonwet_pressure_matrix_index, overall_mol_frac_air_matrix_index);
    auto Kwt = local_K.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Kap =
        local_K.template block<overall_mol_frac_air_size, nonwet_pressure_size>(
            overall_mol_frac_air_matrix_index, nonwet_pressure_matrix_index);
    auto Kaa = local_K.template block<overall_mol_frac_air_size,
                                       overall_mol_frac_air_size>(
        overall_mol_frac_air_matrix_index, overall_mol_frac_air_matrix_index);
    auto Kat =
        local_K.template block<overall_mol_frac_air_size, temperature_size>(
            overall_mol_frac_air_matrix_index, temperature_matrix_index);

    auto Kep = local_K.template block<temperature_size, nonwet_pressure_size>(
        temperature_matrix_index, nonwet_pressure_matrix_index);
    auto Kea =
        local_K.template block<temperature_size, overall_mol_frac_air_size>(
            temperature_matrix_index, overall_mol_frac_air_matrix_index);
    auto Ket = local_K.template block<temperature_size, temperature_size>(
        temperature_matrix_index, temperature_matrix_index);

    auto Bw = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);

    auto Ba = local_b.template segment<overall_mol_frac_air_size>(
        overall_mol_frac_air_matrix_index);

    auto Be = local_b.template segment<temperature_size>(
        temperature_matrix_index);

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& gas_phase = medium.phase("Gas");
    auto const& solid_phase = medium.phase("Solid");

    // components
    auto const& liquid_water = liquid_phase.component("wasser");
    auto const& dissolved_air = liquid_phase.component("air");    
    auto const& water_vapor = gas_phase.component("wasser");
    auto const& gaseous_air = gas_phase.component("air");

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const num_nodes = ShapeFunction::NPOINTS;
    // only needed for advective form
    /* auto const pg_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes);
    auto const Xa_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);
    auto const Xc_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[2 * num_nodes], num_nodes);
    auto const T_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[3 * num_nodes], num_nodes);
    */

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip); // TODO: necessary?
        auto const& sm = _shape_matrices[ip];

        double pg_int_pt = 0.;
        double Xa_int_pt =
            0.;  // total molar fraction of the light component air
        double T_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pg_int_pt,
                                         Xa_int_pt,
                                         T_int_pt);

        double& Sw = _ip_data[ip].sw;
        double& x_water_wet = _ip_data[ip].x_w_L;
        double& x_air_wet = _ip_data[ip].x_a_L;
        double& dsw_dpg = _ip_data[ip].dsw_dpg;
        double& dxwG_dpg = _ip_data[ip].dxwG_dpg;
        double& dxaG_dpg = _ip_data[ip].dxaG_dpg;
        double& dsw_dXa = _ip_data[ip].dsw_dXa;
        double& dxwG_dXa = _ip_data[ip].dxwG_dXa;
        double& dxaG_dXa = _ip_data[ip].dxaG_dXa;
        double& dsw_dT = _ip_data[ip].dsw_dT;
        double& dxwG_dT = _ip_data[ip].dxwG_dT;
        double& dxaG_dT = _ip_data[ip].dxaG_dT;
        double& dxwL_dpg = _ip_data[ip].dxwL_dpg;
        double& dxaL_dpg = _ip_data[ip].dxaL_dpg;
        double& dxwL_dXa = _ip_data[ip].dxwL_dXa;
        double& dxaL_dXa = _ip_data[ip].dxaL_dXa;
        double& dxwL_dT = _ip_data[ip].dxwL_dT;
        double& dxaL_dT = _ip_data[ip].dxaL_dT;

        // get reference
        auto& capillary_pressure_model =
            medium.property(MPL::PropertyType::capillary_pressure);

        MPL::VariableArray variables;
        variables[static_cast<int>(MPL::Variable::phase_pressure)] = pg_int_pt;
        variables[static_cast<int>(MPL::Variable::temperature)] = T_int_pt;

        // if water density depends on T, NCP needs to be modified
        double const density_water =
            liquid_water.property(MPL::PropertyType::density)
                .template value<double>(variables, pos, t, dt);

        // TODO: should be able to read H_ref and delta from input file.
        double const henry_air =
            dissolved_air.property(MPL::PropertyType::henry_constant)
                .template value<double>(variables, pos, t, dt);

         const double error_tolerance = _process_data.error_tolerance;
        // Calculate Sw, x_w_L, x_a_L, x_c_L and various derivatives from PVs
        if (!_process_data.material->computeConstitutiveRelation(        
                t,
                dt,
                pos,
                capillary_pressure_model,
                error_tolerance,
                density_water,
                henry_air,
                pg_int_pt,
                Xa_int_pt,
                T_int_pt,
                Sw,
                x_water_wet,
                x_air_wet,
                dsw_dpg, 
                dxwG_dpg,
                dxaG_dpg,
                dsw_dXa, 
                dxwG_dXa,
                dxaG_dXa,
                dsw_dT,
                dxwG_dT, 
                dxaG_dT,
                dxwL_dpg,
                dxaL_dpg,
                dxwL_dXa,
                dxaL_dXa,
                dxwL_dT, 
                dxaL_dT))
        {
            Sw = 0.0;
            x_water_wet = 0.99;
            x_air_wet = 0.01;
            if (!_process_data.material->computeConstitutiveRelation(
                    t,
                    dt,
                    pos,
                    capillary_pressure_model,
                    error_tolerance,
                    density_water,
                    henry_air,
                    pg_int_pt,
                    Xa_int_pt,
                    T_int_pt,
                    Sw,
                    x_water_wet,
                    x_air_wet,
                    dsw_dpg,
                    dxwG_dpg,
                    dxaG_dpg,
                    dsw_dXa,
                    dxwG_dXa,
                    dxaG_dXa,
                    dsw_dT,
                    dxwG_dT,
                    dxaG_dT,
                    dxwL_dpg,
                    dxaL_dpg,
                    dxwL_dXa,
                    dxaL_dXa,
                    dxwL_dT,
                    dxaL_dT))
            {
                OGS_FATAL("Computation of local constitutive relation failed.");
            }
        }

        variables[static_cast<int>(MPL::Variable::liquid_saturation)] =
            std::clamp(Sw, 0., 1.);

         auto pc = capillary_pressure_model.template value<double>(variables,
                                                                  pos, t, dt);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] = pc;

        _capillary_pressure[ip] = pc;

        _saturation[ip] = Sw;
        _pressure_wetting[ip] = pg_int_pt - pc;
        _liquid_molar_fraction_air[ip] = x_air_wet;
        // TODO: remaining SVs

        double const mol_density_water = density_water / water_mol_mass;
        double const mol_density_wet = mol_density_water;
        double const density_air_wet = mol_density_wet * air_mol_mass * x_air_wet;
        double const density_wet = density_air_wet + density_water;

        double const p_vapor_nonwet =
            _process_data.material->calculateVaporPressureNonwet(pc, T_int_pt,
                                                                 density_water);
        double const x_water_nonwet = p_vapor_nonwet / pg_int_pt * x_water_wet;
        _gas_molar_fraction_water[ip] = x_water_nonwet;

        double const x_air_nonwet =
            mol_density_wet / henry_air / pg_int_pt * x_air_wet;
        double const mol_density_nonwet =
            pg_int_pt / IdealGasConstant / T_int_pt;
        double const d_mol_density_nonwet_dpg =
            1. / IdealGasConstant / T_int_pt;
        double const d_mol_density_nonwet_dT =
            -pg_int_pt / IdealGasConstant / T_int_pt / T_int_pt;
        double const density_water_nonwet =
            mol_density_nonwet * water_mol_mass * x_water_nonwet;
        double const density_air_nonwet =
            mol_density_nonwet * air_mol_mass * x_air_nonwet;
        double const density_nonwet =
            density_air_nonwet + density_water_nonwet;

        double const mol_density_tot =
            Sw * mol_density_wet + (1 - Sw) * mol_density_nonwet;
        double const d_mol_density_tot_dpg =
            (mol_density_wet - mol_density_nonwet) * dsw_dpg +
            (1 - Sw) * d_mol_density_nonwet_dpg;
        double const d_mol_density_tot_dXa =
            (mol_density_wet - mol_density_nonwet) * dsw_dXa;
        double const d_mol_density_tot_dT =
            (mol_density_wet - mol_density_nonwet) * dsw_dT +
            (1 - Sw) * d_mol_density_nonwet_dT;

        auto dPC_dSw = capillary_pressure_model
                .template dValue<double>(
                    variables, MPL::Variable::liquid_saturation, pos, t, dt);

        auto const porosity =
            medium.property(MPL::PropertyType::porosity)
                .template value<double>(variables, pos, t, dt);

        // Assemble M matrix
        Mwp.noalias() +=
            porosity *
            (mol_density_wet * (x_water_wet * dsw_dpg + Sw * dxwL_dpg) +
             (1 - Sw) * x_water_nonwet * d_mol_density_nonwet_dpg +
             mol_density_nonwet *
                 ((1 - Sw) * dxwG_dpg - x_water_nonwet * dsw_dpg)) *
            _ip_data[ip].mass_operator;
        Mwa.noalias() +=
            porosity *
            (mol_density_wet * (x_water_wet * dsw_dXa + Sw * dxwL_dXa) +
             mol_density_nonwet *
                 ((1 - Sw) * dxwG_dXa - x_water_nonwet * dsw_dXa)) *
            _ip_data[ip].mass_operator;
        Mwt.noalias() +=
            porosity *
            (mol_density_wet * (x_water_wet * dsw_dT + Sw * dxwL_dT) +
             d_mol_density_nonwet_dT * (1 - Sw) * x_water_nonwet +
             mol_density_nonwet *
                 ((1 - Sw) * dxwG_dT - x_water_nonwet * dsw_dT)) *
            _ip_data[ip].mass_operator;
        Map.noalias() += porosity * Xa_int_pt * d_mol_density_tot_dpg *
                         _ip_data[ip].mass_operator;
        Maa.noalias() += porosity *
                         (Xa_int_pt * d_mol_density_tot_dXa + mol_density_tot) *
                         _ip_data[ip].mass_operator;     
        Mat.noalias() += porosity * Xa_int_pt * d_mol_density_tot_dT *
                         _ip_data[ip].mass_operator;

        double k_rel_wet = std::clamp(std::pow(Sw, 3.), 0., 1.);
        double k_rel_nonwet = std::clamp(std::pow(1 - Sw, 3.), 0., 1.);
        
        auto const mu_nonwet =
            gas_phase.property(MPL::PropertyType::viscosity)
                .template value<double>(variables, pos, t, dt);
        double const lambda_nonwet = k_rel_nonwet / mu_nonwet;
        auto const mu_wet = liquid_phase.property(MPL::PropertyType::viscosity)
                                .template value<double>(variables, pos, t, dt);
        double const lambda_wet = k_rel_wet / mu_wet;

        auto const K_int = MPL::formEigenTensor<GlobalDim>(
            medium.property(MPL::PropertyType::permeability)
                .template value<double>(variables, pos, t, dt));

        // only needed for advective form
        /* 
        GlobalDimVectorType const velocity_nonwet =
            -lambda_nonwet * K_int *
            (sm.dNdx * pg_nodal_values);
        GlobalDimVectorType const velocity_wet =
            -lambda_wet * K_int *
            ((1 - dPC_dSw * dsw_dpg) * sm.dNdx * pg_nodal_values -
             dPC_dSw * dsw_dXa * sm.dNdx * Xa_nodal_values -
             dPC_dSw * dsw_dXc * sm.dNdx * Xc_nodal_values -
             dPC_dSw * dsw_dT * sm.dNdx * T_nodal_values);
        */

        // diffusion coefficients
        double const diffusion_coeff_wet =
            liquid_phase.property(MPL::PropertyType::molecular_diffusion)
                .template value<double>(variables, pos, t, dt);
        double const diffusion_coeff_nonwet =
            gas_phase.property(MPL::PropertyType::molecular_diffusion)
                .template value<double>(variables, pos, t, dt);

        // heat capacity
        double const heat_capacity_water =
            water_vapor.property(MPL::PropertyType::heat_capacity)
                .template value<double>(variables, pos, t, dt);
        double const heat_capacity_air =
            gaseous_air.property(MPL::PropertyType::heat_capacity)
                .template value<double>(variables, pos, t, dt);
        double const heat_capacity_solid =
            solid_phase.property(MPL::PropertyType::heat_capacity)
                .template value<double>(variables, pos, t, dt);

        auto const effective_thermal_conductivity =
            medium.property(MPL::PropertyType::thermal_conductivity)
                .template value<double>(variables, pos, t, dt);

        double const density_solid =
            solid_phase.property(MPL::PropertyType::density)
                .template value<double>(variables, pos, t, dt);

        double const latent_heat_evaporation = 2258000.;

        // enthalpy
        double const enthalpy_water_nonwet =
            _process_data.material->getWaterVaporEnthalpySimple(
                T_int_pt, heat_capacity_water, pg_int_pt,
                latent_heat_evaporation);
        double const enthalpy_air_nonwet =
            _process_data.material->getAirEnthalpySimple(
                T_int_pt, heat_capacity_air, pg_int_pt);
        double const mol_mass_nonwet =
            x_water_nonwet * water_mol_mass + x_air_nonwet * air_mol_mass;
        double const X_water_nonwet =
            x_water_nonwet * water_mol_mass / mol_mass_nonwet;
        double const X_air_nonwet = 1 - X_water_nonwet;
        double const enthalpy_nonwet =
            X_water_nonwet * enthalpy_water_nonwet +
            X_air_nonwet * enthalpy_air_nonwet;
        double const internal_energy_nonwet =
            enthalpy_nonwet - pg_int_pt / density_nonwet;

        double const enthalpy_wet =
            _process_data.material->getLiquidWaterEnthalpySimple(
                T_int_pt, heat_capacity_water, _pressure_wetting[ip]);
        double const internal_energy_wet = enthalpy_wet;
        double const volumetric_enthalpy_nonwet =
            density_nonwet * enthalpy_nonwet;
        double const d_vol_enthalpy_nonwet_dpg =
            (x_water_nonwet * water_mol_mass * enthalpy_water_nonwet +
             x_air_nonwet * air_mol_mass * enthalpy_air_nonwet +
             pg_int_pt * (dxwG_dpg * water_mol_mass * enthalpy_water_nonwet +
                          dxaG_dpg * air_mol_mass * enthalpy_air_nonwet)) /
            IdealGasConstant / T_int_pt;
        double const d_vol_enthalpy_nonwet_dXa =
            mol_density_nonwet *
            (dxwG_dXa * water_mol_mass * enthalpy_water_nonwet +
             dxaG_dXa * air_mol_mass * enthalpy_air_nonwet);
        double const d_vol_enthalpy_nonwet_dT =
            mol_density_nonwet *
            (-(x_water_nonwet * water_mol_mass * enthalpy_water_nonwet +
               x_air_nonwet * air_mol_mass * enthalpy_air_nonwet) /
                 T_int_pt +
             water_mol_mass * (heat_capacity_water * x_water_nonwet +
                               enthalpy_water_nonwet * dxwG_dT) +
             (air_mol_mass * heat_capacity_air + IdealGasConstant) *
                 x_air_nonwet +
             air_mol_mass * enthalpy_air_nonwet * dxaG_dT);

        laplace_operator.noalias() = sm.dNdx.transpose() * K_int *
                                     sm.dNdx * _ip_data[ip].integration_weight;

        // Assemble energy equation M matrix

        Mep.noalias() +=
            porosity *
            ((density_wet * internal_energy_wet -
              density_nonwet * enthalpy_nonwet + pg_int_pt) *
                 dsw_dpg +
             (1 - Sw) * (d_vol_enthalpy_nonwet_dpg - 1)) *
            _ip_data[ip].mass_operator;
        Mea.noalias() +=
            porosity * (density_wet * internal_energy_wet * dsw_dXa +
                          (1 - Sw) * d_vol_enthalpy_nonwet_dXa) *
                         _ip_data[ip].mass_operator;
        Met.noalias() += porosity * (density_wet * (internal_energy_wet * dsw_dT +
                                       Sw * heat_capacity_water) +
                              (1 - Sw) * d_vol_enthalpy_nonwet_dT) *
                _ip_data[ip].mass_operator +
            (1 - porosity) *
                density_solid * heat_capacity_solid *
                _ip_data[ip].mass_operator;

        // Assemble K matrix
        Kwp.noalias() += (mol_density_wet * x_water_wet * lambda_wet *
                              (1 - dPC_dSw * dsw_dpg) +
                          mol_density_nonwet * x_water_nonwet * lambda_nonwet) *
                             laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxwL_dpg +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxwG_dpg) *
                _ip_data[ip].diffusion_operator;
        Kwa.noalias() +=
            (-mol_density_wet * x_water_wet * lambda_wet * dPC_dSw *
             dsw_dXa)*
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxwL_dXa +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxwG_dXa) *
                _ip_data[ip].diffusion_operator;
        Kwt.noalias() +=
            (-mol_density_wet * x_water_wet * lambda_wet * dPC_dSw * dsw_dT) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxwL_dT +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxwG_dT) *
                _ip_data[ip].diffusion_operator;
        Kap.noalias() +=
            (mol_density_wet * x_air_wet * lambda_wet *
                 (1 - dPC_dSw * dsw_dpg) +
             mol_density_nonwet * x_air_nonwet * lambda_nonwet) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxaL_dpg +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxaG_dpg) *
                _ip_data[ip].diffusion_operator;
        Kaa.noalias() +=
            (-mol_density_wet * x_air_wet * lambda_wet * dPC_dSw * dsw_dXa) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxaL_dXa +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxaG_dXa) *
                _ip_data[ip].diffusion_operator;
        Kat.noalias() +=
            (-mol_density_wet * x_air_wet * lambda_wet * dPC_dSw * dsw_dT) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxaL_dT +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxaG_dT) *
                _ip_data[ip].diffusion_operator;
        Kep.noalias() += (lambda_nonwet * density_nonwet * enthalpy_nonwet +
                          lambda_wet * density_wet * enthalpy_wet * (1 - dPC_dSw * dsw_dpg)) *
                         laplace_operator;
        Kea.noalias() +=
            (-lambda_wet * density_wet * enthalpy_wet * dPC_dSw * dsw_dXa) *
            laplace_operator;
        Ket.noalias() +=
            (-lambda_wet * density_wet * enthalpy_wet * dPC_dSw * dsw_dT) *
                laplace_operator +
            sm.dNdx.transpose() * effective_thermal_conductivity *
                sm.dNdx * _ip_data[ip].integration_weight;    
        
        if (_process_data.has_gravity)
        {
            auto const& b = _process_data.specific_body_force;
            NodalVectorType gravity_operator = sm.dNdx.transpose() *
                                               K_int * b *
                                               _ip_data[ip].integration_weight;
            Bw.noalias() +=
                (mol_density_wet * x_water_wet * lambda_wet * density_wet +
                 mol_density_nonwet * x_water_nonwet * lambda_nonwet *
                     density_nonwet) * gravity_operator;
            Ba.noalias() +=
                (mol_density_wet * x_air_wet * lambda_wet * density_wet +
                 mol_density_nonwet * x_air_nonwet * lambda_nonwet *
                     density_nonwet) *
                gravity_operator;
            Be.noalias() +=
                (lambda_nonwet * density_nonwet * density_nonwet *
                     enthalpy_nonwet +
                 lambda_wet * density_wet * density_wet * enthalpy_wet) *
                gravity_operator;
        }  // end of has gravity   
    }
    if (_process_data.has_mass_lumping)
    {
        for (unsigned row = 0; row < Mwp.cols(); row++)
        {
            for (unsigned column = 0; column < Mwp.cols(); column++)
            {
                if (row != column)
                {
                    Mwp(row, row) += Mwp(row, column);
                    Mwp(row, column) = 0.0;
                    Mwa(row, row) += Mwa(row, column);
                    Mwa(row, column) = 0.0;
                    Mwt(row, row) += Mwt(row, column);
                    Mwt(row, column) = 0.0;
                    Map(row, row) += Map(row, column);
                    Map(row, column) = 0.0;
                    Maa(row, row) += Maa(row, column);
                    Maa(row, column) = 0.0;
                    Mat(row, row) += Mat(row, column);
                    Mat(row, column) = 0.0;
                    Mep(row, row) += Mep(row, column);
                    Mep(row, column) = 0.0;
                    Mea(row, row) += Mea(row, column);
                    Mea(row, column) = 0.0;
                    Met(row, row) += Met(row, column);
                    Met(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping

}

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
