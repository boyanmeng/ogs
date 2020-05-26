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
    auto Mwc = local_M.template block<nonwet_pressure_size,
                                       overall_mol_frac_contaminant_size>(
        nonwet_pressure_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Mwt = local_M.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Map =
        local_M.template block<overall_mol_frac_air_size, nonwet_pressure_size>(
            overall_mol_frac_air_matrix_index, nonwet_pressure_matrix_index);
    auto Maa = local_M.template block<overall_mol_frac_air_size,
                                       overall_mol_frac_air_size>(
        overall_mol_frac_air_matrix_index, overall_mol_frac_air_matrix_index);
    auto Mac = local_M.template block<overall_mol_frac_air_size,
                                       overall_mol_frac_contaminant_size>(
        overall_mol_frac_air_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Mat =
        local_M.template block<overall_mol_frac_air_size, temperature_size>(
            overall_mol_frac_air_matrix_index, temperature_matrix_index);

    auto Mcp = local_M.template block<overall_mol_frac_contaminant_size,
                                      nonwet_pressure_size>(
        overall_mol_frac_contaminant_matrix_index,
        nonwet_pressure_matrix_index);
    auto Mca = local_M.template block<overall_mol_frac_contaminant_size,
                                       overall_mol_frac_air_size>(
        overall_mol_frac_contaminant_matrix_index,
        overall_mol_frac_air_matrix_index);
    auto Mcc = local_M.template block<overall_mol_frac_contaminant_size,
                                       overall_mol_frac_contaminant_size>(
        overall_mol_frac_contaminant_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Mct = local_M.template block<overall_mol_frac_contaminant_size,
                                      temperature_size>(
        overall_mol_frac_contaminant_matrix_index, temperature_matrix_index);

    auto Mep = local_M.template block<temperature_size,
                                      nonwet_pressure_size>(
        temperature_matrix_index,
        nonwet_pressure_matrix_index);
    auto Mea =
        local_M.template block<temperature_size,
                                       overall_mol_frac_air_size>(
            temperature_matrix_index,
        overall_mol_frac_air_matrix_index);
    auto Mec = local_M.template block<temperature_size,
                                       overall_mol_frac_contaminant_size>(
        temperature_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
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
    auto Kwc = local_K.template block<nonwet_pressure_size,
                                       overall_mol_frac_contaminant_size>(
        nonwet_pressure_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Kwt = local_K.template block<nonwet_pressure_size, temperature_size>(
        nonwet_pressure_matrix_index, temperature_matrix_index);

    auto Kap =
        local_K.template block<overall_mol_frac_air_size, nonwet_pressure_size>(
            overall_mol_frac_air_matrix_index, nonwet_pressure_matrix_index);
    auto Kaa = local_K.template block<overall_mol_frac_air_size,
                                       overall_mol_frac_air_size>(
        overall_mol_frac_air_matrix_index, overall_mol_frac_air_matrix_index);
    auto Kac = local_K.template block<overall_mol_frac_air_size,
                                       overall_mol_frac_contaminant_size>(
        overall_mol_frac_air_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Kat =
        local_K.template block<overall_mol_frac_air_size, temperature_size>(
            overall_mol_frac_air_matrix_index, temperature_matrix_index);

    auto Kcp = local_K.template block<overall_mol_frac_contaminant_size,
                                      nonwet_pressure_size>(
        overall_mol_frac_contaminant_matrix_index,
        nonwet_pressure_matrix_index);
    auto Kca = local_K.template block<overall_mol_frac_contaminant_size,
                                       overall_mol_frac_air_size>(
        overall_mol_frac_contaminant_matrix_index,
        overall_mol_frac_air_matrix_index);
    auto Kcc = local_K.template block<overall_mol_frac_contaminant_size,
                                       overall_mol_frac_contaminant_size>(
        overall_mol_frac_contaminant_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Kct = local_K.template block<overall_mol_frac_contaminant_size,
                                      temperature_size>(
        overall_mol_frac_contaminant_matrix_index, temperature_matrix_index);

    auto Kep = local_K.template block<temperature_size, nonwet_pressure_size>(
        temperature_matrix_index, nonwet_pressure_matrix_index);
    auto Kea =
        local_K.template block<temperature_size, overall_mol_frac_air_size>(
            temperature_matrix_index, overall_mol_frac_air_matrix_index);
    auto Kec = local_K.template block<temperature_size,
                                       overall_mol_frac_contaminant_size>(
        temperature_matrix_index, overall_mol_frac_contaminant_matrix_index);
    auto Ket = local_K.template block<temperature_size, temperature_size>(
        temperature_matrix_index, temperature_matrix_index);

    auto Bw = local_b.template segment<nonwet_pressure_size>(
        nonwet_pressure_matrix_index);

    auto Ba = local_b.template segment<overall_mol_frac_air_size>(
        overall_mol_frac_air_matrix_index);

    auto Bc = local_b.template segment<overall_mol_frac_contaminant_size>(
        overall_mol_frac_contaminant_matrix_index);

    auto Be = local_b.template segment<temperature_size>(
        temperature_matrix_index);

    auto const& medium = *_process_data.media_map->getMedium(_element.getID());
    auto const& liquid_phase = medium.phase("AqueousLiquid");
    auto const& gas_phase = medium.phase("Gas");
    auto const& solid_phase = medium.phase("Solid");

    // components
    auto const& liquid_water = liquid_phase.component("wasser");
    auto const& dissolved_air = liquid_phase.component("air");    
    auto const& dissolved_contaminant = liquid_phase.component("contaminant");
    auto const& water_vapor = gas_phase.component("wasser");
    auto const& gaseous_air = gas_phase.component("air");
    auto const& contaminant_vapor = gas_phase.component("contaminant");

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const num_nodes = ShapeFunction::NPOINTS;
    auto const pg_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip); // TODO: necessary?
        auto const& sm = _shape_matrices[ip];

        double pg_int_pt = 0.;
        double Xa_int_pt =
            0.;  // total molar fraction of the light component air
        double Xc_int_pt =
            0.;  // total molar fraction of the light component contaminant
        double T_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pg_int_pt,
                                         Xa_int_pt, Xc_int_pt,
                                         T_int_pt);

        double& Sw = _ip_data[ip].sw;
        double& x_water_wet = _ip_data[ip].x_w_L;
        double& x_air_wet = _ip_data[ip].x_a_L;
        double& x_contaminant_wet = _ip_data[ip].x_c_L;
        double& dsw_dpg = _ip_data[ip].dsw_dpg;
        double& dxwG_dpg = _ip_data[ip].dxwG_dpg;
        double& dxaG_dpg = _ip_data[ip].dxaG_dpg;
        double& dxcG_dpg = _ip_data[ip].dxcG_dpg;
        double& dsw_dXa = _ip_data[ip].dsw_dXa;
        double& dxwG_dXa = _ip_data[ip].dxwG_dXa;
        double& dxaG_dXa = _ip_data[ip].dxaG_dXa;
        double& dxcG_dXa = _ip_data[ip].dxcG_dXa;
        double& dsw_dXc = _ip_data[ip].dsw_dXc;
        double& dxwG_dXc = _ip_data[ip].dxwG_dXc;
        double& dxaG_dXc = _ip_data[ip].dxaG_dXc;
        double& dxcG_dXc = _ip_data[ip].dxcG_dXc;
        double& dsw_dT = _ip_data[ip].dsw_dT;
        double& dxwG_dT = _ip_data[ip].dxwG_dT;
        double& dxaG_dT = _ip_data[ip].dxaG_dT;
        double& dxcG_dT = _ip_data[ip].dxcG_dT;
        double& dxwL_dpg = _ip_data[ip].dxwL_dpg;
        double& dxaL_dpg = _ip_data[ip].dxaL_dpg;
        double& dxcL_dpg = _ip_data[ip].dxcL_dpg;
        double& dxwL_dXa = _ip_data[ip].dxwL_dXa;
        double& dxaL_dXa = _ip_data[ip].dxaL_dXa;
        double& dxcL_dXa = _ip_data[ip].dxcL_dXa;
        double& dxwL_dXc = _ip_data[ip].dxwL_dXc;
        double& dxaL_dXc = _ip_data[ip].dxaL_dXc;
        double& dxcL_dXc = _ip_data[ip].dxcL_dXc;
        double& dxwL_dT = _ip_data[ip].dxwL_dT;
        double& dxaL_dT = _ip_data[ip].dxaL_dT;
        double& dxcL_dT = _ip_data[ip].dxcL_dT;

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

        // Calculate Sw, x_w_L, x_a_L, x_c_L and various derivatives from PVs
        if (!_process_data.material->computeConstitutiveRelation(        
                t,
                dt,
                pos,
                capillary_pressure_model,
                density_water,
                pg_int_pt,
                Xa_int_pt,
                Xc_int_pt,
                T_int_pt,
                Sw,
                x_water_wet,
                x_air_wet,
                x_contaminant_wet,
                dsw_dpg, 
                dxwG_dpg,
                dxaG_dpg,
                dxcG_dpg,
                dsw_dXa, 
                dxwG_dXa,
                dxaG_dXa,
                dxcG_dXa,
                dsw_dXc, 
                dxwG_dXc,
                dxaG_dXc,
                dxcG_dXc,
                dsw_dT,
                dxwG_dT, 
                dxaG_dT,
                dxcG_dT,
                dxwL_dpg,
                dxaL_dpg,
                dxcL_dpg,
                dxwL_dXa,
                dxaL_dXa,
                dxcL_dXa,
                dxwL_dXc,
                dxaL_dXc,
                dxcL_dXc,
                dxwL_dT, 
                dxaL_dT,
                dxcL_dT))
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        variables[static_cast<int>(MPL::Variable::liquid_saturation)] = Sw;

         auto pc = capillary_pressure_model.template value<double>(variables,
                                                                  pos, t, dt);

        variables[static_cast<int>(MPL::Variable::capillary_pressure)] = pc;

        _capillary_pressure[ip] = pc;

        _saturation[ip] = Sw;
        _pressure_wetting[ip] = pg_int_pt - pc;
        _liquid_molar_fraction_air[ip] = x_air_wet;
        _liquid_molar_fraction_contaminant[ip] = x_contaminant_wet;
        // TODO: remaining SVs

        double const contaminant_mol_mass =
            contaminant_vapor.property(MPL::PropertyType::molar_mass)
                .template value<double>(variables, pos, t, dt);
        double const mol_density_water = density_water / water_mol_mass;
        double const mol_density_wet = mol_density_water;
        double const density_air_wet = mol_density_wet * air_mol_mass * x_air_wet;
        double const density_contaminant_wet =
            mol_density_wet * contaminant_mol_mass * x_contaminant_wet;
        double const density_wet = density_air_wet + density_water + density_contaminant_wet;

        double const p_vapor_nonwet =
            _process_data.material->calculateVaporPressureNonwet(pc, T_int_pt,
                                                                 density_water);
        double const x_water_nonwet = p_vapor_nonwet / pg_int_pt * x_water_wet;
        _gas_molar_fraction_water[ip] = x_water_nonwet;

        // TODO: should be able to read H_ref and delta from input file.
        double const henry_air =
            _process_data.material->calculateHenryConstant(T_int_pt, 6.4e-6, 1600);
        double const henry_contaminant = _process_data.material->calculateHenryConstant(
            T_int_pt, 6.2e-4, 4500);

        double const x_air_nonwet =
            mol_density_wet / henry_air / pg_int_pt * x_air_wet;
        double const x_contaminant_nonwet =
            mol_density_wet / henry_contaminant / pg_int_pt * x_contaminant_wet;
        _gas_molar_fraction_contaminant[ip] = x_contaminant_nonwet;
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
        double const density_contaminant_nonwet =
            mol_density_nonwet * contaminant_mol_mass * x_contaminant_nonwet;
        double const density_nonwet =
            density_air_nonwet + density_water_nonwet + density_contaminant_nonwet;

        double const mol_density_tot =
            Sw * mol_density_wet + (1 - Sw) * mol_density_nonwet;
        double const d_mol_density_tot_dpg =
            (mol_density_wet - mol_density_nonwet) * dsw_dpg +
            (1 - Sw) * d_mol_density_nonwet_dpg;
        double const d_mol_density_tot_dXa =
            (mol_density_wet - mol_density_nonwet) * dsw_dXa;
        double const d_mol_density_tot_dXc =
            (mol_density_wet - mol_density_nonwet) * dsw_dXc;
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
        Mwc.noalias() +=
            porosity *
            (mol_density_wet * (x_water_wet * dsw_dXc + Sw * dxwL_dXc) +
             mol_density_nonwet *
                 ((1 - Sw) * dxwG_dXc - x_water_nonwet * dsw_dXc)) *
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
        Mac.noalias() += porosity * Xa_int_pt * d_mol_density_tot_dXc *
                         _ip_data[ip].mass_operator;
        Mat.noalias() += porosity * Xa_int_pt * d_mol_density_tot_dT *
                         _ip_data[ip].mass_operator;
        Mcp.noalias() += porosity * Xc_int_pt * d_mol_density_tot_dpg *
                         _ip_data[ip].mass_operator;
        Mca.noalias() += porosity * Xc_int_pt * d_mol_density_tot_dXa *
                         _ip_data[ip].mass_operator;
        Mcc.noalias() += porosity *
                         (Xc_int_pt * d_mol_density_tot_dXc + mol_density_tot) *
                         _ip_data[ip].mass_operator;
        Mct.noalias() += porosity * Xc_int_pt * d_mol_density_tot_dT *
                         _ip_data[ip].mass_operator;

        auto const k_rel =
            medium.property(MPL::PropertyType::relative_permeability)
                .template value<Eigen::Vector2d>(variables, pos, t, dt);

        auto const k_rel_wet = k_rel[0];
        auto const k_rel_nonwet = k_rel[1];
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

        GlobalDimVectorType const velocity_nonwet =
            -lambda_nonwet * K_int *
            (sm.dNdx * pg_nodal_values);
        /* GlobalDimVectorType const velocity_wet =
            -lambda_wet * K *
            (sm.dNdx * (pg_nodal_values - pc_nodal_values)); */

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
        double const heat_capacity_contaminant =
            contaminant_vapor.property(MPL::PropertyType::heat_capacity)
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
        double const enthalpy_contaminant_nonwet =
            _process_data.material->getContaminantEnthalpySimple(
                T_int_pt, heat_capacity_contaminant, contaminant_mol_mass, pg_int_pt);
        double const mol_mass_nonwet =
            x_water_nonwet * water_mol_mass + x_air_nonwet * air_mol_mass +
            x_contaminant_nonwet * contaminant_mol_mass;
        double const X_water_nonwet =
            x_water_nonwet * water_mol_mass / mol_mass_nonwet;
        double const X_air_nonwet =
            x_air_nonwet * air_mol_mass / mol_mass_nonwet;
        double const X_contaminant_nonwet = 1 - X_water_nonwet - X_air_nonwet;
        double const enthalpy_nonwet =
            X_water_nonwet * enthalpy_water_nonwet +
            X_air_nonwet * enthalpy_air_nonwet +
            X_contaminant_nonwet * enthalpy_contaminant_nonwet;
        double const internal_energy_nonwet =
            enthalpy_nonwet - pg_int_pt / density_nonwet;

        double const enthalpy_wet =
            _process_data.material->getLiquidWaterEnthalpySimple(
                T_int_pt, heat_capacity_water, _pressure_wetting[ip]);
        double const internal_energy_wet = enthalpy_wet;

        laplace_operator.noalias() = sm.dNdx.transpose() * K_int *
                                     sm.dNdx * _ip_data[ip].integration_weight;

        // Assemble energy equation M matrix

        Mep.noalias() +=
            porosity *
            ((density_wet * internal_energy_wet -
              density_nonwet * enthalpy_nonwet + pg_int_pt) *
                 dsw_dpg +
             (1 - Sw) *
                 ((x_water_nonwet * water_mol_mass * enthalpy_water_nonwet +
                   x_air_nonwet * air_mol_mass * enthalpy_air_nonwet +
                   x_contaminant_nonwet * contaminant_mol_mass *
                       enthalpy_contaminant_nonwet +
                   pg_int_pt *
                       (dxwG_dpg * water_mol_mass *
                            enthalpy_water_nonwet +
                        dxaG_dpg * air_mol_mass * enthalpy_air_nonwet +
                        dxcG_dpg * contaminant_mol_mass *
                            enthalpy_contaminant_nonwet)) /
                  IdealGasConstant / T_int_pt - 1)) *
            _ip_data[ip].mass_operator;
        Mea.noalias() +=
            porosity * (density_wet * internal_energy_wet * dsw_dXa +
             (1 - Sw) * mol_density_nonwet *
                 (dxwG_dXa * water_mol_mass * enthalpy_water_nonwet +
                  dxaG_dXa * air_mol_mass * enthalpy_air_nonwet +
                  dxcG_dXa * contaminant_mol_mass *
                      enthalpy_contaminant_nonwet)) *
                         _ip_data[ip].mass_operator;
        Mec.noalias() += porosity *
            (density_wet * internal_energy_wet * dsw_dXc +
             (1 - Sw) * mol_density_nonwet *
                 (dxwG_dXc * water_mol_mass * enthalpy_water_nonwet +
                  dxaG_dXc * air_mol_mass * enthalpy_air_nonwet +
                  dxcG_dXc * contaminant_mol_mass *
                      enthalpy_contaminant_nonwet)) *
                         _ip_data[ip].mass_operator;
        Met.noalias() += porosity * (density_wet * (internal_energy_wet * dsw_dT +
                                       Sw * heat_capacity_water) +
                 (1 - Sw) * mol_density_nonwet *
                     (-(x_water_nonwet * water_mol_mass *
                            enthalpy_water_nonwet +
                        x_air_nonwet * air_mol_mass * enthalpy_air_nonwet +
                        x_contaminant_nonwet * contaminant_mol_mass *
                            enthalpy_contaminant_nonwet) /
                      T_int_pt + water_mol_mass *
                     (heat_capacity_water * x_water_nonwet +
                      enthalpy_water_nonwet * dxwG_dT) +
                 (air_mol_mass * heat_capacity_air + IdealGasConstant) * x_air_nonwet +
                 air_mol_mass * enthalpy_air_nonwet * dxaG_dT +
                 (contaminant_mol_mass * heat_capacity_contaminant +
                  IdealGasConstant) *
                     x_contaminant_nonwet +
                 contaminant_mol_mass * enthalpy_contaminant_nonwet * dxcG_dT)) *
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
        Kwc.noalias() +=
            (-mol_density_wet * x_water_wet * lambda_wet * dPC_dSw * dsw_dXc) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxwL_dXc +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxwG_dXc) *
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
        Kac.noalias() +=
            (-mol_density_wet * x_air_wet * lambda_wet * dPC_dSw * dsw_dXc) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxaL_dXc +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxaG_dXc) *
                _ip_data[ip].diffusion_operator;
        Kat.noalias() +=
            (-mol_density_wet * x_air_wet * lambda_wet * dPC_dSw * dsw_dT) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxaL_dT +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxaG_dT) *
                _ip_data[ip].diffusion_operator;
        Kcp.noalias() +=
            (mol_density_wet * x_contaminant_wet * lambda_wet *
                 (1 - dPC_dSw * dsw_dpg) +
             mol_density_nonwet * x_contaminant_nonwet * lambda_nonwet) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxcL_dpg +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxcG_dpg) *
                _ip_data[ip].diffusion_operator;
        Kca.noalias() +=
            (-mol_density_wet * x_contaminant_wet * lambda_wet * dPC_dSw *
             dsw_dXa) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxcL_dXa +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxcG_dXa) *
                _ip_data[ip].diffusion_operator;
        Kcc.noalias() +=
            (-mol_density_wet * x_contaminant_wet * lambda_wet * dPC_dSw *
             dsw_dXc) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxcL_dXc +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxcG_dXc) *
                _ip_data[ip].diffusion_operator;
        Kct.noalias() +=
            (-mol_density_wet * x_contaminant_wet * lambda_wet * dPC_dSw *
             dsw_dT) *
                laplace_operator +
            porosity *
                (Sw * mol_density_wet * diffusion_coeff_wet * dxcL_dT +
                 (1 - Sw) * mol_density_nonwet * diffusion_coeff_nonwet *
                     dxcG_dT) *
                _ip_data[ip].diffusion_operator;
        Kep.noalias() += (lambda_nonwet * density_nonwet * enthalpy_nonwet +
                          lambda_wet * density_wet * enthalpy_wet * (1 - dPC_dSw * dsw_dpg)) *
                         laplace_operator;
        Kea.noalias() +=
            (-lambda_wet * density_wet * enthalpy_wet * dPC_dSw * dsw_dXa) *
            laplace_operator;
        Kec.noalias() +=
            (-lambda_wet * density_wet * enthalpy_wet * dPC_dSw * dsw_dXc) *
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
            Bc.noalias() += (mol_density_wet * x_contaminant_wet * lambda_wet *
                                 density_wet +
                             mol_density_nonwet * x_contaminant_nonwet *
                                 lambda_nonwet * density_nonwet) *
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
                    Mwc(row, row) += Mwc(row, column);
                    Mwc(row, column) = 0.0;
                    Mwt(row, row) += Mwt(row, column);
                    Mwt(row, column) = 0.0;
                    Map(row, row) += Map(row, column);
                    Map(row, column) = 0.0;
                    Maa(row, row) += Maa(row, column);
                    Maa(row, column) = 0.0;
                    Mac(row, row) += Mac(row, column);
                    Mac(row, column) = 0.0;
                    Mat(row, row) += Mat(row, column);
                    Mat(row, column) = 0.0;
                    Mcp(row, row) += Mcp(row, column);
                    Mcp(row, column) = 0.0;
                    Mca(row, row) += Mca(row, column);
                    Mca(row, column) = 0.0;
                    Mcc(row, row) += Mcc(row, column);
                    Mcc(row, column) = 0.0;
                    Mct(row, row) += Mct(row, column);
                    Mct(row, column) = 0.0;
                    Mep(row, row) += Mep(row, column);
                    Mep(row, column) = 0.0;
                    Mea(row, row) += Mea(row, column);
                    Mea(row, column) = 0.0;
                    Mec(row, row) += Mec(row, column);
                    Mec(row, column) = 0.0;
                    Met(row, row) += Met(row, column);
                    Met(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping

}

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
