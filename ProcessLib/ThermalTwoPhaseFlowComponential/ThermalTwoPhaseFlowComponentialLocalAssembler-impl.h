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
#include <iostream>

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

    auto Mwpc =
        local_M.template block<capillary_pressure_size, capillary_pressure_size>(
            capillary_pressure_matrix_index, capillary_pressure_matrix_index);
    auto Mwp =
        local_M.template block<capillary_pressure_size, gas_pressure_size>(
            capillary_pressure_matrix_index, gas_pressure_matrix_index);
    auto Mwx = local_M.template block<capillary_pressure_size,
                                       overall_mol_frac_contaminant_size>(
        capillary_pressure_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Mwt = local_M.template block<capillary_pressure_size, temperature_size>(
        capillary_pressure_matrix_index, temperature_matrix_index);

    auto Mapc =
        local_M.template block<gas_pressure_size, capillary_pressure_size>(
            gas_pressure_matrix_index, capillary_pressure_matrix_index);
    auto Map = local_M.template block<gas_pressure_size,
                                       gas_pressure_size>(
        gas_pressure_matrix_index, gas_pressure_matrix_index);
    auto Max = local_M.template block<gas_pressure_size,
                                       overall_mol_frac_contaminant_size>(
        gas_pressure_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Mat =
        local_M.template block<gas_pressure_size, temperature_size>(
            gas_pressure_matrix_index, temperature_matrix_index);

    auto Mcpc = local_M.template block<overall_mol_frac_contaminant_size,
                                      capillary_pressure_size>(
        overall_mol_frac_contaminant_matrix_index,
        capillary_pressure_matrix_index);
    auto Mcp = local_M.template block<overall_mol_frac_contaminant_size,
                                       gas_pressure_size>(
        overall_mol_frac_contaminant_matrix_index,
        gas_pressure_matrix_index);
    auto Mcx = local_M.template block<overall_mol_frac_contaminant_size,
                                       overall_mol_frac_contaminant_size>(
        overall_mol_frac_contaminant_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Mct = local_M.template block<overall_mol_frac_contaminant_size,
                                      temperature_size>(
        overall_mol_frac_contaminant_matrix_index, temperature_matrix_index);

    auto Mepc = local_M.template block<temperature_size,
                                      capillary_pressure_size>(
        temperature_matrix_index,
        capillary_pressure_matrix_index);
    auto Mep =
        local_M.template block<temperature_size,
                                       gas_pressure_size>(
            temperature_matrix_index,
        gas_pressure_matrix_index);
    auto Mex = local_M.template block<temperature_size,
                                       overall_mol_frac_contaminant_size>(
        temperature_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Met = local_M.template block<temperature_size,
                                      temperature_size>(
        temperature_matrix_index, temperature_matrix_index);
    
    NodalMatrixType laplace_operator =
        NodalMatrixType::Zero(ShapeFunction::NPOINTS, ShapeFunction::NPOINTS);

    auto Kwpc =
        local_K.template block<capillary_pressure_size, capillary_pressure_size>(
            capillary_pressure_matrix_index, capillary_pressure_matrix_index);
    auto Kwp =
        local_K.template block<capillary_pressure_size, gas_pressure_size>(
            capillary_pressure_matrix_index, gas_pressure_matrix_index);
    auto Kwx = local_K.template block<capillary_pressure_size,
                                       overall_mol_frac_contaminant_size>(
        capillary_pressure_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Kwt = local_K.template block<capillary_pressure_size, temperature_size>(
        capillary_pressure_matrix_index, temperature_matrix_index);

    auto Kapc =
        local_K.template block<gas_pressure_size, capillary_pressure_size>(
            gas_pressure_matrix_index, capillary_pressure_matrix_index);
    auto Kap = local_K.template block<gas_pressure_size,
                                       gas_pressure_size>(
        gas_pressure_matrix_index, gas_pressure_matrix_index);
    auto Kax = local_K.template block<gas_pressure_size,
                                       overall_mol_frac_contaminant_size>(
        gas_pressure_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Kat =
        local_K.template block<gas_pressure_size, temperature_size>(
            gas_pressure_matrix_index, temperature_matrix_index);

    auto Kcpc = local_K.template block<overall_mol_frac_contaminant_size,
                                      capillary_pressure_size>(
        overall_mol_frac_contaminant_matrix_index,
        capillary_pressure_matrix_index);
    auto Kcp = local_K.template block<overall_mol_frac_contaminant_size,
                                       gas_pressure_size>(
        overall_mol_frac_contaminant_matrix_index,
        gas_pressure_matrix_index);
    auto Kcx = local_K.template block<overall_mol_frac_contaminant_size,
                                       overall_mol_frac_contaminant_size>(
        overall_mol_frac_contaminant_matrix_index,
        overall_mol_frac_contaminant_matrix_index);
    auto Kct = local_K.template block<overall_mol_frac_contaminant_size,
                                      temperature_size>(
        overall_mol_frac_contaminant_matrix_index, temperature_matrix_index);

    auto Kepc = local_K.template block<temperature_size, capillary_pressure_size>(
        temperature_matrix_index, capillary_pressure_matrix_index);
    auto Kep =
        local_K.template block<temperature_size, gas_pressure_size>(
            temperature_matrix_index, gas_pressure_matrix_index);
    auto Kex = local_K.template block<temperature_size,
                                       overall_mol_frac_contaminant_size>(
        temperature_matrix_index, overall_mol_frac_contaminant_matrix_index);
    auto Ket = local_K.template block<temperature_size, temperature_size>(
        temperature_matrix_index, temperature_matrix_index);

    auto Bw = local_b.template segment<capillary_pressure_size>(
        capillary_pressure_matrix_index);

    auto Ba = local_b.template segment<gas_pressure_size>(
        gas_pressure_matrix_index);

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
    auto const& dissolved_contaminant = liquid_phase.component("contaminant");
    auto const& water_vapor = gas_phase.component("wasser");
    auto const& gaseous_air = gas_phase.component("air");
    auto const& contaminant_vapor = gas_phase.component("contaminant");

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto const& b = _process_data.specific_body_force;
    GlobalDimMatrixType const& I(
        GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

    auto const num_nodes = ShapeFunction::NPOINTS;

    auto const pc_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[0], num_nodes);
    auto const pg_nodal_values =
        Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        pos.setIntegrationPoint(ip); // TODO: necessary?
        auto const& sm = _shape_matrices[ip];

        double pc_int_pt = 0.;
        double pg_int_pt = 0.;
        double Xc_int_pt =
            0.;  // total molar fraction of the light component contaminant
        double T_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pc_int_pt, pg_int_pt,
                                         Xc_int_pt,
                                         T_int_pt);

        MPL::VariableArray variables;

        variables[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T_int_pt;
        variables[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = pc_int_pt;
        variables[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            pg_int_pt;

        double const density_water =
            liquid_water.property(MPL::PropertyType::density)
                .template value<double>(variables, pos, t, dt);

        double const pl = pg_int_pt - pc_int_pt;
        _pressure_wetting[ip] = pl;

        auto& saturation_model =
            medium.property(MPL::PropertyType::saturation);
        double const Sw =
            saturation_model.template value<double>(
            variables, pos, t, dt);
        _saturation[ip] = Sw;
        variables[static_cast<int>(
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
                .template value<double>(variables, pos, t, dt);

        double const p_vapor_nonwet =
            _process_data.material->calculateVaporPressureNonwet(pc_int_pt, T_int_pt,
                                                                 density_water);

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

        _gas_molar_fraction_water[ip] = x_water_nonwet;
        _liquid_molar_fraction_contaminant[ip] = x_contaminant_wet;
        _gas_molar_fraction_contaminant[ip] = x_contaminant_nonwet;

        auto dSw_dPc = saturation_model.template dValue<double>(
            variables, MPL::Variable::capillary_pressure, pos, t, dt);

        double const d_mol_density_nonwet_dpg =
            1 / ideal_gas_constant_times_T_int_pt;
        double const d_mol_density_nonwet_dT = -mol_density_nonwet / T_int_pt;
        double const d_mol_density_tot_dpc =
            (mol_density_wet - mol_density_nonwet) * dSw_dPc;
        double const d_mol_density_tot_dpg =
            (1 - Sw) * d_mol_density_nonwet_dpg;
        double const d_mol_density_tot_dT = (1 - Sw) * d_mol_density_nonwet_dT;

        double const d_kc_dpg = henry_contaminant / mol_density_wet;
        double d_henry_contaminant_dT =
            dissolved_contaminant.property(MPL::PropertyType::henry_constant)
                .template dValue<double>(variables, MPL::Variable::temperature,
                                         pos, t, dt);

        double const d_kc_dT =
            pg_int_pt * d_henry_contaminant_dT / mol_density_wet;

        double const dxcG_dpc =
            Xc_int_pt * (d_mol_density_tot_dpc / ntot_c -
                         (mol_density_wet * k_c - mol_density_nonwet) *
                             dSw_dPc * mol_density_tot / ntot_c / ntot_c);
        double const dxcG_dpg =
            Xc_int_pt * (d_mol_density_tot_dpg / ntot_c -
                         (Sw * mol_density_wet * d_kc_dpg +
                          (1 - Sw) * d_mol_density_nonwet_dpg) *
                             mol_density_tot / ntot_c / ntot_c);
        double const dxcG_dXc = mol_density_tot / ntot_c;
        double const dxcG_dT =
            Xc_int_pt * (d_mol_density_tot_dT / ntot_c -
                         (Sw * mol_density_wet * d_kc_dT +
                          (1 - Sw) * d_mol_density_nonwet_dT) *
                             mol_density_tot / ntot_c / ntot_c);

        double const dxcL_dpc = k_c * dxcG_dpc;
        double const dxcL_dpg =
            k_c * dxcG_dpg + d_kc_dpg * x_contaminant_nonwet;
        double const dxcL_dXc = k_c * dxcG_dXc;
        double const dxcL_dT =
            k_c * dxcG_dT + d_kc_dT * x_contaminant_nonwet;

        double const dxwL_dpc = -dxcL_dpc;
        double const dxwL_dpg = -dxcL_dpg;
        double const dxwL_dXc = -dxcL_dXc;
        double const dxwL_dT = -dxcL_dT;

        double const d_p_vapor_nonwet_d_pc =
            _process_data.material->calculateDerivativedPgwdPC(
                pc_int_pt, T_int_pt, density_water);
        double const d_p_vapor_nonwet_d_T =
            _process_data.material->calculateDerivativedPgwdT(
                pc_int_pt, T_int_pt, density_water);

        double const dxwG_dpc =
            (d_p_vapor_nonwet_d_pc * x_water_wet + p_vapor_nonwet * dxwL_dpc) /
            pg_int_pt;
        double const dxwG_dpg =
            p_vapor_nonwet *
            (dxwL_dpg / pg_int_pt - x_water_wet / pg_int_pt / pg_int_pt);
        double const dxwG_dXc = dxwL_dXc / k_w;
        double const dxwG_dT =
            (d_p_vapor_nonwet_d_T * x_water_wet + p_vapor_nonwet * dxwL_dT) /
            pg_int_pt;

        double const dxaG_dpc = -dxwG_dpc - dxcG_dpc;
        double const dxaG_dpg = -dxwG_dpg - dxcG_dpg;
        double const dxaG_dXc = -dxwG_dXc - dxcG_dXc;
        double const dxaG_dT = -dxwG_dT - dxcG_dT;

        auto const porosity =
            medium.property(MPL::PropertyType::porosity)
                .template value<double>(variables, pos, t, dt);

        // Assemble mass blance equations
        Mwpc.noalias() +=
            porosity *
            (mol_density_wet * (x_water_wet * dSw_dPc + Sw * dxwL_dpc) +
             mol_density_nonwet *
                 ((1 - Sw) * dxwG_dpc - x_water_nonwet * dSw_dPc)) *
            _ip_data[ip].mass_operator;
        Mwp.noalias() += porosity *
                         (mol_density_wet * Sw * dxwL_dpg +
                          (1 - Sw) * x_water_nonwet * d_mol_density_nonwet_dpg +
                          mol_density_nonwet * (1 - Sw) * dxwG_dpg) *
                         _ip_data[ip].mass_operator;
        Mwx.noalias() += porosity *
                         (mol_density_wet * Sw * dxwL_dXc +
                          mol_density_nonwet * (1 - Sw) * dxwG_dXc) *
                         _ip_data[ip].mass_operator;
        Mwt.noalias() += porosity *
                         (mol_density_wet * Sw * dxwL_dT +
                          (1 - Sw) * x_water_nonwet * d_mol_density_nonwet_dT +
                          mol_density_nonwet * (1 - Sw) * dxwG_dT) *
                         _ip_data[ip].mass_operator;
        Mapc.noalias() += porosity * mol_density_nonwet *
                          ((1 - Sw) * dxaG_dpc - x_air_nonwet * dSw_dPc) *
                          _ip_data[ip].mass_operator;
        Map.noalias() += porosity * (1 - Sw) *
                         (x_air_nonwet * d_mol_density_nonwet_dpg +
                          mol_density_nonwet * dxaG_dpg) *
                         _ip_data[ip].mass_operator;
        Max.noalias() += porosity * mol_density_nonwet * (1 - Sw) * dxaG_dXc *
                         _ip_data[ip].mass_operator;
        Mat.noalias() += porosity * (1 - Sw) *
                         (x_air_nonwet * d_mol_density_nonwet_dT +
                          mol_density_nonwet * dxaG_dT) *
                         _ip_data[ip].mass_operator;
        Mcpc.noalias() +=
            porosity *
            (mol_density_wet * (x_contaminant_wet * dSw_dPc + Sw * dxcL_dpc) +
             mol_density_nonwet *
                 ((1 - Sw) * dxcG_dpc - x_contaminant_nonwet * dSw_dPc)) *
            _ip_data[ip].mass_operator;
        Mcp.noalias() += porosity *
                         (mol_density_wet * Sw * dxcL_dpg +
                          (1 - Sw) * x_contaminant_nonwet * d_mol_density_nonwet_dpg +
                          mol_density_nonwet * (1 - Sw) * dxcG_dpg) *
                         _ip_data[ip].mass_operator;
        Mcx.noalias() += porosity *
                         (mol_density_wet * Sw * dxcL_dXc +
                          mol_density_nonwet * (1 - Sw) * dxcG_dXc) *
                         _ip_data[ip].mass_operator;
        Mct.noalias() += porosity *
                         (mol_density_wet * Sw * dxcL_dT +
                          (1 - Sw) * x_contaminant_nonwet * d_mol_density_nonwet_dT +
                          mol_density_nonwet * (1 - Sw) * dxcG_dT) *
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

        // diffusion coefficients
        double diffusion_coeff_wet =
            liquid_phase.property(MPL::PropertyType::molecular_diffusion)
                .template value<double>(variables, pos, t, dt);
        double diffusion_coeff_vapor_nonwet =
            water_vapor.property(MPL::PropertyType::molecular_diffusion)
                .template value<double>(variables, pos, t, dt);
        double diffusion_coeff_contaminant_nonwet =
            contaminant_vapor.property(MPL::PropertyType::molecular_diffusion)
                .template value<double>(variables, pos, t, dt);
        // multiply tortuosity factor
        diffusion_coeff_wet *=
            std::pow(Sw, 7. / 3) * std::pow(porosity, 1. / 3);
        diffusion_coeff_vapor_nonwet *=
            std::pow(1 - Sw, 7. / 3) * std::pow(porosity, 1. / 3);
        diffusion_coeff_contaminant_nonwet *=
            std::pow(1 - Sw, 7. / 3) * std::pow(porosity, 1. / 3);

        auto const solute_dispersivity_transverse =
            medium.template value<double>(
                MaterialPropertyLib::PropertyType::transversal_dispersivity);
        auto const solute_dispersivity_longitudinal =
            medium.template value<double>(
                MaterialPropertyLib::PropertyType::longitudinal_dispersivity);

        double const contaminant_mol_mass =
            contaminant_vapor.property(MPL::PropertyType::molar_mass)
                .template value<double>(variables, pos, t, dt);
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

        GlobalDimVectorType const velocity_wet =
            _process_data.has_gravity
                ? GlobalDimVectorType(
                      -lambda_wet * K_int *
                      (sm.dNdx * (pg_nodal_values - pc_nodal_values) -
                       density_wet * b))
                : GlobalDimVectorType(-lambda_wet * K_int * sm.dNdx *
                                      (pg_nodal_values - pc_nodal_values));

        double const velocity_wet_magnitude = velocity_wet.norm();
        GlobalDimMatrixType const hydrodynamic_dispersion =
            velocity_wet_magnitude != 0.0
                ? GlobalDimMatrixType((porosity * Sw * diffusion_coeff_wet +
                                       solute_dispersivity_transverse *
                                           velocity_wet_magnitude) *
                                          I +
                                      (solute_dispersivity_longitudinal -
                                       solute_dispersivity_transverse) /
                                          velocity_wet_magnitude *
                                          velocity_wet *
                                          velocity_wet.transpose())
                : GlobalDimMatrixType((porosity * Sw * diffusion_coeff_wet +
                                       solute_dispersivity_transverse *
                                           velocity_wet_magnitude) *
                                      I);
        
        laplace_operator.noalias() = sm.dNdx.transpose() * K_int * sm.dNdx *
                                     _ip_data[ip].integration_weight;

        auto const dispersion_operator = sm.dNdx.transpose() *
                                         hydrodynamic_dispersion * sm.dNdx *
                                         _ip_data[ip].integration_weight;

        Kwpc.noalias() +=
            (-mol_density_wet * x_water_wet * lambda_wet) * laplace_operator +
            mol_density_wet * dispersion_operator * dxwL_dpc +
            porosity * (1 - Sw) * mol_density_nonwet *
                diffusion_coeff_vapor_nonwet * dxwG_dpc *
                _ip_data[ip].diffusion_operator;
        Kwp.noalias() += (mol_density_wet * x_water_wet * lambda_wet +
                          mol_density_nonwet * x_water_nonwet * lambda_nonwet) *
                             laplace_operator +
                         mol_density_wet * dispersion_operator * dxwL_dpg +
                         porosity * (1 - Sw) * mol_density_nonwet *
                             diffusion_coeff_vapor_nonwet * dxwG_dpg *
                             _ip_data[ip].diffusion_operator;
        Kwx.noalias() += mol_density_wet * dispersion_operator * dxwL_dXc +
                         porosity * (1 - Sw) * mol_density_nonwet *
                             diffusion_coeff_vapor_nonwet * dxwG_dXc *
                             _ip_data[ip].diffusion_operator;
        Kwt.noalias() += mol_density_wet * dispersion_operator * dxwL_dT +
                         porosity * (1 - Sw) * mol_density_nonwet *
                             diffusion_coeff_vapor_nonwet * dxwG_dT *
                             _ip_data[ip].diffusion_operator;
        Kapc.noalias() += porosity * (1 - Sw) * mol_density_nonwet *
                          diffusion_coeff_vapor_nonwet * dxaG_dpc *
                          _ip_data[ip].diffusion_operator;
        Kap.noalias() += (mol_density_nonwet * x_air_nonwet * lambda_nonwet) *
                             laplace_operator +
                         porosity * (1 - Sw) * mol_density_nonwet *
                             diffusion_coeff_vapor_nonwet * dxaG_dpg *
                             _ip_data[ip].diffusion_operator;
        Kax.noalias() += porosity * (1 - Sw) * mol_density_nonwet *
                         diffusion_coeff_vapor_nonwet * dxaG_dXc *
                         _ip_data[ip].diffusion_operator;
        Kat.noalias() += porosity * (1 - Sw) * mol_density_nonwet *
                         diffusion_coeff_vapor_nonwet * dxaG_dT *
                         _ip_data[ip].diffusion_operator;
        Kcpc.noalias() += (-mol_density_wet * x_contaminant_wet * lambda_wet) *
                              laplace_operator +
                          mol_density_wet * dispersion_operator * dxcL_dpc +
                          porosity * (1 - Sw) * mol_density_nonwet *
                              diffusion_coeff_contaminant_nonwet * dxcG_dpc *
                              _ip_data[ip].diffusion_operator;
        Kcp.noalias() +=
            (mol_density_wet * x_contaminant_wet * lambda_wet +
             mol_density_nonwet * x_contaminant_nonwet * lambda_nonwet) *
                laplace_operator +
            mol_density_wet * dispersion_operator * dxcL_dpg +
            porosity * (1 - Sw) * mol_density_nonwet *
                diffusion_coeff_contaminant_nonwet * dxcG_dpg *
                _ip_data[ip].diffusion_operator;
        Kcx.noalias() += mol_density_wet * dispersion_operator * dxcL_dXc +
                         porosity * (1 - Sw) * mol_density_nonwet *
                             diffusion_coeff_contaminant_nonwet * dxcG_dXc *
                             _ip_data[ip].diffusion_operator;
        Kct.noalias() += mol_density_wet * dispersion_operator * dxcL_dT +
                         porosity * (1 - Sw) * mol_density_nonwet *
                             diffusion_coeff_contaminant_nonwet * dxcG_dT *
                             _ip_data[ip].diffusion_operator;

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
                T_int_pt, heat_capacity_contaminant, contaminant_mol_mass,
                pg_int_pt);

        double const enthalpy_nonwet =
            X_water_nonwet * enthalpy_water_nonwet +
            X_air_nonwet * enthalpy_air_nonwet +
            X_contaminant_nonwet * enthalpy_contaminant_nonwet;
        double const internal_energy_nonwet =
            enthalpy_nonwet - pg_int_pt / density_nonwet;

        double const enthalpy_wet =
            _process_data.material->getLiquidWaterEnthalpySimple(
                T_int_pt, heat_capacity_water, pl);
        double const internal_energy_wet = enthalpy_wet;
        double const volumetric_enthalpy_nonwet =
            density_nonwet * enthalpy_nonwet;

        double const d_enthalpy_water_nonwet_dT = heat_capacity_water;
        double const d_enthalpy_air_nonwet_dT =
            heat_capacity_air + IdealGasConstant / air_mol_mass;
        double const d_enthalpy_contaminant_nonwet_dT =
            heat_capacity_contaminant + IdealGasConstant / contaminant_mol_mass;

        double const d_vol_enthalpy_nonwet_dpc =
            mol_density_nonwet *
            (dxwG_dpc * water_mol_mass * enthalpy_water_nonwet +
             dxaG_dpc * air_mol_mass * enthalpy_air_nonwet +
             dxcG_dpc * contaminant_mol_mass * enthalpy_contaminant_nonwet);
        double const d_vol_enthalpy_nonwet_dpg =
            d_mol_density_nonwet_dpg *
                (x_water_nonwet * water_mol_mass * enthalpy_water_nonwet +
                 x_air_nonwet * air_mol_mass * enthalpy_air_nonwet +
                 x_contaminant_nonwet * contaminant_mol_mass *
                     enthalpy_contaminant_nonwet) +
            mol_density_nonwet *
                (dxwG_dpg * water_mol_mass * enthalpy_water_nonwet +
                 dxaG_dpg * air_mol_mass * enthalpy_air_nonwet +
                 dxcG_dpg * contaminant_mol_mass * enthalpy_contaminant_nonwet);
        double const d_vol_enthalpy_nonwet_dXc =
            mol_density_nonwet *
            (dxwG_dXc * water_mol_mass * enthalpy_water_nonwet +
             dxaG_dXc * air_mol_mass * enthalpy_air_nonwet +
             dxcG_dXc * contaminant_mol_mass * enthalpy_contaminant_nonwet);
        double const d_vol_enthalpy_nonwet_dT =
            d_mol_density_nonwet_dT *
                (x_water_nonwet * water_mol_mass * enthalpy_water_nonwet +
                 x_air_nonwet * air_mol_mass * enthalpy_air_nonwet +
                 x_contaminant_nonwet * contaminant_mol_mass *
                     enthalpy_contaminant_nonwet) +
            mol_density_nonwet *
                (water_mol_mass *
                     (dxwG_dT * enthalpy_water_nonwet +
                      x_water_nonwet * d_enthalpy_water_nonwet_dT) +
                 air_mol_mass * (dxaG_dT * enthalpy_air_nonwet +
                                 x_air_nonwet * d_enthalpy_air_nonwet_dT) +
                 contaminant_mol_mass *
                     (dxcG_dT * enthalpy_contaminant_nonwet +
                      x_contaminant_nonwet * d_enthalpy_contaminant_nonwet_dT));
           
        // Assemble energy blance equation
        Mepc.noalias() += porosity *
                              (density_wet * internal_energy_wet -
                               density_nonwet * internal_energy_nonwet) *
                              dSw_dPc * _ip_data[ip].mass_operator +
                          porosity * (1 - Sw) * d_vol_enthalpy_nonwet_dpc *
                              _ip_data[ip].mass_operator;
        Mep.noalias() += porosity * (1 - Sw) * (d_vol_enthalpy_nonwet_dpg - 1) *
                         _ip_data[ip].mass_operator;
        Mex.noalias() += porosity * (1 - Sw) * d_vol_enthalpy_nonwet_dXc *
                         _ip_data[ip].mass_operator;
        Met.noalias() += ((1 - porosity) * density_solid * heat_capacity_solid +
                          porosity * ((1 - Sw) * d_vol_enthalpy_nonwet_dT +
                                      Sw * density_wet * heat_capacity_water)) *
                         _ip_data[ip].mass_operator;

        Kepc.noalias() +=
            -lambda_wet * enthalpy_wet * density_wet * laplace_operator;
        Kep.noalias() += (lambda_nonwet * density_nonwet * enthalpy_nonwet +
                          lambda_wet * density_wet * enthalpy_wet) *
                         laplace_operator;
        Ket.noalias() += sm.dNdx.transpose() * effective_thermal_conductivity *
                         sm.dNdx * _ip_data[ip].integration_weight;

        if (_process_data.has_gravity)
        {
            NodalVectorType gravity_operator = sm.dNdx.transpose() * K_int * b *
                                               _ip_data[ip].integration_weight;
            Bw.noalias() +=
                (mol_density_wet * x_water_wet * lambda_wet * density_wet +
                 mol_density_nonwet * x_water_nonwet * lambda_nonwet *
                     density_nonwet) *
                gravity_operator;
            Ba.noalias() += (mol_density_nonwet * x_air_nonwet * lambda_nonwet *
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
                    Mwpc(row, row) += Mwpc(row, column);
                    Mwpc(row, column) = 0.0;
                    Mwp(row, row) += Mwp(row, column);
                    Mwp(row, column) = 0.0;
                    Mwx(row, row) += Mwx(row, column);
                    Mwx(row, column) = 0.0;
                    Mwt(row, row) += Mwt(row, column);
                    Mwt(row, column) = 0.0;
                    Mapc(row, row) += Mapc(row, column);
                    Mapc(row, column) = 0.0;
                    Map(row, row) += Map(row, column);
                    Map(row, column) = 0.0;
                    Max(row, row) += Max(row, column);
                    Max(row, column) = 0.0;
                    Mat(row, row) += Mat(row, column);
                    Mat(row, column) = 0.0;
                    Mcpc(row, row) += Mcpc(row, column);
                    Mcpc(row, column) = 0.0;
                    Mcp(row, row) += Mcp(row, column);
                    Mcp(row, column) = 0.0;
                    Mcx(row, row) += Mcx(row, column);
                    Mcx(row, column) = 0.0;
                    Mct(row, row) += Mct(row, column);
                    Mct(row, column) = 0.0;
                    Mepc(row, row) += Mepc(row, column);
                    Mepc(row, column) = 0.0;
                    Mep(row, row) += Mep(row, column);
                    Mep(row, column) = 0.0;
                    Mex(row, row) += Mex(row, column);
                    Mex(row, column) = 0.0;
                    Met(row, row) += Met(row, column);
                    Met(row, column) = 0.0;
                }
            }
        }
    }  // end of mass-lumping

}

}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
