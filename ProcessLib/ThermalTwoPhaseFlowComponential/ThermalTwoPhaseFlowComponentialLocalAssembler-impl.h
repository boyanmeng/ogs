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

        double pc_int_pt = 0.;
        double pg_int_pt = 0.;
        double Xc_int_pt =
            0.;  // total molar fraction of the light component contaminant
        double T_int_pt = 0.;
        NumLib::shapeFunctionInterpolate(local_x, sm.N, pc_int_pt, pg_int_pt,
                                         Xc_int_pt,
                                         T_int_pt);


        
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
