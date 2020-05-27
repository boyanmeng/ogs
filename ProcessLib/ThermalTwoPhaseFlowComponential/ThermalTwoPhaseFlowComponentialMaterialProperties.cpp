/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermalTwoPhaseFlowComponentialMaterialProperties.h"
#include <utility>  // ?
#include <cmath>
#include "BaseLib/Logging.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"
#include "MaterialLib/MPL/Property.h"

namespace ProcessLib
{
using MaterialLib::PhysicalConstant::IdealGasConstant;
using MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;

namespace ThermalTwoPhaseFlowComponential
{
ThermalTwoPhaseFlowComponentialMaterialProperties::ThermalTwoPhaseFlowComponentialMaterialProperties(
        std::unique_ptr<MaterialLib::Fluid::WaterVaporProperties>&&
            water_vapor_properties)
    : _water_vapor_properties(std::move(water_vapor_properties))
{
    DBUG("Create material properties for ThermalTwoPhaseFlowComponential process model.");
}

double
ThermalTwoPhaseFlowComponentialMaterialProperties::calculateHenryConstant(
    const double T, const double H_ref, const double delta) const
{
    double const T_ref = 298.15;
    return H_ref * std::exp((1 / T - 1 / T_ref) * delta);
}

double
ThermalTwoPhaseFlowComponentialMaterialProperties::calculateDerivativedHdT(
    const double T, const double H_ref, const double delta) const
{
    double const henry = calculateHenryConstant(T, H_ref, delta);
    return 0;
}

  double
  ThermalTwoPhaseFlowComponentialMaterialProperties::calculateSaturatedVaporPressure(
      const double T) const
  {
      return _water_vapor_properties->calculateSaturatedVaporPressure(T);
  }
  double
  ThermalTwoPhaseFlowComponentialMaterialProperties::calculateVaporPressureNonwet(
      const double pc, const double T, const double mass_density_water) const
  {
      return _water_vapor_properties->calculateVaporPressureNonwet(
          pc, T, mass_density_water);
  }
  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateDerivativedPsatdT(
      const double T) const
  {
      return _water_vapor_properties->calculateDerivativedPsatdT(T);
  }
  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateDerivativedPgwdT(
      const double pc, const double T, const double mass_density_water) const
  {
      return _water_vapor_properties->calculateDerivativedPgwdT(
          pc, T, mass_density_water);
  }
  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateDerivativedPgwdPC(
      const double pc, const double T, const double mass_density_water) const
  {
      return _water_vapor_properties->calculateDerivativedPgwdPC(
          pc, T, mass_density_water);
  }
  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculatedDensityNonwetdT(
      const double p_air_nonwet, const double p_vapor_nonwet, const double pc,
      const double T, const double mass_density_water) const
  {
      return _water_vapor_properties->calculatedDensityNonwetdT(
          p_air_nonwet, p_vapor_nonwet, pc, T, mass_density_water);
  }
  
  double ThermalTwoPhaseFlowComponentialMaterialProperties::getWaterVaporEnthalpySimple(
      const double temperature, const double heat_capacity_water_vapor,
      const double pg, const double latent_heat_evaporation) const
  {
      return _water_vapor_properties->getWaterVaporEnthalpySimple(
          temperature, heat_capacity_water_vapor, pg, latent_heat_evaporation);
  }
  
  double ThermalTwoPhaseFlowComponentialMaterialProperties::getAirEnthalpySimple(
      const double temperature,
      const double heat_capacity_dry_air,
      const double /*pg*/) const
  {
      return heat_capacity_dry_air * (temperature - CelsiusZeroInKelvin) +
             IdealGasConstant * temperature /
                 _air_mol_mass;
  }

  double ThermalTwoPhaseFlowComponentialMaterialProperties::
      getContaminantEnthalpySimple(const double temperature,
                                   const double contaminant_heat_capacity,
                                   const double contaminant_molar_mass,
                                   const double) const
  {
      return contaminant_heat_capacity * (temperature - CelsiusZeroInKelvin) +
             IdealGasConstant * temperature /
                 contaminant_molar_mass;
  }
  
  double
  ThermalTwoPhaseFlowComponentialMaterialProperties::getLiquidWaterEnthalpySimple(
      const double temperature,
      const double heat_capacity_liquid_water,
      const double /*pl*/) const
  {
      return heat_capacity_liquid_water * (temperature - CelsiusZeroInKelvin);
  }

  bool ThermalTwoPhaseFlowComponentialMaterialProperties::
      computeConstitutiveRelation(
          double const t, double const dt,
          ParameterLib::SpatialPosition const& x_position,
          MaterialPropertyLib::Property const& pc_model, double const rho_w,
          double const pg,
                                  double const Xa,
                                  double const Xc, double const T, double& Sw, double& x_w_L, double& x_a_L, double& x_c_L,
          double& dsw_dpg, double& dxwG_dpg, double& dxaG_dpg, double& dxcG_dpg,
          double& dsw_dXa, double& dxwG_dXa, double& dxaG_dXa, double& dxcG_dXa,
          double& dsw_dXc, double& dxwG_dXc, double& dxaG_dXc, double& dxcG_dXc,
          double& dsw_dT, double& dxwG_dT, double& dxaG_dT, double& dxcG_dT,
          double& dxwL_dpg, double& dxaL_dpg, double& dxcL_dpg,
          double& dxwL_dXa, double& dxaL_dXa, double& dxcL_dXa,
          double& dxwL_dXc, double& dxaL_dXc, double& dxcL_dXc, double& dxwL_dT,
          double& dxaL_dT, double& dxcL_dT)
  {
      {  // Local Newton solver
          using LocalJacobianMatrix =
              Eigen::Matrix<double, 4, 4, Eigen::RowMajor>;
          using LocalResidualVector = Eigen::Matrix<double, 4, 1>;
          using LocalUnknownVector = Eigen::Matrix<double, 4, 1>;
          LocalJacobianMatrix J_loc;

          Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(4);
          auto const update_residual = [&](LocalResidualVector& residual) {
              calculateResidual(t, dt, x_position, pc_model, rho_w, pg, Xa, Xc,
                                T, Sw,
                                x_w_L, x_a_L, x_c_L, residual);
          };

          auto const update_jacobian = [&](LocalJacobianMatrix& jacobian) {
              calculateJacobian(t, dt, x_position, pc_model, rho_w, pg, Xa, Xc, T, jacobian, Sw, x_w_L, x_a_L,
                                x_c_L);  // for solution dependent Jacobians
          };

          auto const update_solution =
              [&](LocalUnknownVector const& increment) {
                  // increment solution vectors
                  Sw += increment[0];
                  x_w_L += increment[1];
                  x_a_L += increment[2];
                  x_c_L += increment[3];
              };

          // TODO Make the following choice of maximum iterations and
          // convergence criteria available from the input file configuration.
          // See Ehlers material model implementation for the example.
          const int maximum_iterations(20);
          const double tolerance(1.e-14);

          auto newton_solver = NumLib::NewtonRaphson<
              decltype(linear_solver), LocalJacobianMatrix,
              decltype(update_jacobian), LocalResidualVector,
              decltype(update_residual), decltype(update_solution)>(
              linear_solver, update_jacobian, update_residual, update_solution,
              {maximum_iterations, tolerance});

          auto const success_iterations = newton_solver.solve(J_loc);

          if (!success_iterations)
          {
              return false;
          }
      }
      calculateDerivatives(t, dt, x_position, pc_model, rho_w, pg, Xa, Xc, T,
                           Sw, x_w_L, x_a_L, x_c_L, dsw_dpg, dxwG_dpg, dxaG_dpg, dxcG_dpg, dsw_dXa,
          dxwG_dXa, dxaG_dXa, dxcG_dXa, dsw_dXc,
          dxwG_dXc, dxaG_dXc, dxcG_dXc, dsw_dT,
          dxwG_dT, dxaG_dT, dxcG_dT, dxwL_dpg, dxaL_dpg,
          dxcL_dpg, dxwL_dXa, dxaL_dXa, dxcL_dXa,
          dxwL_dXc, dxaL_dXc, dxcL_dXc, dxwL_dT,
          dxaL_dT, dxcL_dT);
      return true;
  }

  void ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResidual(
      double const t, double const dt,
      ParameterLib::SpatialPosition const& x_position,
      MaterialPropertyLib::Property const& pc_model, double const rho_w,
      double const pg, double const Xa, double const Xc, double const T,
      double Sw, double x_w_L, double x_a_L, double x_c_L, ResidualVector& res)
  {
      MaterialPropertyLib::VariableArray vars;
      vars[static_cast<int>(MaterialPropertyLib::Variable::liquid_saturation)] =
          Sw;
      double const pc = pc_model.template value<double>(vars, x_position, t, dt);
      // use water density for simplicity
      double const p_vap =
          calculateVaporPressureNonwet(pc, T, rho_w);
      double const x_w_G = p_vap / pg * x_w_L;
      double const H_a = 8.5906e-6;
      double const H_c = calculateHenryConstant(T, 6.2e-4, 4500);
      double const N_w = rho_w / _water_mol_mass;
      double const x_a_G = N_w / H_a / pg * x_a_L;
      double const x_c_G = N_w / H_c / pg * x_c_L;
      // calculating residual
      res(0) = calculateResEq1(Sw, x_w_L, x_a_L, x_c_L);
      res(1) = calculateResEq2(Sw, x_w_G, x_a_G, x_c_G);
      res(2) = calculateResEq3(Sw, Xa, x_a_L, x_a_G, rho_w, pg, T);
      res(3) = calculateResEq4(Sw, Xc, x_c_L, x_c_G, rho_w, pg, T);
  }

  void ThermalTwoPhaseFlowComponentialMaterialProperties::calculateJacobian(
      double const t, double const dt,
      ParameterLib::SpatialPosition const& x_position,
      MaterialPropertyLib::Property const& pc_model, double const rho_w,
      double const pg, double const Xa, double const Xc,
      double const T,
      JacobianMatrix& Jac, double Sw, double x_w_L, double x_a_L, double x_c_L)
  {
      MaterialPropertyLib::VariableArray vars;
      vars[static_cast<int>(MaterialPropertyLib::Variable::liquid_saturation)] =
          Sw;
      double const pc =
          pc_model.template value<double>(vars, x_position, t, dt);
      // use water density for simplicity
      double const p_vap = calculateVaporPressureNonwet(pc, T, rho_w);
      double const x_w_G = p_vap / pg * x_w_L;
      double const H_a = 8.5906e-6;
      double const H_c = calculateHenryConstant(T, 6.2e-4, 4500);
      double const N_w = rho_w / _water_mol_mass;
      double const N_G = pg / IdealGasConstant / T;
      double const x_a_G = N_w / H_a / pg * x_a_L;
      double const x_c_G = N_w / H_c / pg * x_c_L;
      double const dpc_dSw = pc_model.template dValue<double>(
          vars, MaterialPropertyLib::Variable::liquid_saturation, x_position, t,
          dt);
      double const dpvap_dSw = -p_vap * dpc_dSw / N_w / IdealGasConstant / T;
      Jac.setZero();
      if (Sw <= (1 - x_w_L - x_a_L - x_c_L))
      {
          Jac(0, 0) = 1;
      }
      else
      {
          Jac(0, 1) = -1;
          Jac(0, 2) = -1;
          Jac(0, 3) = -1;
      }
      if ((1 - Sw) <= (1 - x_w_G - x_a_G - x_c_G))
      {
          Jac(1, 0) = -1;
      }
      else
      {
          Jac(1, 0) = -x_w_L * dpvap_dSw / pg;
          Jac(1, 1) = -p_vap / pg;
          Jac(1, 2) = -N_w / pg / H_a;
          Jac(1, 3) = -N_w / pg / H_c;
      }

      Jac(2, 0) = N_w * (Xa - x_a_L) - N_G * (Xa - x_a_G);
      Jac(2, 2) = -(Sw + (1 - Sw) / H_a / IdealGasConstant / T) * N_w;
      Jac(3, 0) = N_w * (Xc - x_c_L) - N_G * (Xc - x_c_G);
      Jac(3, 3) = -(Sw + (1 - Sw) / H_c / IdealGasConstant / T) * N_w;
  }

  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResEq1(
      double Sw, double x_w_L, double x_a_L, double x_c_L) const
  {
      return std::min(Sw, 1 - x_w_L - x_a_L - x_c_L);
  }

  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResEq2(
      double Sw, double x_w_G, double x_a_G, double x_c_G) const
  {
      return std::min(1 - Sw, 1 - x_w_G - x_a_G - x_c_G);
  }

  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResEq3(
      double Sw, double Xa, double x_a_L, double x_a_G, double rho_w, double pg,
      double T) const
  {
      double const N_w = rho_w / _water_mol_mass;
      double const N_G = pg / IdealGasConstant / T;
      return Sw * N_w * (Xa - x_a_L) + (1 - Sw) * N_G * (Xa - x_a_G);
  }

  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResEq4(
      double Sw, double Xc, double x_c_L, double x_c_G, double rho_w, double pg,
      double T) const
  {
      double const N_w = rho_w / _water_mol_mass;
      double const N_G = pg / IdealGasConstant / T;
      return Sw * N_w * (Xc - x_c_L) + (1 - Sw) * N_G * (Xc - x_c_G);
  }

  void ThermalTwoPhaseFlowComponentialMaterialProperties::calculateDerivatives(
      double const t, double const dt,
      ParameterLib::SpatialPosition const& x_position,
      MaterialPropertyLib::Property const& pc_model, double const rho_w,
      double const pg, double const Xa, double const Xc, double const T,
      double Sw, double x_w_L, double x_a_L, double x_c_L, double& dsw_dpg,
      double& dxwG_dpg, double& dxaG_dpg, double& dxcG_dpg, double& dsw_dXa,
      double& dxwG_dXa, double& dxaG_dXa, double& dxcG_dXa, double& dsw_dXc,
      double& dxwG_dXc, double& dxaG_dXc, double& dxcG_dXc, double& dsw_dT,
      double& dxwG_dT, double& dxaG_dT, double& dxcG_dT, double& dxwL_dpg,
      double& dxaL_dpg, double& dxcL_dpg, double& dxwL_dXa, double& dxaL_dXa,
      double& dxcL_dXa, double& dxwL_dXc, double& dxaL_dXc, double& dxcL_dXc,
      double& dxwL_dT, double& dxaL_dT, double& dxcL_dT)
  {
      MaterialPropertyLib::VariableArray vars;
      vars[static_cast<int>(MaterialPropertyLib::Variable::liquid_saturation)] =
          Sw;
      double const pc =
          pc_model.template value<double>(vars, x_position, t, dt);
      // use water density for simplicity
      double const p_vap = calculateVaporPressureNonwet(pc, T, rho_w);
      double const x_w_G = p_vap / pg * x_w_L;
      double const H_a = 8.5906e-6;
      double const H_c = calculateHenryConstant(T, 6.2e-4, 4500);
      double const N_w = rho_w / _water_mol_mass;
      double const N_L = N_w;
      double const N_G = pg / IdealGasConstant / T;
      double const x_a_G = N_w / H_a / pg * x_a_L;
      double const x_c_G = N_w / H_c / pg * x_c_L;
      double const dpvap_dT = calculateDerivativedPgwdT(pc, T, rho_w);
      double const dHa_dT = calculateDerivativedHdT(T, 6.4e-6, 1600);
      double const dHc_dT = calculateDerivativedHdT(T, 6.2e-4, 4500);
      double const dpvap_dpc = calculateDerivativedPgwdPC(pc, T, rho_w);
      double const dpc_dsw = pc_model.template dValue<double>(
          vars, MaterialPropertyLib::Variable::liquid_saturation, x_position, t,
          dt);
      if (Sw <= (1 - x_w_L - x_a_L - x_c_L))
      {
          /*dsw_dpg = 0;
          dxwG_dpg = 0;
          dxaG_dpg = 0;
          dxcG_dpg = 0;*/      // zero by ctor
          dxwG_dXa = -1;
          dxaG_dXa = 1;
          dxwG_dXc = -1;
          dxcG_dXc = 1;
      }
      else if ((1 - Sw) <= (1 - x_w_G - x_a_G - x_c_G))
      {
          // dsw_dpg = 0;
          dxwG_dpg = -x_w_G / pg;
          dxaG_dpg = -x_a_G / pg;
          dxcG_dpg = -x_c_G / pg;
          dxwG_dXa = -p_vap / pg;
          dxaG_dXa = N_L / pg / H_a;
          dxwG_dXc = -p_vap / pg;
          dxcG_dXc = N_L / pg / H_c;
          dxwG_dT = x_w_G * dpvap_dT / p_vap;
          dxaG_dT = -dHa_dT * x_a_G / H_a;
          dxcG_dT = -dHc_dT * x_c_G / H_c;
      }
      else
      {
          // derivatives w.r.t. pg, Xa, Xc, T
          Eigen::Matrix4d Kp, Ka, Kc, Kt;
          Eigen::Vector4d bp, ba, bc, bt;
          Kp.setZero();
          Ka.setZero();
          Kc.setZero();
          Kt.setZero();

          double const coeff_a_1 = (N_L - N_G) * Xa - N_L * x_a_L + N_G * x_a_G;
          double const coeff_a_2 = -(H_a * Sw * pg + (1 - Sw) * N_G);
          double const coeff_a_3 = H_a * pg / N_L;
          double const coeff_c_1 = (N_L - N_G) * Xc - N_L * x_c_L + N_G * x_c_G;
          double const coeff_c_2 = -(H_c * Sw * pg + (1 - Sw) * N_G);
          double const coeff_c_3 = H_c * pg / N_L;
          double const coeff_w_1 = -dpvap_dpc * dpc_dsw * x_w_G / p_vap / p_vap;
          double const pg_over_pvap = pg / p_vap;
          double const N_tot = Sw * N_L + (1 - Sw) * N_G;

          Kp(0, 0) = coeff_a_1;
          Kp(0, 2) = coeff_a_2;
          Kp(1, 0) = coeff_c_1;
          Kp(1, 3) = coeff_c_2;
          Kp(2, 1) = 1;
          Kp(2, 2) = 1;
          Kp(2, 3) = 1;
          Kp(3, 0) = coeff_w_1 * pg;
          Kp(3, 1) = pg_over_pvap;
          Kp(3, 2) = coeff_a_3;
          Kp(3, 3) = coeff_c_3;
          bp << H_a * Sw * x_a_G -
                     (1 - Sw) * (Xa - x_a_G) / IdealGasConstant / T,
              H_c * Sw * x_c_G - (1 - Sw) * (Xc - x_c_G) / IdealGasConstant / T,
              0, -(x_w_G / p_vap + (H_a * x_a_G + H_c * x_c_G) / N_L);

          Ka(0, 0) = coeff_a_1;
          Ka(0, 2) = coeff_a_2;
          Ka(1, 0) = coeff_c_1;
          Ka(1, 3) = coeff_c_2;
          Ka(2, 1) = 1;
          Ka(2, 2) = 1;
          Ka(2, 3) = 1;
          Ka(3, 0) = coeff_w_1;
          Ka(3, 1) = pg_over_pvap;
          Ka(3, 2) = coeff_a_3;
          Ka(3, 3) = coeff_c_3;
          ba << -N_tot, 0, 0, 0;

          Kc(0, 0) = coeff_a_1;
          Kc(0, 2) = coeff_a_2;
          Kc(1, 0) = coeff_c_1;
          Kc(1, 3) = coeff_c_2;
          Kc(2, 1) = 1;
          Kc(2, 2) = 1;
          Kc(2, 3) = 1;
          Kc(3, 0) = coeff_w_1;
          Kc(3, 1) = pg_over_pvap;
          Kc(3, 2) = coeff_a_3;
          Kc(3, 3) = coeff_c_3;
          bc << 0, -N_tot, 0, 0;

          Kt(0, 0) = coeff_a_1;
          Kt(0, 2) = coeff_a_2;
          Kt(1, 0) = coeff_c_1;
          Kt(1, 3) = coeff_c_2;
          Kt(2, 1) = 1;
          Kt(2, 2) = 1;
          Kt(2, 3) = 1;
          Kt(3, 0) = coeff_w_1 * pg;
          Kt(3, 1) = pg_over_pvap;
          Kt(3, 2) = coeff_a_3;
          Kt(3, 3) = coeff_c_3;
          bt << pg * Sw * dHa_dT * x_a_G, pg * Sw * dHc_dT * x_c_G, 0,
              pg * (x_w_G * dpvap_dT / p_vap / p_vap -
                    (dHa_dT * x_a_G + dHc_dT * x_c_G) / N_L);

          Eigen::Vector4d x_dpg = Kp.partialPivLu().solve(bp);
          Eigen::Vector4d x_dXa = Ka.partialPivLu().solve(ba);
          Eigen::Vector4d x_dXc = Kc.partialPivLu().solve(bc);
          Eigen::Vector4d x_dT = Kt.partialPivLu().solve(bt);
          dsw_dpg = x_dpg[0];
          dxwG_dpg = x_dpg[1]; 
          dxaG_dpg = x_dpg[2];
          dxcG_dpg = x_dpg[3];
          dsw_dXa = x_dXa[0];
          dxwG_dXa = x_dXa[1];
          dxaG_dXa = x_dXa[2];
          dxcG_dXa = x_dXa[3];
          dsw_dXc = x_dXc[0];
          dxwG_dXc = x_dXc[1];
          dxaG_dXc = x_dXc[2];
          dxcG_dXc = x_dXc[3];
          dsw_dT = x_dT[0];
          dxwG_dT= x_dT[1];
          dxaG_dT= x_dT[2];
          dxcG_dT= x_dT[3];
      }
      dxwL_dpg = (x_w_G + pg * dxwG_dpg) / p_vap -
                 dpvap_dpc * dpc_dsw * pg * x_w_G * dsw_dpg / p_vap / p_vap;
      dxaL_dpg = (x_a_G + pg * dxaG_dpg) * H_a / N_L;
      dxcL_dpg = (x_c_G + pg * dxcG_dpg) * H_c / N_L;
      dxwL_dXa = pg * dxwG_dXa / p_vap -
                 dpvap_dpc * dpc_dsw * x_w_G * dsw_dXa / p_vap / p_vap;
      dxaL_dXa = pg * H_a * dxaG_dXa / N_L;
      dxcL_dXa = pg * H_c * dxcG_dXa / N_L;
      dxwL_dXc = pg * dxwG_dXc / p_vap -
                 dpvap_dpc * dpc_dsw * x_w_G * dsw_dXc / p_vap / p_vap;
      dxaL_dXc = pg * H_a * dxaG_dXc / N_L;
      dxcL_dXc = pg * H_c * dxcG_dXc / N_L;
      dxwL_dT = pg * dxwG_dT / p_vap -
                pg * x_w_G * (dpvap_dpc * dpc_dsw * dsw_dT + dpvap_dT) / p_vap /
                    p_vap;
      dxaL_dT = pg * (dHa_dT * x_a_G + H_a * dxaG_dT) / N_L;
      dxcL_dT = pg * (dHc_dT * x_c_G + H_c * dxcG_dT) / N_L;
  }
      
}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
