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
      return 1;
  }
  double
  ThermalTwoPhaseFlowComponentialMaterialProperties::calculateVaporPressureNonwet(
      const double pc, const double T, const double mass_density_water) const
  {
      return 1;
  }
  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateDerivativedPsatdT(
      const double T) const
  {
      return 0;
  }
  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateDerivativedPgwdT(
      const double pc, const double T, const double mass_density_water) const
  {
      return 0;
  }
  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateDerivativedPgwdPC(
      const double pc, const double T, const double mass_density_water) const
  {
      return 0;
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
          MaterialPropertyLib::Property const& pc_model, double const err_tol, double const rho_w,
          double const H_a,
          double const pl,
                                  double const Xa,
                                  double const T, double& Sw, double& x_w_L, double& x_a_L, double& dsw_dpl, double& dxwG_dpl, double& dxaG_dpl,
          double& dsw_dXa, double& dxwG_dXa, double& dxaG_dXa, double& dxwL_dpl,
          double& dxaL_dpl, double& dxwL_dXa, double& dxaL_dXa)
  {
        // Local Newton solver
          using LocalJacobianMatrix =
              Eigen::Matrix<double, 3, 3, Eigen::RowMajor>;
          using LocalResidualVector = Eigen::Matrix<double, 3, 1>;
          using LocalUnknownVector = Eigen::Matrix<double, 3, 1>;
          LocalJacobianMatrix J_loc;

          Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(3);
          auto const update_residual = [&](LocalResidualVector& residual) {
              calculateResidual(t, dt, x_position, pc_model, rho_w, H_a, pl, Xa,
                                T, Sw,
                                x_w_L, x_a_L, residual);
          };

          auto const update_jacobian = [&](LocalJacobianMatrix& jacobian) {
              calculateJacobian(t, dt, x_position, pc_model, rho_w, H_a, pl, Xa,  T, jacobian, Sw, x_w_L, x_a_L
                                );  // for solution dependent Jacobians
          };

          auto const update_solution =
              [&](LocalUnknownVector const& increment) {
                  // increment solution vectors
                  Sw += increment[0];
                  x_w_L += increment[1];
                  x_a_L += increment[2];
              };

          // TODO Make the following choice of maximum iterations and
          // convergence criteria available from the input file configuration.
          // See Ehlers material model implementation for the example.
          const int maximum_iterations(20);
          const double tolerance(err_tol);

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
      
      calculateDerivatives(t, dt, x_position, pc_model, rho_w, H_a, pl, Xa, T, J_loc,
                           Sw, x_w_L, x_a_L, dsw_dpl, dxwG_dpl, dxaG_dpl, dsw_dXa,
          dxwG_dXa, dxaG_dXa, dxwL_dpl, dxaL_dpl,
          dxwL_dXa, dxaL_dXa);
      return true;
  }

  void ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResidual(
      double const t, double const dt,
      ParameterLib::SpatialPosition const& x_position,
      MaterialPropertyLib::Property const& pc_model, double const rho_w,
      double const H_a,
      double const pl, double const Xa, double const T,
      double Sw, double x_w_L, double x_a_L, ResidualVector& res)
  {
      MaterialPropertyLib::VariableArray vars;
      vars[static_cast<int>(MaterialPropertyLib::Variable::liquid_saturation)] =
          Sw;
      double const pc = pc_model.template value<double>(vars, x_position, t, dt);
      // use water density for simplicity
      double const pg = pl + pc;
      double const p_vap =
          calculateVaporPressureNonwet(pc, T, rho_w);
      double const x_w_G = p_vap / pg * x_w_L;
      double const N_w = rho_w / _water_mol_mass;
      double const x_a_G = N_w / H_a / pg * x_a_L;
      // calculating residual
      res(0) = calculateResEq1(Sw, x_w_L, x_a_L);
      res(1) = calculateResEq2(Sw, x_w_G, x_a_G);
      res(2) = calculateResEq3(Sw, Xa, x_a_L, x_a_G, rho_w, pg, T);
  }

  void ThermalTwoPhaseFlowComponentialMaterialProperties::calculateJacobian(
      double const t, double const dt,
      ParameterLib::SpatialPosition const& x_position,
      MaterialPropertyLib::Property const& pc_model, double const rho_w,
      double const H_a,
      double const pl, double const Xa,
      double const T,
      JacobianMatrix& Jac, double Sw, double x_w_L, double x_a_L)
  {
      MaterialPropertyLib::VariableArray vars;
      vars[static_cast<int>(MaterialPropertyLib::Variable::liquid_saturation)] =
          std::clamp(Sw, 0., 1.);
      double const pc =
          pc_model.template value<double>(vars, x_position, t, dt);
      double const pg = pl + pc;
      // use water density for simplicity
      double const p_vap = calculateVaporPressureNonwet(pc, T, rho_w);
      double const x_w_G = p_vap / pg * x_w_L;
      double const N_w = rho_w / _water_mol_mass;
      double const N_G = pg / IdealGasConstant / T;
      double const x_a_G = N_w / H_a / pg * x_a_L;
      double const dpc_dSw = pc_model.template dValue<double>(
          vars, MaterialPropertyLib::Variable::liquid_saturation, x_position, t,
          dt);
      double const dpvap_dpc = calculateDerivativedPgwdPC(pc, T, rho_w);
      double const dpvap_dSw = dpvap_dpc * dpc_dSw;
      double const dxaG_dpc = -N_w * x_a_L / H_a / pg / pg;
      Jac.setZero();
      if (Sw <= (1 - x_w_L - x_a_L))
      {
          Jac(0, 0) = 1;
      }
      else
      {
          Jac(0, 1) = -1;
          Jac(0, 2) = -1;
      }
      if ((1 - Sw) <= (1 - x_w_G - x_a_G))
      {
          Jac(1, 0) = -1;
      }
      else
      {
          Jac(1, 0) = -x_w_L * (dpvap_dSw / pg - p_vap * dpc_dSw / pg / pg) -
                      dxaG_dpc * dpc_dSw;
          Jac(1, 1) = -p_vap / pg;
          Jac(1, 2) = -N_w / pg / H_a;
      }

      Jac(2, 0) =
          N_w * (Xa - x_a_L) - N_G * (Xa - x_a_G) +
          (1 - Sw) * ((Xa - x_a_G) / IdealGasConstant / T - N_G * dxaG_dpc) *
              dpc_dSw;           
      Jac(2, 2) = -(Sw + (1 - Sw) / H_a / IdealGasConstant / T) * N_w;
  }

  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResEq1(
      double Sw, double x_w_L, double x_a_L) const
  {
      return std::min(Sw, 1 - x_w_L - x_a_L);
  }

  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResEq2(
      double Sw, double x_w_G, double x_a_G) const
  {
      return std::min(1 - Sw, 1 - x_w_G - x_a_G);
  }

  double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResEq3(
      double Sw, double Xa, double x_a_L, double x_a_G, double rho_w, double pg,
      double T) const
  {
      double const N_w = rho_w / _water_mol_mass;
      double const N_G = pg / IdealGasConstant / T;
      return Sw * N_w * (Xa - x_a_L) + (1 - Sw) * N_G * (Xa - x_a_G);
  }

  void ThermalTwoPhaseFlowComponentialMaterialProperties::calculateDerivatives(
      double const t, double const dt,
      ParameterLib::SpatialPosition const& x_position,
      MaterialPropertyLib::Property const& pc_model, double const rho_w,
      double const H_a, double const pl, double const Xa, double const T,
      JacobianMatrix Jac,
      double Sw, double x_w_L, double x_a_L, double& dsw_dpl, double& dxwG_dpl,
      double& dxaG_dpl, double& dsw_dXa, double& dxwG_dXa, double& dxaG_dXa,
      double& dxwL_dpl, double& dxaL_dpl, double& dxwL_dXa, double& dxaL_dXa)
  {
      MaterialPropertyLib::VariableArray vars;
      vars[static_cast<int>(MaterialPropertyLib::Variable::liquid_saturation)] =
          std::clamp(Sw, 0., 1.);
      double const pc =
          pc_model.template value<double>(vars, x_position, t, dt);
      double const pg = pl + pc;
      // use water density for simplicity
      double const p_vap = calculateVaporPressureNonwet(pc, T, rho_w);
      double const x_w_G = p_vap / pg * x_w_L;
      double const N_w = rho_w / _water_mol_mass;
      double const N_L = N_w;
      double const N_G = pg / IdealGasConstant / T;
      double const x_a_G = N_w / H_a / pg * x_a_L;
      double const dpvap_dT = calculateDerivativedPgwdT(pc, T, rho_w);
      double const dpvap_dpc = calculateDerivativedPgwdPC(pc, T, rho_w);
      double const dpc_dSw = pc_model.template dValue<double>(
          vars, MaterialPropertyLib::Variable::liquid_saturation, x_position, t,
          dt);
      double const dpvap_dSw = dpvap_dpc * dpc_dSw;

      Eigen::MatrixXd dF_dPV(3, 2), dSV_dPV(3, 2);
      dF_dPV.setZero();

      if ((1 - Sw) > (1 - x_w_G - x_a_G))
      {
          dF_dPV(1, 0) = (p_vap * x_w_L + N_L * x_a_L / H_a) / pg / pg;
      }
      dF_dPV(2, 0) = (1 - Sw) * (Xa - x_a_G + N_L * x_a_L / H_a / pg) /
                     IdealGasConstant / T;
      dF_dPV(2, 1) = Sw * N_L + (1 - Sw) * N_G;

      dSV_dPV = -Jac.partialPivLu().solve(dF_dPV);
      dsw_dpl = dSV_dPV(0, 0);
      dsw_dXa = dSV_dPV(0, 1);
      dxwL_dpl = dSV_dPV(1, 0);
      dxwL_dXa = dSV_dPV(1, 1);
      dxaL_dpl = dSV_dPV(2, 0);
      dxaL_dXa = dSV_dPV(2, 1);

      dxwG_dpl = (dpvap_dSw * dsw_dpl * x_w_L + p_vap * dxwL_dpl) / pg -
                 p_vap * x_w_L * (1 + dpc_dSw * dsw_dpl) / pg / pg;
      dxwG_dXa = (dpvap_dSw * dsw_dXa * x_w_L + p_vap * dxwL_dXa) / pg -
                 p_vap * x_w_L * dpc_dSw * dsw_dXa / pg / pg;
      dxaG_dpl = N_L / H_a *
                 (dxaL_dpl / pg - x_a_L * (1 + dpc_dSw * dsw_dpl) / pg / pg);
      dxaG_dXa =
          N_L / H_a * (dxaL_dXa / pg - x_a_L * dpc_dSw * dsw_dXa / pg / pg);
      }
      
}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
