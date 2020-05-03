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
             IdealGasConstant * (temperature - CelsiusZeroInKelvin) /
                 _air_mol_mass;
  }

  double ThermalTwoPhaseFlowComponentialMaterialProperties::
      getContaminantEnthalpySimple(const double temperature,
                                   const double contaminant_heat_capacity,
                                   const double contaminant_molar_mass,
                                   const double) const
  {
      return contaminant_heat_capacity * (temperature - CelsiusZeroInKelvin) +
             IdealGasConstant * (temperature - CelsiusZeroInKelvin) /
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
                                  double const Xc, double const T, double& Sw,
                                  double& x_w_L, double& x_a_L, double& x_c_L)
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
      // dsw_dpg = calculatedSwdP(pg, Sw, X_m, T, material_id);
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
      double const H_a = calculateHenryConstant(T, 6.4e-6, 1600);
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
      double const H_a = calculateHenryConstant(T, 6.4e-6, 1600);
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
      
}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
