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
#include "BaseLib/Logging.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "NumLib/NewtonRaphson.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/SpatialPosition.h"

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

//bool ThermalTwoPhaseFlowComponentialMaterialProperties::computeConstitutiveRelation(
//    double const t,
//    ParameterLib::SpatialPosition const& x,
//    int const material_id,
//    double const pg,
//    double const X,
//    double const T,
//    double& Sw,
//    double& X_m,
//    double& dsw_dpg,
//    double& dsw_dX,
//    double& dxm_dpg,
//    double& dxm_dX)
//{
//    {  // Local Newton solver
//        using LocalJacobianMatrix =
//            Eigen::Matrix<double, 2, 2, Eigen::RowMajor>;
//        using LocalResidualVector = Eigen::Matrix<double, 2, 1>;
//        using LocalUnknownVector = Eigen::Matrix<double, 2, 1>;
//        LocalJacobianMatrix J_loc;
//
//        Eigen::PartialPivLU<LocalJacobianMatrix> linear_solver(2);
//        auto const update_residual = [&](LocalResidualVector& residual) {
//            calculateResidual(material_id, pg, X, T, Sw, X_m, residual);
//        };
//
//        auto const update_jacobian = [&](LocalJacobianMatrix& jacobian) {
//            calculateJacobian(material_id, t, x, pg, X, T, jacobian, Sw,
//                              X_m);  // for solution dependent Jacobians
//        };
//
//        auto const update_solution = [&](LocalUnknownVector const& increment) {
//            // increment solution vectors
//            Sw += increment[0];
//            X_m += increment[1];
//        };
//
//        // TODO Make the following choice of maximum iterations and convergence
//        // criteria available from the input file configuration. See Ehlers
//        // material model implementation for the example.
//        const int maximum_iterations(20);
//        const double tolerance(1.e-14);
//
//        auto newton_solver = NumLib::NewtonRaphson<
//            decltype(linear_solver), LocalJacobianMatrix,
//            decltype(update_jacobian), LocalResidualVector,
//            decltype(update_residual), decltype(update_solution)>(
//            linear_solver, update_jacobian, update_residual, update_solution,
//            {maximum_iterations, tolerance});
//
//        auto const success_iterations = newton_solver.solve(J_loc);
//
//        if (!success_iterations)
//        {
//            return false;
//        }
//    }
//    dsw_dpg = calculatedSwdP(pg, Sw, X_m, T, material_id);
//    dsw_dX = calculatedSwdX(pg, X, Sw, X_m, T, material_id);
//    dxm_dpg = calculatedXmdP(pg, Sw, X_m, dsw_dpg, material_id);
//    dxm_dX = calculatedXmdX(pg, Sw, X_m, dsw_dX, material_id);
//    return true;
//}
//void ThermalTwoPhaseFlowComponentialMaterialProperties::calculateResidual(
//    const int material_id, double const pl, double const X, double const T,
//    double Sw, double rho_h2_wet, ResidualVector& res)
//{
//    const double pg =
//        pl + _capillary_pressure_models[material_id]->getCapillaryPressure(Sw);
//    const double rho_h2_nonwet = pg * H2 / IdealGasConstant / T;
//
//    // calculating residual
//    res(0) = calculateEquilibiumRhoWetLight(pg, Sw, rho_h2_wet);
//    res(1) = calculateSaturation(pl, X, Sw, rho_h2_wet, rho_h2_nonwet, T);
//}
//
//void ThermalTwoPhaseFlowComponentialMaterialProperties::calculateJacobian(
//    const int material_id, double const /*t*/,
//    ParameterLib::SpatialPosition const& /*x*/, double const pl,
//    double const /*X*/, double const T, JacobianMatrix& Jac, double Sw,
//    double rho_h2_wet)
//{
//    const double pg =
//        pl + _capillary_pressure_models[material_id]->getCapillaryPressure(Sw);
//    const double rho_h2_nonwet = pg * H2 / IdealGasConstant / T;
//    double const rho_equili_h2_wet = pg * HenryConstantH2 * H2;
//    double const dPC_dSw =
//        _capillary_pressure_models[material_id]->getdPcdS(Sw);
//    double const drhoh2wet_dpg = HenryConstantH2 * H2;
//    Jac.setZero();
//    if ((1 - Sw) < (rho_equili_h2_wet - rho_h2_wet))
//    {
//        Jac(0, 0) = -1;
//    }
//    else
//    {
//        Jac(0, 0) = drhoh2wet_dpg * dPC_dSw;
//        Jac(0, 1) = -1;
//    }
//
//    Jac(1, 0) = rho_h2_nonwet - rho_h2_wet;
//    Jac(1, 1) = -Sw;
//}
//
///** Complementary condition 1
//* for calculating molar fraction of light component in the liquid phase
//*/
//double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateEquilibiumRhoWetLight(
//    double const pg, double const Sw, double const rho_wet_h2) const
//{
//    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
//    return std::min(1 - Sw, rho_equilibrium_wet_h2 - rho_wet_h2);
//}
//
///** Complementary condition 2
//* for calculating the saturation
//*/
//double ThermalTwoPhaseFlowComponentialMaterialProperties::calculateSaturation(
//    double /*PL*/, double X, double Sw, double rho_wet_h2, double rho_nonwet_h2,
//    double /*T*/) const
//{
//    return X - (Sw * rho_wet_h2 + (1 - Sw) * rho_nonwet_h2);
//}
//
///**
//* Calculate the derivatives using the analytical way
//*/
//double ThermalTwoPhaseFlowComponentialMaterialProperties::calculatedSwdP(
//    double pl, double S, double rho_wet_h2, double const T,
//    int current_material_id) const
//{
//    const double pg =
//        pl +
//        _capillary_pressure_models[current_material_id]->getCapillaryPressure(
//            S);
//    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
//    if ((1 - S) < (rho_equilibrium_wet_h2 - rho_wet_h2))
//    {
//        return 0.0;
//    }
//    double const drhoh2wet_dpg = HenryConstantH2 * H2;
//    double const drhoh2nonwet_dpg = H2 / IdealGasConstant / T;
//    double const alpha =
//        ((drhoh2nonwet_dpg - drhoh2wet_dpg) * (1 - S) + drhoh2wet_dpg);
//    double const beta = (drhoh2nonwet_dpg - drhoh2wet_dpg) *
//                        pg;  // NOTE here should be PG^h, but we ignore vapor
//    double const dPC_dSw =
//        _capillary_pressure_models[current_material_id]->getdPcdS(S);
//    return alpha / (beta - alpha * dPC_dSw);
//}
///**
//* Calculate the derivatives using the analytical way
//*/
//double ThermalTwoPhaseFlowComponentialMaterialProperties::calculatedSwdX(
//    double const pl, const double /*X*/, const double S,
//    const double rho_wet_h2, double const T, int current_material_id) const
//{
//    const double pg =
//        pl +
//        _capillary_pressure_models[current_material_id]->getCapillaryPressure(
//            S);
//    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
//    if ((1 - S) < (rho_equilibrium_wet_h2 - rho_wet_h2))
//    {
//        return 0.0;
//    }
//    double const drhoh2wet_dpg = HenryConstantH2 * H2;
//    double const drhoh2nonwet_dpg = H2 / IdealGasConstant / T;
//    double const alpha =
//        ((drhoh2nonwet_dpg - drhoh2wet_dpg) * (1 - S) + drhoh2wet_dpg);
//    double const beta = (drhoh2nonwet_dpg - drhoh2wet_dpg) *
//                        pg;  // NOTE here should be PG^h, but we ignore vapor
//    double const dPC_dSw =
//        _capillary_pressure_models[current_material_id]->getdPcdS(S);
//    return -1 / (beta - alpha * dPC_dSw);
//}
///**
//* Calculate the derivatives using the analytical way
//*/
//double ThermalTwoPhaseFlowComponentialMaterialProperties::calculatedXmdX(
//    double pl, double Sw, double rho_wet_h2, double dSwdX,
//    int current_material_id) const
//{
//    const double pg =
//        pl +
//        _capillary_pressure_models[current_material_id]->getCapillaryPressure(
//            Sw);
//    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
//    double const dPC_dSw =
//        _capillary_pressure_models[current_material_id]->getdPcdS(Sw);
//    if ((1 - Sw) < (rho_equilibrium_wet_h2 - rho_wet_h2))
//    {
//        return 1.0;
//    }
//    return HenryConstantH2 * H2 * dPC_dSw * dSwdX;
//}
///**
//* Calculate the derivatives using the analytical way
//*/
//double ThermalTwoPhaseFlowComponentialMaterialProperties::calculatedXmdP(
//    double pl, double Sw, double rho_wet_h2, double dSwdP,
//    int current_material_id) const
//{
//    const double pg =
//        pl +
//        _capillary_pressure_models[current_material_id]->getCapillaryPressure(
//            Sw);
//    double const rho_equilibrium_wet_h2 = pg * HenryConstantH2 * H2;
//    double const dPC_dSw =
//        _capillary_pressure_models[current_material_id]->getdPcdS(Sw);
//    if ((1 - Sw) < (rho_equilibrium_wet_h2 - rho_wet_h2))
//    {
//        return 0.0;
//    }
//    return HenryConstantH2 * H2 * (1 + dPC_dSw * dSwdP);
//}
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
  
  double
  ThermalTwoPhaseFlowComponentialMaterialProperties::getLiquidWaterEnthalpySimple(
      const double temperature,
      const double heat_capacity_liquid_water,
      const double /*pl*/) const
  {
      return heat_capacity_liquid_water * (temperature - CelsiusZeroInKelvin);
  }
      
}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
