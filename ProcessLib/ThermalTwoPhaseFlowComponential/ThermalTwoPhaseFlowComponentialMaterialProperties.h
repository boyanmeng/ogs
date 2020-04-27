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

#include <memory>
#include <vector>
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/Fluid/WaterVaporProperties/WaterVaporProperties.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib::ThermalTwoPhaseFlowComponential
{
class ThermalTwoPhaseFlowComponentialMaterialProperties
{
public:
    ThermalTwoPhaseFlowComponentialMaterialProperties(
        std::unique_ptr<MaterialLib::Fluid::WaterVaporProperties>&&
            water_vapor_properties);

    /// water vapor saturation pressure
    double calculateSaturatedVaporPressure(const double T) const;
    /// partial water vapor pressure in nonwetting phase
    /// Kelvin equation
    double calculateVaporPressureNonwet(const double pc, const double T,
                                        const double mass_density_water) const;
    /// Derivative of SaturatedVaporPressure in terms of T
    double calculateDerivativedPsatdT(const double T) const;
    /// Derivative of partial vapor pressure in terms of T
    double calculateDerivativedPgwdT(const double pc, const double T,
                                     const double mass_density_water) const;
    /// Derivative of partial vapor pressure in terms of PC
    double calculateDerivativedPgwdPC(const double pc, const double T,
                                      const double mass_density_water) const;
    ///
    double calculatedDensityNonwetdT(const double p_air_nonwet,
                                     const double p_vapor_nonwet,
                                     const double pc, const double T,
                                     const double mass_density_water) const;
    /// Specific enthalpy of water vapor
    double getWaterVaporEnthalpySimple(
        const double temperature,
        const double heat_capacity_water_vapor,
        const double pg,
        const double latent_heat_evaporation) const;
    /// Specific enthalpy of air
    double getAirEnthalpySimple(const double temperature,
                                const double heat_capacity_dry_air,
                                const double /*pg*/) const;
    /// Specific enthalpy of liquid water
    double getLiquidWaterEnthalpySimple(const double temperature,
                                        const double heat_capacity_liquid_water,
                                        const double /*pl*/) const;
    /*
    bool computeConstitutiveRelation(
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        const int material_id,
        double const pg,
        double const X,
        double const T,
        double& Sw,
        double& X_m,
        double& dsw_dpg,
        double& dsw_dX,
        double& dxm_dpg,
        double& dxm_dX);
    */
private:
    double const& _air_mol_mass = MaterialLib::PhysicalConstant::MolarMass::Air;
    std::unique_ptr<MaterialLib::Fluid::WaterVaporProperties> const
        _water_vapor_properties;
    /*
    static int const jacobian_residual_size = 2;
    using ResidualVector = Eigen::Matrix<double, jacobian_residual_size, 1>;
    using JacobianMatrix =
        Eigen::Matrix<double, jacobian_residual_size, jacobian_residual_size,
                      Eigen::RowMajor>;
    using UnknownVector = Eigen::Matrix<double, jacobian_residual_size, 1>;

private:
    
    // Calculates the residual vector.
  
    void calculateResidual(const int material_id, double const pl,
                           double const X, double const T, double Sw,
                           double rho_h2_wet, ResidualVector& res);
    
    // Calculates the Jacobian.
    
    void calculateJacobian(const int material_id, double const t,
                           ParameterLib::SpatialPosition const& x,
                           double const pl, double const X, double const T,
                           JacobianMatrix& Jac, double Sw, double rho_h2_wet);
    // Complementary condition 1
    // for calculating molar fraction of light component in the liquid phase
    
    double calculateEquilibiumRhoWetLight(double const pg, double const Sw,
                                          double const rho_wet_h2) const;
    // Complementary condition 2
    // for calculating the saturation
    
    double calculateSaturation(double , double X, double Sw,
                               double rho_wet_h2, double rho_nonwet_h2,
                               double ) const;
    
    // Calculate the derivatives using the analytical way
    
    double calculatedSwdP(double pl, double S, double rho_wet_h2,
    
    double calculatedSwdX(double const pl, const double , const double S,
                          const double rho_wet_h2, double const T,
                          int current_material_id) const;
    double calculatedXmdX(double pl, double Sw, double rho_wet_h2, double dSwdX,
                          int current_material_id) const;
  
    double calculatedXmdP(double pl, double Sw, double rho_wet_h2, double dSwdP,
                          int current_material_id) const;
    */
};

}  // namespace ProcessLib::ThermalTwoPhaseFlowComponential
