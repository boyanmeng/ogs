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
    return - henry * delta / T / T;
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
      
}  // namespace ThermalTwoPhaseFlowComponential
}  // namespace ProcessLib
