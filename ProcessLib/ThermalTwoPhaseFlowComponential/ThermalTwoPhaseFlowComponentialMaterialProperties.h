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
#include "ParameterLib/SpatialPosition.h"
#include "MaterialLib/MPL/Property.h"

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

    /// calculate Henry constants
    double calculateHenryConstant(const double T, const double H_ref, const double delta) const;
    /// calculate Henry constant derivative
    double calculateDerivativedHdT(const double T, const double H_ref,
                                  const double delta) const;
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
    /// Specific enthalpy of contaminant
    double getContaminantEnthalpySimple(const double temperature,
                                const double contaminant_heat_capacity,
                                const double contaminant_molar_mass,
                                const double /*pg*/) const;
    /// Specific enthalpy of liquid water
    double getLiquidWaterEnthalpySimple(const double temperature,
                                        const double heat_capacity_liquid_water,
                                        const double /*pl*/) const;

private:
    double const& _air_mol_mass = MaterialLib::PhysicalConstant::MolarMass::Air;
    double const& _water_mol_mass = MaterialLib::PhysicalConstant::MolarMass::Water;
    std::unique_ptr<MaterialLib::Fluid::WaterVaporProperties> const
        _water_vapor_properties;
};

}  // namespace ProcessLib::ThermalTwoPhaseFlowComponential
