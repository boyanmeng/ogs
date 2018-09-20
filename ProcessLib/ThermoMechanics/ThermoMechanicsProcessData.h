/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <utility>

#include <Eigen/Eigen>

#include "ProcessLib/Parameter/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}
namespace ProcessLib
{
namespace ThermoMechanics
{
template <int DisplacementDim>
struct ThermoMechanicsProcessData
{
    ThermoMechanicsProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Parameter<double> const& reference_solid_density_,
        Parameter<double> const& linear_thermal_expansion_coefficient_,
        Parameter<double> const& specific_heat_capacity_,
        Parameter<double> const& thermal_conductivity_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_)
        : material{std::move(material_)},
          reference_solid_density(reference_solid_density_),
          linear_thermal_expansion_coefficient(
              linear_thermal_expansion_coefficient_),
          specific_heat_capacity(specific_heat_capacity_),
          thermal_conductivity(thermal_conductivity_),
          specific_body_force(specific_body_force_)
    {
    }

    ThermoMechanicsProcessData(ThermoMechanicsProcessData&& other)
        : material{std::move(other.material)},
          reference_solid_density(other.reference_solid_density),
          linear_thermal_expansion_coefficient(
              other.linear_thermal_expansion_coefficient),
          specific_heat_capacity(other.specific_heat_capacity),
          thermal_conductivity(other.thermal_conductivity),
          specific_body_force(other.specific_body_force),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    ThermoMechanicsProcessData(ThermoMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermoMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(ThermoMechanicsProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    Parameter<double> const& reference_solid_density;
    Parameter<double> const& linear_thermal_expansion_coefficient;
    Parameter<double> const& specific_heat_capacity;
    Parameter<double> const&
        thermal_conductivity;  // TODO To be changed as a matrix type variable.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt = 0;
    double t = 0;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoMechanics
}  // namespace ProcessLib
