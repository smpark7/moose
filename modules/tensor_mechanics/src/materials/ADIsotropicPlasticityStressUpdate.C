//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADIsotropicPlasticityStressUpdate.h"

#include "Function.h"
#include "ElasticityTensorTools.h"

registerMooseObject("TensorMechanicsApp", ADIsotropicPlasticityStressUpdate);

InputParameters
ADIsotropicPlasticityStressUpdate::validParams()
{
  InputParameters params = ADRadialReturnStressUpdate::validParams();
  params.addClassDescription("This class uses the discrete material in a radial return isotropic "
                             "plasticity model.  This class is one of the basic radial return "
                             "constitutive models, yet it can be used in conjunction with other "
                             "creep and plasticity materials for more complex simulations.");
  // Linear strain hardening parameters
  params.addParam<FunctionName>("yield_stress_function",
                                "Yield stress as a function of temperature");
  params.addParam<Real>(
      "yield_stress", 0.0, "The point at which plastic strain begins accumulating");
  params.addParam<FunctionName>("hardening_function",
                                "True stress as a function of plastic strain");
  params.addParam<Real>("hardening_constant", 0.0, "Hardening slope");
  params.addCoupledVar("temperature", 0.0, "Coupled Temperature");
  params.addDeprecatedParam<std::string>(
      "plastic_prepend",
      "",
      "String that is prepended to the plastic_strain Material Property",
      "This has been replaced by the 'base_name' parameter");
  params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";
  return params;
}

ADIsotropicPlasticityStressUpdate::ADIsotropicPlasticityStressUpdate(
    const InputParameters & parameters)
  : ADRadialReturnStressUpdate(parameters),
    _plastic_prepend(getParam<std::string>("plastic_prepend")),
    _yield_stress_function(
        isParamValid("yield_stress_function") ? &getFunction("yield_stress_function") : NULL),
    _yield_stress(getParam<Real>("yield_stress")),
    _hardening_constant(getParam<Real>("hardening_constant")),
    _hardening_function(isParamValid("hardening_function") ? &getFunction("hardening_function")
                                                           : NULL),
    _yield_condition(-1.0), // set to a non-physical value to catch uninitalized yield condition
    _hardening_slope(0.0),
    _plastic_strain(
        declareADProperty<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_strain")),
    _plastic_strain_old(
        getMaterialPropertyOld<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_strain")),
    _hardening_variable(declareADProperty<Real>(_base_name + "hardening_variable")),
    _hardening_variable_old(getMaterialPropertyOld<Real>(_base_name + "hardening_variable")),
    _temperature(adCoupledValue("temperature"))
{
  if (parameters.isParamSetByUser("yield_stress") && _yield_stress <= 0.0)
    mooseError("Yield stress must be greater than zero");

  if (_yield_stress_function == NULL && !parameters.isParamSetByUser("yield_stress"))
    mooseError("Either yield_stress or yield_stress_function must be given");

  if (!parameters.isParamSetByUser("hardening_constant") && !isParamValid("hardening_function"))
    mooseError("Either hardening_constant or hardening_function must be defined");

  if (parameters.isParamSetByUser("hardening_constant") && isParamValid("hardening_function"))
    mooseError(
        "Only the hardening_constant or only the hardening_function can be defined but not both");
}

void
ADIsotropicPlasticityStressUpdate::initQpStatefulProperties()
{
  _hardening_variable[_qp] = 0.0;
  _plastic_strain[_qp].zero();
}

void
ADIsotropicPlasticityStressUpdate::propagateQpStatefulProperties()
{
  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];

  propagateQpStatefulPropertiesRadialReturn();
}

void
ADIsotropicPlasticityStressUpdate::computeStressInitialize(
    const ADReal & effective_trial_stress, const ADRankFourTensor & elasticity_tensor)
{
  computeYieldStress(elasticity_tensor);

  _yield_condition = effective_trial_stress - _hardening_variable_old[_qp] - _yield_stress;
  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
}

ADReal
ADIsotropicPlasticityStressUpdate::computeResidual(const ADReal & effective_trial_stress,
                                                   const ADReal & scalar)
{
  ADReal residual = 0.0;

  mooseAssert(_yield_condition != -1.0,
              "the yield stress was not updated by computeStressInitialize");

  if (_yield_condition > 0.0)
  {
    _hardening_slope = computeHardeningDerivative(scalar);
    _hardening_variable[_qp] = computeHardeningValue(scalar);

    residual =
        (effective_trial_stress - _hardening_variable[_qp] - _yield_stress) / _three_shear_modulus -
        scalar;
  }
  return residual;
}

ADReal
ADIsotropicPlasticityStressUpdate::computeDerivative(const ADReal & /*effective_trial_stress*/,
                                                     const ADReal & /*scalar*/)
{
  if (_yield_condition > 0.0)
    return -1.0 - _hardening_slope / _three_shear_modulus;

  return 1.0;
}

void
ADIsotropicPlasticityStressUpdate::iterationFinalize(ADReal scalar)
{
  if (_yield_condition > 0.0)
    _hardening_variable[_qp] = computeHardeningValue(scalar);
}

void
ADIsotropicPlasticityStressUpdate::computeStressFinalize(
    const ADRankTwoTensor & plastic_strain_increment)
{
  _plastic_strain[_qp] += plastic_strain_increment;
}

ADReal
ADIsotropicPlasticityStressUpdate::computeHardeningValue(const ADReal & scalar)
{
  if (_hardening_function)
  {
    const Real strain_old = _effective_inelastic_strain_old[_qp];
    const Point p;
    const Real t = strain_old + MetaPhysicL::raw_value(scalar);

    DualReal hardening_function_value = _hardening_function->value(t, p);
    hardening_function_value.derivatives() =
        (strain_old + scalar).derivatives() * _hardening_function->timeDerivative(t, p);

    return hardening_function_value - _yield_stress;
  }

  return _hardening_variable_old[_qp] + _hardening_slope * scalar;
}

ADReal
ADIsotropicPlasticityStressUpdate::computeHardeningDerivative(const ADReal & /*scalar*/)
{
  if (_hardening_function)
  {
    const Real strain_old = _effective_inelastic_strain_old[_qp];
    const Point p; // Always (0,0,0)

    return _hardening_function->timeDerivative(strain_old, p);
  }

  return _hardening_constant;
}

void
ADIsotropicPlasticityStressUpdate::computeYieldStress(
    const ADRankFourTensor & /*elasticity_tensor*/)
{
  if (_yield_stress_function)
  {
    const Point p;

    _yield_stress = _yield_stress_function->value(MetaPhysicL::raw_value(_temperature[_qp]), p);
    _yield_stress.derivatives() =
        _temperature[_qp].derivatives() *
        _yield_stress_function->timeDerivative(MetaPhysicL::raw_value(_temperature[_qp]), p);

    if (_yield_stress <= 0.0)
      mooseException("In ",
                     _name,
                     ": The calculated yield stress (",
                     _yield_stress.value(),
                     ") is less than zero");
  }
}
