/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSMomentumTurbulentViscosityLaplaceForm.h"

template <>
InputParameters
validParams<INSMomentumTurbulentViscosityLaplaceForm>()
{
  InputParameters params = validParams<Kernel>();
  params += validParams<ScalarTransportBase>();

  params.addClassDescription("This class provides the methods for the turbulent viscosity term in "
                             "the INS momentum component equations.");
  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("kin", "The turbulent kinetic energy.");
  params.addRequiredCoupledVar("epsilon", "The turbulent dissipation.");

  // Required parameters
  params.addRequiredParam<unsigned>(
      "component",
      "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");

  // Optional parameters
  params.addParam<MaterialPropertyName>("mu_name", "mu", "The name of the dynamic viscosity");
  params.addParam<MaterialPropertyName>("rho_name", "rho", "The name of the density");

  return params;
}

INSMomentumTurbulentViscosityLaplaceForm::INSMomentumTurbulentViscosityLaplaceForm(
    const InputParameters & parameters)
  : Kernel(parameters),
    ScalarTransportBase(parameters),

    // Coupled variables
    _u_vel(coupledValue("u")),
    _v_vel(coupledValue("v")),
    _w_vel(coupledValue("w")),
    _kin(coupledValue("kin")),
    _epsilon(coupledValue("epsilon")),

    // Gradients
    _grad_u_vel(coupledGradient("u")),
    _grad_v_vel(coupledGradient("v")),
    _grad_w_vel(coupledGradient("w")),

    // Variable numberings
    _u_vel_var_number(coupled("u")),
    _v_vel_var_number(coupled("v")),
    _w_vel_var_number(coupled("w")),
    _kin_var_number(coupled("kin")),
    _epsilon_var_number(coupled("epsilon")),

    // Required parameters
    _component(getParam<unsigned>("component")),

    // Material properties
    _mu(getMaterialProperty<Real>("mu_name")),
    _rho(getMaterialProperty<Real>("rho_name"))
{
}

Real
INSMomentumTurbulentViscosityLaplaceForm::computeQpResidual()
{
  Real Cmu = 0.09;
  Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
             _epsilon[_qp];

  return muT * (_grad_u[_qp] * _grad_test[_i][_qp]);
}

Real
INSMomentumTurbulentViscosityLaplaceForm::computeQpJacobian()
{
  // Viscous part, full stress tensor.  The extra contribution comes from the "2"
  // on the diagonal of the viscous stress tensor.
  Real Cmu = 0.09;
  Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
             _epsilon[_qp];

  return muT * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

Real
INSMomentumTurbulentViscosityLaplaceForm::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _kin_var_number)
  {
    Real Cmu = 0.09;
    // Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin,
    // _qp) / _epsilon[_qp];
    Real d_muT_d_kin = _rho[_qp] * Cmu * 2. * computeConcentration(_kin, _qp) *
                       computeConcentrationDerivative(_kin, _phi, _j, _qp) / _epsilon[_qp];

    return d_muT_d_kin * (_grad_u[_qp] * _grad_test[_i][_qp]);
  }

  else if (jvar == _epsilon_var_number)
  {
    Real Cmu = 0.09;
    // Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin,
    // _qp) / _epsilon[_qp];
    Real d_muT_d_epsilon = -_rho[_qp] * Cmu * computeConcentration(_kin, _qp) *
                           computeConcentration(_kin, _qp) * _phi[_j][_qp] /
                           (_epsilon[_qp] * _epsilon[_qp]);

    return d_muT_d_epsilon * (_grad_u[_qp] * _grad_test[_i][_qp]);
  }

  else
    return 0;
}
