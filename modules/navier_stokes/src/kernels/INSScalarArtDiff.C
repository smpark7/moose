/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSScalarArtDiff.h"

template <>
InputParameters
validParams<INSScalarArtDiff>()
{
  InputParameters params = validParams<Kernel>();
  params += validParams<ScalarTransportBase>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D

  // Optional parameters
  params.addParam<MaterialPropertyName>("mu_name", "mu", "The name of the dynamic viscosity");
  params.addParam<MaterialPropertyName>("rho_name", "rho", "The name of the density");

  return params;
}

INSScalarArtDiff::INSScalarArtDiff(const InputParameters & parameters)
  : Kernel(parameters),
    ScalarTransportBase(parameters),

    // Coupled variables
    _u_vel(coupledValue("u")),
    _v_vel(coupledValue("v")),
    _w_vel(coupledValue("w")),

    // Gradients
    _grad_u_vel(coupledGradient("u")),
    _grad_v_vel(coupledGradient("v")),
    _grad_w_vel(coupledGradient("w")),

    // Variable numberings
    _u_vel_var_number(coupled("u")),
    _v_vel_var_number(coupled("v")),
    _w_vel_var_number(coupled("w")),

    // Material properties
    _mu(getMaterialProperty<Real>("mu_name")),
    _rho(getMaterialProperty<Real>("rho_name"))
{
}

Real
INSScalarArtDiff::computeQpResidual()
{
  // The convection part, rho * (u.grad) * u_component * v.
  // Note: _grad_u is the gradient of the _component entry of the velocity vector.

  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  // Real Pe = _rho[_qp] * _current_elem->hmax() * U.norm() / (2. * _mu[_qp]);
  // Real delta = 1. / std::tanh(Pe) - 1. / Pe;
  // Real alpha0 = delta * U.norm() * _current_elem->hmax() / 2.;
  // Real PG_test = alpha0 / (U.norm() * U.norm()) * U * _grad_test[_i][_qp];
  // Real convective_part = _rho[_qp] * U * computeConcentrationGradient(_u, _grad_u, _qp) *
  // PG_test;

  return _current_elem->hmax() * U.norm() / 2. * _grad_test[_i][_qp] *
         computeConcentrationGradient(_u, _grad_u, _qp);
}

Real
INSScalarArtDiff::computeQpJacobian()
{
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
  // Real Pe = _rho[_qp] * _current_elem->hmax() * U.norm() / (2. * _mu[_qp]);
  // Real delta = 1. / std::tanh(Pe) - 1. / Pe;
  // Real alpha0 = delta * U.norm() * _current_elem->hmax() / 2.;
  // Real PG_test = alpha0 / (U.norm() * U.norm()) * U * _grad_test[_i][_qp];

  // // Convective part
  // Real convective_part =
  //     _rho[_qp] * U *
  //     computeConcentrationGradientDerivative(_u, _grad_u, _phi, _grad_phi, _j, _qp) * PG_test;

  return _current_elem->hmax() * U.norm() / 2. * _grad_test[_i][_qp] *
         computeConcentrationGradientDerivative(_u, _grad_u, _phi, _grad_phi, _j, _qp);
}

Real
INSScalarArtDiff::computeQpOffDiagJacobian(unsigned jvar)
{
  // In Stokes/Laplacian version, off-diag Jacobian entries wrt u,v,w are zero
  if (jvar == _u_vel_var_number)
  {
    RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
    // Real Pe = _rho[_qp] * _current_elem->hmax() * U.norm() / (2. * _mu[_qp]);
    // Real delta = 1. / std::tanh(Pe) - 1. / Pe;
    // Real alpha0 = delta * U.norm() * _current_elem->hmax() / 2.;
    // Real PG_test = alpha0 / (U.norm() * U.norm()) * U * _grad_test[_i][_qp];

    // Real convective_part =
    //     _phi[_j][_qp] * computeConcentrationGradient(_u, _grad_u, _qp)(0) * PG_test;

    return _current_elem->hmax() * _u_vel[_qp] * _phi[_j][_qp] / U.norm() / 2. *
           _grad_test[_i][_qp] * computeConcentrationGradient(_u, _grad_u, _qp);
  }

  else if (jvar == _v_vel_var_number)
  {
    RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
    // Real Pe = _rho[_qp] * _current_elem->hmax() * U.norm() / (2. * _mu[_qp]);
    // Real delta = 1. / std::tanh(Pe) - 1. / Pe;
    // Real alpha0 = delta * U.norm() * _current_elem->hmax() / 2.;
    // Real PG_test = alpha0 / (U.norm() * U.norm()) * U * _grad_test[_i][_qp];

    // Real convective_part =
    //     _phi[_j][_qp] * computeConcentrationGradient(_u, _grad_u, _qp)(1) * PG_test;

    return _current_elem->hmax() * _v_vel[_qp] * _phi[_j][_qp] / U.norm() / 2. *
           _grad_test[_i][_qp] * computeConcentrationGradient(_u, _grad_u, _qp);
  }

  else if (jvar == _w_vel_var_number)
  {
    RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);
    // Real Pe = _rho[_qp] * _current_elem->hmax() * U.norm() / (2. * _mu[_qp]);
    // Real delta = 1. / std::tanh(Pe) - 1. / Pe;
    // Real alpha0 = delta * U.norm() * _current_elem->hmax() / 2.;
    // Real PG_test = alpha0 / (U.norm() * U.norm()) * U * _grad_test[_i][_qp];

    // Real convective_part =
    //     _phi[_j][_qp] * computeConcentrationGradient(_u, _grad_u, _qp)(2) * PG_test;

    return _current_elem->hmax() * _w_vel[_qp] * _phi[_j][_qp] / U.norm() / 2. *
           _grad_test[_i][_qp] * computeConcentrationGradient(_u, _grad_u, _qp);
  }

  else
    return 0;
}
