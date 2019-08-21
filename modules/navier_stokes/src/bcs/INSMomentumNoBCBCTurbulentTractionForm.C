/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSMomentumNoBCBCTurbulentTractionForm.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<INSMomentumNoBCBCTurbulentTractionForm>()
{
  InputParameters params = validParams<INSMomentumNoBCBCTractionForm>();
  params += validParams<ScalarTransportBase>();

  params.addClassDescription("This class implements the 'No BC' boundary condition for the "
                             "k-epsilon turbulence model based on the "
                             "'traction' form of the viscous stress tensor.");
  params.addRequiredCoupledVar("kin", "The turbulent kinetic energy");
  params.addRequiredCoupledVar("epsilon", "The turbulent dissipation");
  return params;
}

INSMomentumNoBCBCTurbulentTractionForm::INSMomentumNoBCBCTurbulentTractionForm(
    const InputParameters & parameters)
  : INSMomentumNoBCBCTractionForm(parameters),
    ScalarTransportBase(parameters),
    _kin(coupledValue("kin")),
    _epsilon(coupledValue("epsilon")),

    _kin_var_number(coupled("kin")),
    _epsilon_var_number(coupled("epsilon"))
{
}

Real
INSMomentumNoBCBCTurbulentTractionForm::computeQpResidual()
{
  // Compute n . sigma . v, where n is unit normal and v is the test function.
  Real Cmu = 0.09;
  Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
             _epsilon[_qp];
  RealTensorValue sigmaT;

  // First row
  sigmaT(0, 0) = 2. * muT * _grad_u_vel[_qp](0);
  sigmaT(0, 1) = muT * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
  sigmaT(0, 2) = muT * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));

  // Second row
  sigmaT(1, 0) = muT * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1));
  sigmaT(1, 1) = 2. * muT * _grad_v_vel[_qp](1);
  sigmaT(1, 2) = muT * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));

  // Third row
  sigmaT(2, 0) = muT * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2));
  sigmaT(2, 1) = muT * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2));
  sigmaT(2, 2) = 2. * muT * _grad_w_vel[_qp](2);

  // Set up test function
  RealVectorValue test;
  test(_component) = _test[_i][_qp];

  return -_normals[_qp] * (sigmaT * test) + INSMomentumNoBCBCTractionForm::computeQpResidual();
}

Real
INSMomentumNoBCBCTurbulentTractionForm::computeQpJacobian()
{
  Real Cmu = 0.09;
  Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
             _epsilon[_qp];

  // The extra contribution comes from the "2" on the diagonal of the viscous stress tensor
  return -muT * (_grad_phi[_j][_qp] * _normals[_qp] +
                 _grad_phi[_j][_qp](_component) * _normals[_qp](_component)) *
             _test[_i][_qp] +
         INSMomentumNoBCBCTractionForm::computeQpJacobian();
}

Real
INSMomentumNoBCBCTurbulentTractionForm::computeQpOffDiagJacobian(unsigned jvar)
{
  Real Cmu = 0.09;
  Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
             _epsilon[_qp];

  if (jvar == _u_vel_var_number)
    return -muT * _grad_phi[_j][_qp](_component) * _normals[_qp](0) * _test[_i][_qp] +
           INSMomentumNoBCBCTractionForm::computeQpOffDiagJacobian(jvar);

  else if (jvar == _v_vel_var_number)
    return -muT * _grad_phi[_j][_qp](_component) * _normals[_qp](1) * _test[_i][_qp] +
           INSMomentumNoBCBCTractionForm::computeQpOffDiagJacobian(jvar);

  else if (jvar == _w_vel_var_number)
    return -muT * _grad_phi[_j][_qp](_component) * _normals[_qp](2) * _test[_i][_qp] +
           INSMomentumNoBCBCTractionForm::computeQpOffDiagJacobian(jvar);

  else if (jvar == _p_var_number)
    return INSMomentumNoBCBCTractionForm::computeQpOffDiagJacobian(jvar);

  else if (jvar == _kin_var_number)
  {
    // Compute n . sigma . v, where n is unit normal and v is the test function.
    Real Cmu = 0.09;
    // Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin,
    // _qp) / _epsilon[_qp];
    Real d_muT_d_kin =
        _rho[_qp] * Cmu * 2. * computeConcentration(_kin, _qp) * _phi[_j][_qp] / _epsilon[_qp];
    RealTensorValue d_sigmaT_d_kin;

    // First row
    d_sigmaT_d_kin(0, 0) = 2. * d_muT_d_kin * _grad_u_vel[_qp](0);
    d_sigmaT_d_kin(0, 1) = d_muT_d_kin * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
    d_sigmaT_d_kin(0, 2) = d_muT_d_kin * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));

    // Second row
    d_sigmaT_d_kin(1, 0) = d_muT_d_kin * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1));
    d_sigmaT_d_kin(1, 1) = 2. * d_muT_d_kin * _grad_v_vel[_qp](1);
    d_sigmaT_d_kin(1, 2) = d_muT_d_kin * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));

    // Third row
    d_sigmaT_d_kin(2, 0) = d_muT_d_kin * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2));
    d_sigmaT_d_kin(2, 1) = d_muT_d_kin * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2));
    d_sigmaT_d_kin(2, 2) = 2. * d_muT_d_kin * _grad_w_vel[_qp](2);

    // Set up test function
    RealVectorValue test;
    test(_component) = _test[_i][_qp];

    return -_normals[_qp] * (d_sigmaT_d_kin * test);
  }

  else if (jvar == _epsilon_var_number)
  {
    // Compute n . sigma . v, where n is unit normal and v is the test function.
    Real Cmu = 0.09;
    // Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin,
    // _qp) / _epsilon[_qp];
    Real d_muT_d_epsilon = -_rho[_qp] * Cmu * 2. * computeConcentration(_kin, _qp) *
                           computeConcentrationDerivative(_kin, _phi, _j, _qp) /
                           (_epsilon[_qp] * _epsilon[_qp]);
    RealTensorValue d_sigmaT_d_epsilon;

    // First row
    d_sigmaT_d_epsilon(0, 0) = 2. * d_muT_d_epsilon * _grad_u_vel[_qp](0);
    d_sigmaT_d_epsilon(0, 1) = d_muT_d_epsilon * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
    d_sigmaT_d_epsilon(0, 2) = d_muT_d_epsilon * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));

    // Second row
    d_sigmaT_d_epsilon(1, 0) = d_muT_d_epsilon * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1));
    d_sigmaT_d_epsilon(1, 1) = 2. * d_muT_d_epsilon * _grad_v_vel[_qp](1);
    d_sigmaT_d_epsilon(1, 2) = d_muT_d_epsilon * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));

    // Third row
    d_sigmaT_d_epsilon(2, 0) = d_muT_d_epsilon * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2));
    d_sigmaT_d_epsilon(2, 1) = d_muT_d_epsilon * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2));
    d_sigmaT_d_epsilon(2, 2) = 2. * d_muT_d_epsilon * _grad_w_vel[_qp](2);

    // Set up test function
    RealVectorValue test;
    test(_component) = _test[_i][_qp];

    return -_normals[_qp] * (d_sigmaT_d_epsilon * test);
  }

  else
    return 0.;
}
