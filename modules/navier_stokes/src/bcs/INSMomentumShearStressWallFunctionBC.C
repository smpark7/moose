/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSMomentumShearStressWallFunctionBC.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<INSMomentumShearStressWallFunctionBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params += validParams<ScalarTransportBase>();

  params.addClassDescription("Tangential momentum BC for k-epsilon turbulence model.");
  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("p", "pressure");
  params.addRequiredCoupledVar("kin", "The turbulent kinetic energy");
  params.addRequiredCoupledVar("epsilon", "The turbulent dissipation");

  // Required parameters
  params.addRequiredParam<unsigned>(
      "component",
      "0,1,2 depending on if we are solving the x,y,z component of the momentum equation");
  params.addParam<bool>("integrate_p_by_parts",
                        true,
                        "Allows simulations to be run with pressure BC if set to false");

  // Optional parameters
  params.addParam<MaterialPropertyName>("mu_name", "mu", "The name of the dynamic viscosity");
  params.addParam<MaterialPropertyName>("rho_name", "rho", "The name of the density");
  params.addParam<bool>(
      "add_iso_art_diff", false, "Whether to add isotropic artificial diffusion.");

  return params;
}

INSMomentumShearStressWallFunctionBC::INSMomentumShearStressWallFunctionBC(
    const InputParameters & parameters)
  : IntegratedBC(parameters),
    ScalarTransportBase(parameters),

    // Coupled variables
    _u_vel(coupledValue("u")),
    _v_vel(coupledValue("v")),
    _w_vel(coupledValue("w")),
    _p(coupledValue("p")),
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
    _p_var_number(coupled("p")),
    _kin_var_number(coupled("kin")),
    _epsilon_var_number(coupled("epsilon")),

    // Required parameters
    _component(getParam<unsigned>("component")),
    _integrate_p_by_parts(getParam<bool>("integrate_p_by_parts")),

    // Material properties
    _mu(getMaterialProperty<Real>("mu_name")),
    _rho(getMaterialProperty<Real>("rho_name")),
    _add_iso_art_diff(getParam<bool>("add_iso_art_diff"))
{
}

Real
INSMomentumShearStressWallFunctionBC::computeQpResidual()
{
  // Compute (n . sigma . n) . n . v, where n is unit normal and v is the test function.
  RealTensorValue rateOfStrain, sigma, sigmaT;

  // First row
  rateOfStrain(0, 0) = _grad_u_vel[_qp](0);
  rateOfStrain(0, 1) = 0.5 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
  rateOfStrain(0, 2) = 0.5 * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));

  // Second row
  rateOfStrain(1, 0) = 0.5 * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1));
  rateOfStrain(1, 1) = _grad_v_vel[_qp](1);
  rateOfStrain(1, 2) = 0.5 * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));

  // Third row
  rateOfStrain(2, 0) = 0.5 * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2));
  rateOfStrain(2, 1) = 0.5 * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2));
  rateOfStrain(2, 2) = _grad_w_vel[_qp](2);

  sigma = rateOfStrain * _mu[_qp] * 2.;
  if (_add_iso_art_diff)
    sigma += rateOfStrain * _current_elem->hmax() *
             (RealVectorValue(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp])).norm() / 2. * 2.;

  // If the pressure term is integrated by parts, it is part of the
  // no-BC-BC, otherwise, it is not.
  if (_integrate_p_by_parts)
  {
    sigma(0, 0) -= _p[_qp];
    sigma(1, 1) -= _p[_qp];
    sigma(2, 2) -= _p[_qp];
  }

  // Set up test function
  RealVectorValue test;
  test(_component) = _test[_i][_qp];

  Real normal_stress_component = -_normals[_qp] * (sigma * _normals[_qp]) * _normals[_qp] * test;

  // Compute tangential stress component
  Real Cmu = 0.09, yStarPlus = 11.06;
  RealVectorValue U = {_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]};
  Real uTau =
      std::max(std::pow(Cmu, 0.25) * std::sqrt(std::max(computeConcentration(_kin, _qp), 0.)),
               U.norm() / yStarPlus);
  Real tangential_stress_component = uTau / yStarPlus * U * test;

  // Compute turbulent stress component
  Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
             _epsilon[_qp];
  sigmaT = rateOfStrain * muT * 2.;

  Real turbulent_stress_component = -_normals[_qp] * (sigmaT * test);

  return normal_stress_component + turbulent_stress_component + tangential_stress_component;
}

Real
INSMomentumShearStressWallFunctionBC::computeQpJacobian()
{
  // normal stress component
  Real normal_stress_component = -2. * _mu[_qp] * _normals[_qp](_component) *
                                 _normals[_qp](_component) * _normals[_qp] * _grad_phi[_j][_qp] *
                                 _test[_i][_qp];

  // tangential stress component
  Real Cmu = 0.09, yStarPlus = 11.06;
  RealVectorValue U = {_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]};
  Real uTau, d_uTau_d_var;
  Real uTauOpt1 = std::pow(Cmu, 0.25) * std::sqrt(std::max(computeConcentration(_kin, _qp), 0.));
  Real uTauOpt2 = U.norm() / yStarPlus;
  if (uTauOpt1 > uTauOpt2)
  {
    uTau = uTauOpt1;
    d_uTau_d_var = 0;
  }
  else
  {
    uTau = uTauOpt2;
    d_uTau_d_var = 1. / yStarPlus * _phi[_j][_qp] * U(_component) / U.norm();
  }
  Real tangential_stress_component = uTau / yStarPlus * _phi[_j][_qp] * _test[_i][_qp] +
                                     d_uTau_d_var / yStarPlus * U(_component) * _test[_i][_qp];

  // turbulent stress component
  // The extra contribution comes from the "2" on the diagonal of the viscous stress tensor
  Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
             _epsilon[_qp];
  Real turbulent_stress_component =
      -muT * _test[_i][_qp] * (_grad_phi[_j][_qp] * _normals[_qp] +
                               _grad_phi[_j][_qp](_component) * _normals[_qp](_component));

  return normal_stress_component + tangential_stress_component + turbulent_stress_component;
}

Real
INSMomentumShearStressWallFunctionBC::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    // normal stress component
    Real normal_stress_component = -2. * _mu[_qp] * _normals[_qp](0) * _normals[_qp](_component) *
                                   _test[_i][_qp] * _normals[_qp] * _grad_phi[_j][_qp];

    // tangential stress component
    Real Cmu = 0.09, yStarPlus = 11.06;
    RealVectorValue U = {_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]};
    Real uTau, d_uTau_d_u_vel;
    Real uTauOpt1 = std::pow(Cmu, 0.25) * std::sqrt(std::max(computeConcentration(_kin, _qp), 0.));
    Real uTauOpt2 = U.norm() / yStarPlus;
    if (uTauOpt1 > uTauOpt2)
    {
      uTau = uTauOpt1;
      d_uTau_d_u_vel = 0;
    }
    else
    {
      uTau = uTauOpt2;
      d_uTau_d_u_vel = 1. / yStarPlus * _phi[_j][_qp] * U(0) / U.norm();
    }

    Real tangential_stress_component = d_uTau_d_u_vel / yStarPlus * U(_component) * _test[_i][_qp];

    // turbulent stress component
    Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
               _epsilon[_qp];
    Real turbulent_stress_component =
        -muT * _grad_phi[_j][_qp](_component) * _normals[_qp](0) * _test[_i][_qp];

    return normal_stress_component + tangential_stress_component + turbulent_stress_component;
  }

  if (jvar == _v_vel_var_number)
  {
    // normal stress component
    Real normal_stress_component = -2. * _mu[_qp] * _normals[_qp](1) * _normals[_qp](_component) *
                                   _test[_i][_qp] * _normals[_qp] * _grad_phi[_j][_qp];

    // tangential stress component
    Real Cmu = 0.09, yStarPlus = 11.06;
    RealVectorValue U = {_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]};
    Real uTau, d_uTau_d_v_vel;
    Real uTauOpt1 = std::pow(Cmu, 0.25) * std::sqrt(std::max(computeConcentration(_kin, _qp), 0.));
    Real uTauOpt2 = U.norm() / yStarPlus;
    if (uTauOpt1 > uTauOpt2)
    {
      uTau = uTauOpt1;
      d_uTau_d_v_vel = 0;
    }
    else
    {
      uTau = uTauOpt2;
      d_uTau_d_v_vel = 1. / yStarPlus * _phi[_j][_qp] * U(1) / U.norm();
    }

    Real tangential_stress_component = d_uTau_d_v_vel / yStarPlus * U(_component) * _test[_i][_qp];

    // turbulent stress component
    Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
               _epsilon[_qp];
    Real turbulent_stress_component =
        -muT * _grad_phi[_j][_qp](_component) * _normals[_qp](1) * _test[_i][_qp];

    return normal_stress_component + tangential_stress_component + turbulent_stress_component;
  }

  if (jvar == _w_vel_var_number)
  {
    // normal stress component
    Real normal_stress_component = -2. * _mu[_qp] * _normals[_qp](2) * _normals[_qp](_component) *
                                   _test[_i][_qp] * _normals[_qp] * _grad_phi[_j][_qp];

    // tangential stress component
    Real Cmu = 0.09, yStarPlus = 11.06;
    RealVectorValue U = {_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]};
    Real uTau, d_uTau_d_w_vel;
    Real uTauOpt1 = std::pow(Cmu, 0.25) * std::sqrt(std::max(computeConcentration(_kin, _qp), 0.));
    Real uTauOpt2 = U.norm() / yStarPlus;
    if (uTauOpt1 > uTauOpt2)
    {
      uTau = uTauOpt1;
      d_uTau_d_w_vel = 0;
    }
    else
    {
      uTau = uTauOpt2;
      d_uTau_d_w_vel = 1. / yStarPlus * _phi[_j][_qp] * U(2) / U.norm();
    }

    Real tangential_stress_component = d_uTau_d_w_vel / yStarPlus * U(_component) * _test[_i][_qp];

    // turbulent stress component
    Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
               _epsilon[_qp];
    Real turbulent_stress_component =
        -muT * _grad_phi[_j][_qp](_component) * _normals[_qp](2) * _test[_i][_qp];

    return normal_stress_component + tangential_stress_component + turbulent_stress_component;
  }

  else if (jvar == _p_var_number)
  {
    if (_integrate_p_by_parts)
      return _phi[_j][_qp] * _normals[_qp](_component) * _test[_i][_qp];
    else
      return 0.;
  }

  else if (jvar == _kin_var_number)
  {
    // Compute tangential stress component
    Real Cmu = 0.09, yStarPlus = 11.06;
    RealVectorValue U = {_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]};
    Real uTau, d_uTau_d_kin;
    Real uTauOpt1 = std::pow(Cmu, 0.25) * std::sqrt(std::max(computeConcentration(_kin, _qp), 0.));
    Real uTauOpt2 = U.norm() / yStarPlus;
    if (uTauOpt1 > uTauOpt2)
    {
      uTau = uTauOpt1;
      d_uTau_d_kin = 0.5 * std::pow(Cmu, 0.25) /
                     std::sqrt(std::max(computeConcentration(_kin, _qp), 0.)) * _phi[_j][_qp];
    }
    else
    {
      uTau = uTauOpt2;
      d_uTau_d_kin = 0.;
    }

    Real tangential_stress_component = d_uTau_d_kin / yStarPlus * U(_component) * _test[_i][_qp];

    // Compute turbulent stress component
    RealTensorValue rateOfStrain;

    // First row
    rateOfStrain(0, 0) = _grad_u_vel[_qp](0);
    rateOfStrain(0, 1) = 0.5 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
    rateOfStrain(0, 2) = 0.5 * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));

    // Second row
    rateOfStrain(1, 0) = 0.5 * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1));
    rateOfStrain(1, 1) = _grad_v_vel[_qp](1);
    rateOfStrain(1, 2) = 0.5 * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));

    // Third row
    rateOfStrain(2, 0) = 0.5 * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2));
    rateOfStrain(2, 1) = 0.5 * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2));
    rateOfStrain(2, 2) = _grad_w_vel[_qp](2);

    // Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin,
    // _qp) / _epsilon[_qp];
    Real d_muT_d_kin = _rho[_qp] * Cmu * 2. * computeConcentration(_kin, _qp) *
                       computeConcentrationDerivative(_kin, _phi, _j, _qp) / _epsilon[_qp];
    RealTensorValue d_sigmaT_d_kin = rateOfStrain * d_muT_d_kin * 2.;

    // Set up test function
    RealVectorValue test;
    test(_component) = _test[_i][_qp];

    Real turbulent_stress_component = -_normals[_qp] * (d_sigmaT_d_kin * test);

    return tangential_stress_component + turbulent_stress_component;
  }

  else if (jvar == _epsilon_var_number)
  {
    // Compute turbulent stress component
    RealTensorValue rateOfStrain;

    // First row
    rateOfStrain(0, 0) = _grad_u_vel[_qp](0);
    rateOfStrain(0, 1) = 0.5 * (_grad_u_vel[_qp](1) + _grad_v_vel[_qp](0));
    rateOfStrain(0, 2) = 0.5 * (_grad_u_vel[_qp](2) + _grad_w_vel[_qp](0));

    // Second row
    rateOfStrain(1, 0) = 0.5 * (_grad_v_vel[_qp](0) + _grad_u_vel[_qp](1));
    rateOfStrain(1, 1) = _grad_v_vel[_qp](1);
    rateOfStrain(1, 2) = 0.5 * (_grad_v_vel[_qp](2) + _grad_w_vel[_qp](1));

    // Third row
    rateOfStrain(2, 0) = 0.5 * (_grad_w_vel[_qp](0) + _grad_u_vel[_qp](2));
    rateOfStrain(2, 1) = 0.5 * (_grad_w_vel[_qp](1) + _grad_v_vel[_qp](2));
    rateOfStrain(2, 2) = _grad_w_vel[_qp](2);

    Real Cmu = 0.09;
    // Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin,
    // _qp) / _epsilon[_qp];
    Real d_muT_d_epsilon = -_rho[_qp] * Cmu * computeConcentration(_kin, _qp) *
                           computeConcentration(_kin, _qp) * _phi[_j][_qp] /
                           (_epsilon[_qp] * _epsilon[_qp]);
    RealTensorValue d_sigmaT_d_epsilon = rateOfStrain * d_muT_d_epsilon * 2.;

    // Set up test function
    RealVectorValue test;
    test(_component) = _test[_i][_qp];

    Real turbulent_stress_component = -_normals[_qp] * (d_sigmaT_d_epsilon * test);

    return turbulent_stress_component;
  }

  else
    return 0.;
}
