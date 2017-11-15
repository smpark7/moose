/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSTurbulenceMaterial.h"

template <>
InputParameters
validParams<INSTurbulenceMaterial>()
{
  InputParameters params = validParams<Material>();

  params.addClassDescription("Material class for INS turbulence modeling.");
  params.addRequiredCoupledVar("k", "The turbulent kinetic energy");
  params.addRequiredCoupledVar("epsilon", "The turbulent dissipation");
  params.addParam<Real>("mu", 1, "The dynamic viscosity");
  params.addParam<Real>("rho", 1, "The density");
  return params;
}

INSTurbulenceMaterial::INSTurbulenceMaterial(const InputParameters & parameters)
  : Material(parameters),
    _k(coupledValue("k")),
    _epsilon(coupledValue("epsilon")),
    _grad_k(coupledGradient("k")),
    _grad_epsilon(coupledGradient("epsilon")),
    _user_mu(getParam<Real>("mu")),
    _user_rho(getParam<Real>("rho")),
    _Cmu(0.09),
    _sigk(1.00),
    _sigeps(1.30),
    _C1eps(1.44),
    _C2eps(1.92),
    _mu(declareProperty<Real>("mu")),
    _rho(declareProperty<Real>("rho")),
    _nu(declareProperty<Real>("nu")),

    _nu_turb(declareProperty<Real>("nu_turb")),
    _d_nu_turb_d_k(declareProperty<Real>("d_nu_turb_d_k")),
    _d_nu_turb_d_epsilon(declareProperty<Real>("d_nu_turb_d_epsilon")),
    _grad_nu_turb(declareProperty<RealVectorValue>("grad_nu_turb")),
    _d_grad_nu_turb_d_k(declareProperty<RealVectorValue>("d_grad_nu_turb_d_k")),
    _d_grad_nu_turb_d_grad_k(declareProperty<Real>("d_grad_nu_turb_d_grad_k")),
    _d_grad_nu_turb_d_epsilon(declareProperty<RealVectorValue>("d_grad_nu_turb_d_epsilon")),
    _d_grad_nu_turb_d_grad_epsilon(declareProperty<Real>("d_grad_nu_turb_d_grad_epsilon")),

    _diff_k(declareProperty<Real>("diff_k")),
    _d_diff_k_d_k(declareProperty<Real>("d_diff_k_d_k")),
    _d_diff_k_d_epsilon(declareProperty<Real>("d_diff_k_d_epsilon")),
    _grad_diff_k(declareProperty<RealVectorValue>("grad_diff_k")),
    _d_grad_diff_k_d_k(declareProperty<RealVectorValue>("d_grad_diff_k_d_k")),
    _d_grad_diff_k_d_grad_k(declareProperty<Real>("d_grad_diff_k_d_grad_k")),
    _d_grad_diff_k_d_epsilon(declareProperty<RealVectorValue>("d_grad_diff_k_d_epsilon")),
    _d_grad_diff_k_d_grad_epsilon(declareProperty<Real>("d_grad_diff_k_d_grad_epsilon")),

    _diff_epsilon(declareProperty<Real>("diff_epsilon")),
    _d_diff_epsilon_d_k(declareProperty<Real>("d_diff_epsilon_d_k")),
    _d_diff_epsilon_d_epsilon(declareProperty<Real>("d_diff_epsilon_d_epsilon")),
    _grad_diff_epsilon(declareProperty<RealVectorValue>("grad_diff_epsilon")),
    _d_grad_diff_epsilon_d_k(declareProperty<RealVectorValue>("d_grad_diff_epsilon_d_k")),
    _d_grad_diff_epsilon_d_grad_k(declareProperty<Real>("d_grad_diff_epsilon_d_grad_k")),
    _d_grad_diff_epsilon_d_epsilon(
        declareProperty<RealVectorValue>("d_grad_diff_epsilon_d_epsilon")),
    _d_grad_diff_epsilon_d_grad_epsilon(declareProperty<Real>("d_grad_diff_epsilon_d_grad_epsilon"))
{
}

void
INSTurbulenceMaterial::computeQpProperties()
{
  _mu[_qp] = _user_mu;
  _rho[_qp] = _user_rho;
  _nu[_qp] = _mu[_qp] / _rho[_qp];

  _nu_turb[_qp] = _Cmu * _k[_qp] * _k[_qp] / _epsilon[_qp];
  _d_nu_turb_d_k[_qp] = _Cmu * 2. * _k[_qp] / _epsilon[_qp];
  _d_nu_turb_d_epsilon[_qp] = -_Cmu * _k[_qp] * _k[_qp] / (_epsilon[_qp] * _epsilon[_qp]);
  _grad_nu_turb[_qp] = _Cmu * _k[_qp] / _epsilon[_qp] *
                       (2 * _grad_k[_qp] - _k[_qp] / _epsilon[_qp] * _grad_epsilon[_qp]);
  _d_grad_nu_turb_d_k[_qp] = 2. * _grad_k[_qp] / _epsilon[_qp] -
                             2. * _k[_qp] / (_epsilon[_qp] * _epsilon[_qp]) * _grad_epsilon[_qp];
  _d_grad_nu_turb_d_grad_k[_qp] = 2. * _k[_qp] / _epsilon[_qp];
  _d_grad_nu_turb_d_epsilon[_qp] =
      -2. * _k[_qp] * _grad_k[_qp] / (_epsilon[_qp] * _epsilon[_qp]) +
      2. * _k[_qp] * _k[_qp] * _grad_epsilon[_qp] / std::pow(_epsilon[_qp], 3);
  _d_grad_nu_turb_d_grad_epsilon[_qp] = -_k[_qp] * _k[_qp] / (_epsilon[_qp] * _epsilon[_qp]);

  _diff_k[_qp] = _nu_turb[_qp] / _sigk;
  _d_diff_k_d_k[_qp] = _d_nu_turb_d_k[_qp] / _sigk;
  _d_diff_k_d_epsilon[_qp] = _d_nu_turb_d_epsilon[_qp] / _sigk;
  _grad_diff_k[_qp] = _grad_nu_turb[_qp] / _sigk;
  _d_grad_diff_k_d_k[_qp] = _d_grad_nu_turb_d_k[_qp] / _sigk;
  _d_grad_diff_k_d_grad_k[_qp] = _d_grad_nu_turb_d_grad_k[_qp] / _sigk;
  _d_grad_diff_k_d_epsilon[_qp] = _d_grad_nu_turb_d_epsilon[_qp] / _sigk;
  _d_grad_diff_k_d_grad_epsilon[_qp] = _d_grad_nu_turb_d_grad_epsilon[_qp] / _sigk;

  _diff_epsilon[_qp] = _nu_turb[_qp] / _sigeps;
  _d_diff_epsilon_d_k[_qp] = _d_nu_turb_d_k[_qp] / _sigeps;
  _d_diff_epsilon_d_epsilon[_qp] = _d_nu_turb_d_epsilon[_qp] / _sigeps;
  _grad_diff_epsilon[_qp] = _grad_nu_turb[_qp] / _sigeps;
  _d_grad_diff_epsilon_d_k[_qp] = _d_grad_nu_turb_d_k[_qp] / _sigeps;
  _d_grad_diff_epsilon_d_grad_k[_qp] = _d_grad_nu_turb_d_grad_k[_qp] / _sigeps;
  _d_grad_diff_epsilon_d_epsilon[_qp] = _d_grad_nu_turb_d_epsilon[_qp] / _sigeps;
  _d_grad_diff_epsilon_d_grad_epsilon[_qp] = _d_grad_nu_turb_d_grad_epsilon[_qp] / _sigeps;
}
