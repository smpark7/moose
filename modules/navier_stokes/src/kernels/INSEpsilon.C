#include "INSEpsilon.h"

template <>
InputParameters
validParams<INSEpsilon>()
{
  InputParameters params = validParams<INSScalar>();
  params.addRequiredCoupledVar("kin", "turbulent kinetic energy");
  return params;
}

INSEpsilon::INSEpsilon(const InputParameters & parameters)
  : INSScalar(parameters),
    _kin(coupledValue("kin")),
    _grad_kin(coupledGradient("kin")),
    _kin_var_number(coupled("kin")),
    _diff_epsilon(getMaterialProperty<Real>("diff_epsilon")),
    _d_diff_epsilon_d_epsilon(getMaterialProperty<Real>("d_diff_epsilon_d_epsilon")),
    _grad_diff_epsilon(getMaterialProperty<RealVectorValue>("grad_diff_epsilon")),
    _d_grad_diff_epsilon_d_epsilon(
        getMaterialProperty<RealVectorValue>("d_grad_diff_epsilon_d_epsilon")),
    _d_grad_diff_epsilon_d_grad_epsilon(
        getMaterialProperty<Real>("d_grad_diff_epsilon_d_grad_epsilon")),
    _nu_turb(getMaterialProperty<Real>("nu_turb")),
    _d_nu_turb_d_epsilon(getMaterialProperty<Real>("d_nu_turb_d_epsilon")),
    _C1eps(1.44),
    _C2eps(1.92)
{
}

Real
INSEpsilon::diffusivity()
{
  return _diff_epsilon[_qp];
}

RealVectorValue
INSEpsilon::grad_diffusivity()
{
  return _grad_diff_epsilon[_qp];
}

Real
INSEpsilon::d_diff_d_uj()
{
  return _d_diff_epsilon_d_epsilon[_qp] * _phi[_j][_qp];
}

RealVectorValue
INSEpsilon::d_grad_diff_d_uj()
{
  return _d_grad_diff_epsilon_d_epsilon[_qp] * _phi[_j][_qp] +
         _d_grad_diff_epsilon_d_grad_epsilon[_qp] * _grad_phi[_j][_qp];
}

Real
INSEpsilon::strongSink()
{
  // Turbulent dissipative sink part
  Real dissipative_part = _C2eps * _u[_qp] * _u[_qp] / _kin[_qp];

  // Source term from velocity gradients
  Real strain_tensor_double_dot_product = 0.;
  for (int m = 0; m < 3; m++)
  {
    for (int n = 0; n < 3; n++)
    {
      Real Em;
      switch (m)
      {
        case 0:
          Em = _grad_u_vel[_qp](n);
          break;
        case 1:
          Em = _grad_v_vel[_qp](n);
          break;
        case 2:
          Em = _grad_w_vel[_qp](n);
          break;
        default:
          mooseError("Unrecognized index");
      }
      Real En;
      switch (n)
      {
        case 0:
          En = _grad_u_vel[_qp](m);
          break;
        case 1:
          En = _grad_v_vel[_qp](m);
          break;
        case 2:
          En = _grad_w_vel[_qp](m);
          break;
        default:
          mooseError("Unrecognized index");
      }
      Real Emn = (Em + En);
      strain_tensor_double_dot_product += Emn * Emn;
    }
  }

  Real source_part =
      -_C1eps * 0.5 * _nu_turb[_qp] * _u[_qp] / _kin[_qp] * strain_tensor_double_dot_product;

  return dissipative_part + source_part;
}

Real
INSEpsilon::strongSinkJac()
{
  // Turbulent dissipative sink part
  Real dissipative_part = 2. * _C2eps * _u[_qp] * _phi[_j][_qp] / _kin[_qp];

  // Source term from velocity gradients
  Real strain_tensor_double_dot_product = 0.;
  for (int m = 0; m < 3; m++)
  {
    for (int n = 0; n < 3; n++)
    {
      Real Em;
      switch (m)
      {
        case 0:
          Em = _grad_u_vel[_qp](n);
          break;
        case 1:
          Em = _grad_v_vel[_qp](n);
          break;
        case 2:
          Em = _grad_w_vel[_qp](n);
          break;
        default:
          mooseError("Unrecognized index");
      }
      Real En;
      switch (n)
      {
        case 0:
          En = _grad_u_vel[_qp](m);
          break;
        case 1:
          En = _grad_v_vel[_qp](m);
          break;
        case 2:
          En = _grad_w_vel[_qp](m);
          break;
        default:
          mooseError("Unrecognized index");
      }
      Real Emn = (Em + En);
      strain_tensor_double_dot_product += Emn * Emn;
    }
  }

  Real source_part =
      -_C1eps * 0.5 / _kin[_qp] * strain_tensor_double_dot_product *
      (_d_nu_turb_d_epsilon[_qp] * _phi[_j][_qp] * _u[_qp] + _nu_turb[_qp] * _phi[_j][_qp]);

  return dissipative_part + source_part;
}
