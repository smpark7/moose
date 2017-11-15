#include "INSScalar.h"
#include "MooseVariable.h"

template <>
InputParameters
validParams<INSScalar>()
{
  InputParameters params = validParams<INSBase>();
  params.addParam<bool>("supg", true, "Whether to perform SUPG stabilization of the k residual");
  return params;
}

INSScalar::INSScalar(const InputParameters & parameters)
  : INSBase(parameters), _second_u(_var.secondSln()), _supg(getParam<bool>("supg"))
{
}

Real
INSScalar::diffusivity()
{
  mooseError("You must implement a derivative class that implements the diffusivity() method, e.g. "
             "it computes or returns a material property.");

  return 0;
}

RealVectorValue
INSScalar::grad_diffusivity()
{
  return RealVectorValue(0, 0, 0);
}

Real
INSScalar::d_diff_d_uj()
{
  return 0;
}

RealVectorValue
INSScalar::d_grad_diff_d_uj()
{
  return RealVectorValue(0, 0, 0);
}

Real
INSScalar::d_tau_d_uj()
{
  return d_tau_d_diff() * d_diff_d_uj();
}

Real
INSScalar::strongConvection()
{
  return a() * _grad_u[_qp];
}

Real
INSScalar::strongConvectionJac()
{
  return a() * _grad_phi[_j][_qp];
}

Real
INSScalar::strongDiffusion()
{
  return -(diffusivity() * _second_u[_qp].tr() + grad_diffusivity() * _grad_u[_qp]);
}

Real
INSScalar::strongDiffusionJac()
{
  return -(d_diff_d_uj() * _second_u[_qp].tr() + diffusivity() * _second_phi[_j][_qp].tr() +
           d_grad_diff_d_uj() * _grad_u[_qp] + grad_diffusivity() * _grad_phi[_j][_qp]);
}

Real
INSScalar::weakDiffusion()
{
  return _grad_test[_i][_qp] * diffusivity() * _grad_u[_qp];
}

Real
INSScalar::weakDiffusionJac()
{
  return _grad_test[_i][_qp] * (d_diff_d_uj() * _grad_u[_qp] + diffusivity() * _grad_phi[_j][_qp]);
}

Real
INSScalar::strongSink()
{
  return 0;
}

Real
INSScalar::strongSinkJac()
{
  return 0;
}

Real
INSScalar::computeQpResidual()
{
  return _test[_i][_qp] * (strongConvection() + strongSink()) + weakDiffusion() +
         (_supg ? computeQpPGResidual() : 0);
}

Real
INSScalar::computeQpPGResidual()
{
  return tau() * a() * _grad_test[_i][_qp] *
         (strongConvection() + strongDiffusion() + strongSink());
}

Real
INSScalar::computeQpJacobian()
{
  return _test[_i][_qp] * (strongConvectionJac() + strongSinkJac()) + weakDiffusionJac() +
         (_supg ? computeQpPGJacobian() : 0);
}

Real
INSScalar::computeQpPGJacobian()
{
  return d_tau_d_uj() * a() * _grad_test[_i][_qp] *
             (strongConvection() + strongDiffusion() + strongSink()) +
         tau() * a() * _grad_test[_i][_qp] *
             (strongConvectionJac() + strongDiffusionJac() + strongSinkJac());
}

Real
INSScalar::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _u_vel_var_number)
  {
    Real jac = _test[_i][_qp] * _phi[_j][_qp] * _grad_u[_qp](0);
    if (_supg)
      jac += computeQpPGOffDiagJacobian(0);
    return jac;
  }

  else if (jvar == _v_vel_var_number)
  {
    Real jac = _test[_i][_qp] * _phi[_j][_qp] * _grad_u[_qp](1);
    if (_supg)
      jac += computeQpPGOffDiagJacobian(1);
    return jac;
  }

  else if (jvar == _w_vel_var_number)
  {
    Real jac = _test[_i][_qp] * _phi[_j][_qp] * _grad_u[_qp](2);
    if (_supg)
      jac += computeQpPGOffDiagJacobian(2);
    return jac;
  }

  else
    return 0.0;
}

Real
INSScalar::computeQpPGOffDiagJacobian(unsigned comp)
{
  RealVectorValue d_a_d_a_comp(0, 0, 0);
  d_a_d_a_comp(comp) = _phi[_j][_qp];

  return dTauDUComp(comp) * a() * _grad_test[_i][_qp] *
             (strongConvection() + strongDiffusion() + strongSink()) +
         tau() * d_a_d_a_comp * _grad_test[_i][_qp] *
             (strongConvection() + strongDiffusion() + strongSink()) +
         tau() * a() * _grad_test[_i][_qp] * _phi[_j][_qp] * _grad_u[_qp](comp);
}
