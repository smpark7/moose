#include "INSScalarTransportTimeDerivative.h"

template <>
InputParameters
validParams<INSScalarTransportTimeDerivative>()
{
  InputParameters params = validParams<TimeKernel>();
  params += validParams<ScalarTransportBase>();
  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
  params.addParam<MaterialPropertyName>("rho_name", "rho", "Name of the density.");
  return params;
}

INSScalarTransportTimeDerivative::INSScalarTransportTimeDerivative(
    const InputParameters & parameters)
  : TimeKernel(parameters),
    ScalarTransportBase(parameters),
    _lumping(getParam<bool>("lumping")),
    _rho(getMaterialProperty<Real>("rho_name"))
{
}

Real
INSScalarTransportTimeDerivative::computeQpResidual()
{
  return _test[_i][_qp] * _rho[_qp] * computeConcentrationDot(_u, _u_dot, _qp);
}

Real
INSScalarTransportTimeDerivative::computeQpJacobian()
{
  return _test[_i][_qp] * _rho[_qp] *
         computeConcentrationDotDerivative(_u, _u_dot, _du_dot_du, _phi, _j, _qp);
}
