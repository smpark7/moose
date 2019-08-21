/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSMomentumTurbulentPressureRZ.h"

template <>
InputParameters
validParams<INSMomentumTurbulentPressureRZ>()
{
  InputParameters params = validParams<INSMomentumTurbulentPressure>();

  params.addClassDescription(
      "Models the turbulent pressure term of the INS momentum equation in RZ coordinates.");
  return params;
}

INSMomentumTurbulentPressureRZ::INSMomentumTurbulentPressureRZ(const InputParameters & parameters)
  : INSMomentumTurbulentPressure(parameters)
{
}

Real
INSMomentumTurbulentPressureRZ::computeQpResidual()
{
  Real res_base = INSMomentumTurbulentPressure::computeQpResidual();

  if (_component == 0 && _integrate_by_parts)
  {
    const Real r = _q_point[_qp](0);
    res_base += -_test[_i][_qp] / r * 2. / 3. * _rho[_qp] * computeConcentration(_kin, _qp);
  }

  return res_base;
}

Real
INSMomentumTurbulentPressureRZ::computeQpJacobian()
{
  return 0.;
}

Real
INSMomentumTurbulentPressureRZ::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _kin_id)
  {
    Real jac_base = INSMomentumTurbulentPressure::computeQpOffDiagJacobian(jvar);

    if (_component == 0 && _integrate_by_parts)
    {
      const Real r = _q_point[_qp](0);
      jac_base += -_test[_i][_qp] / r * 2. / 3. * _rho[_qp] * _phi[_j][_qp];
    }

    return jac_base;
  }

  else
    return 0;
}
