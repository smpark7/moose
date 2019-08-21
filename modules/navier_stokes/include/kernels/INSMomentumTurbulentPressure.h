/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMTURBULENTPRESSURE_H
#define INSMOMENTUMTURBULENTPRESSURE_H

#include "Kernel.h"
#include "ScalarTransportBase.h"

// Foward Declarations
class INSMomentumTurbulentPressure;

template <>
InputParameters validParams<INSMomentumTurbulentPressure>();

class INSMomentumTurbulentPressure : public Kernel, public ScalarTransportBase
{
public:
  INSMomentumTurbulentPressure(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  const VariableValue & _kin;
  const VariableGradient & _grad_kin;
  unsigned _kin_id;

  bool _integrate_by_parts;
  unsigned _component;

  const MaterialProperty<Real> & _rho;
};

#endif /* INSMOMENTUMTURBULENTPRESSURE_H */
