/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMTURBULENTPRESSURERZ_H
#define INSMOMENTUMTURBULENTPRESSURERZ_H

#include "INSMomentumTurbulentPressure.h"

// Foward Declarations
class INSMomentumTurbulentPressureRZ;

template <>
InputParameters validParams<INSMomentumTurbulentPressureRZ>();

class INSMomentumTurbulentPressureRZ : public INSMomentumTurbulentPressure
{
public:
  INSMomentumTurbulentPressureRZ(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;
};

#endif /* INSMOMENTUMTURBULENTPRESSURERZ_H */
