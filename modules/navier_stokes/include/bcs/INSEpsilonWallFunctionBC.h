/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSEPSILONWALLFUNCTIONBC_H
#define INSEPSILONWALLFUNCTIONBC_H

#include "IntegratedBC.h"
#include "ScalarTransportBase.h"

// Forward Declarations
class INSEpsilonWallFunctionBC;

template <>
InputParameters validParams<INSEpsilonWallFunctionBC>();

/**
 * Implements a wall function for the boundary tangential component of the momentum equation. This
 * boundary condition is applicable when using the k-epsilon turbulence model.
 */
class INSEpsilonWallFunctionBC : public IntegratedBC, public ScalarTransportBase
{
public:
  INSEpsilonWallFunctionBC(const InputParameters & parameters);

  virtual ~INSEpsilonWallFunctionBC() {}

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // Coupled variables
  const VariableValue & _kin;

  // Variable numberings
  unsigned _kin_var_number;

  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _rho;
};

#endif
