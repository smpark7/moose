/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMSHEARSTRESSWALLFUNCTIONBC_H
#define INSMOMENTUMSHEARSTRESSWALLFUNCTIONBC_H

#include "IntegratedBC.h"
#include "ScalarTransportBase.h"

// Forward Declarations
class INSMomentumShearStressWallFunctionBC;

template <>
InputParameters validParams<INSMomentumShearStressWallFunctionBC>();

/**
 * Implements a wall function for the boundary tangential component of the momentum equation. This
 * boundary condition is applicable when using the k-epsilon turbulence model.
 */
class INSMomentumShearStressWallFunctionBC : public IntegratedBC, public ScalarTransportBase
{
public:
  INSMomentumShearStressWallFunctionBC(const InputParameters & parameters);

  virtual ~INSMomentumShearStressWallFunctionBC() {}

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _p;
  const VariableValue & _kin;
  const VariableValue & _epsilon;

  // Gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;

  // Variable numberings
  unsigned _u_vel_var_number;
  unsigned _v_vel_var_number;
  unsigned _w_vel_var_number;
  unsigned _p_var_number;
  unsigned _kin_var_number;
  unsigned _epsilon_var_number;

  unsigned _component;
  bool _integrate_p_by_parts;

  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _rho;
  bool _add_iso_art_diff;
};

#endif
