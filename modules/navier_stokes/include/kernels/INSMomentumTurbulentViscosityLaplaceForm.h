/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMTURBULENTVISCOSITYLAPLACEFORM_H
#define INSMOMENTUMTURBULENTVISCOSITYLAPLACEFORM_H

#include "Kernel.h"
#include "ScalarTransportBase.h"

// Forward Declarations
class INSMomentumTurbulentViscosityLaplaceForm;

template <>
InputParameters validParams<INSMomentumTurbulentViscosityLaplaceForm>();

/**
 * This class computes the residual and Jacobians for the turbulent viscosity term in the
 * INS momentum component equations for the k-epsilon turbulence model.
 */
class INSMomentumTurbulentViscosityLaplaceForm : public Kernel, public ScalarTransportBase
{
public:
  INSMomentumTurbulentViscosityLaplaceForm(const InputParameters & parameters);

  virtual ~INSMomentumTurbulentViscosityLaplaceForm() {}

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
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
  unsigned _kin_var_number;
  unsigned _epsilon_var_number;

  // Parameters
  unsigned _component;

  // Material properties
  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _rho;
};

#endif
