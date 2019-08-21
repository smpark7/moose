#ifndef INSEpsilon_H
#define INSEpsilon_H

#include "Kernel.h"
#include "ScalarTransportBase.h"

// Forward Declarations
class INSEpsilon;

template <>
InputParameters validParams<INSEpsilon>();

/**
 * This class computes the epsilon residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation with a standard k-epsilon turbulence model.
 */
class INSEpsilon : public Kernel, public ScalarTransportBase
{
public:
  INSEpsilon(const InputParameters & parameters);

  virtual ~INSEpsilon() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  // Coupled variables
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;
  const VariableValue & _kin;

  // Gradients
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;
  const VariableGradient & _grad_kin;

  // Variable numberings
  unsigned int _u_vel_var_number;
  unsigned int _v_vel_var_number;
  unsigned int _w_vel_var_number;
  unsigned int _kin_var_number;

  // Material properties
  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _rho;

  // Parameters
  Real _Cmu;
  Real _sigk;
  Real _sigeps;
  Real _C1eps;
  Real _C2eps;
};

#endif // INSEpsilon_H
