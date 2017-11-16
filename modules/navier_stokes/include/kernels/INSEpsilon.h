#ifndef INSEpsilon_H
#define INSEpsilon_H

#include "INSScalar.h"

// Forward Declarations
class INSEpsilon;

template <>
InputParameters validParams<INSEpsilon>();

/**
 * This class computes the epsilon residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation with a standard k-epsilon turbulence model.
 */
class INSEpsilon : public INSScalar
{
public:
  INSEpsilon(const InputParameters & parameters);

  virtual ~INSEpsilon() {}

protected:
  virtual Real diffusivity() override;
  virtual RealVectorValue grad_diffusivity() override;
  virtual Real d_diff_d_uj() override;
  virtual RealVectorValue d_grad_diff_d_uj() override;

  virtual Real strongSink() override;
  virtual Real strongSinkJac() override;

  // Coupled variables
  const VariableValue & _kin;

  // Gradients
  const VariableGradient & _grad_kin;

  // Variable numberings
  unsigned int _kin_var_number;

  // Material properties
  const MaterialProperty<Real> & _diff_epsilon;
  const MaterialProperty<Real> & _d_diff_epsilon_d_epsilon;
  const MaterialProperty<RealVectorValue> & _grad_diff_epsilon;
  const MaterialProperty<RealVectorValue> & _d_grad_diff_epsilon_d_epsilon;
  const MaterialProperty<Real> & _d_grad_diff_epsilon_d_grad_epsilon;
  const MaterialProperty<Real> & _nu_turb;
  const MaterialProperty<Real> & _d_nu_turb_d_epsilon;

  // Parameters
  const Real _C1eps;
  const Real _C2eps;
};

#endif // INSEpsilon_H
