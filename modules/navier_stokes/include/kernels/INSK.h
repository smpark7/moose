#ifndef INSK_H
#define INSK_H

#include "INSScalar.h"

// Forward Declarations
class INSK;

template <>
InputParameters validParams<INSK>();

/**
 * This class computes turbulent energy residual and Jacobian
 * contributions for the incompressible Navier-Stokes momentum
 * equation with a standard k-epsilon turbulence model.
 */
class INSK : public INSScalar
{
public:
  INSK(const InputParameters & parameters);

  virtual ~INSK() {}

protected:
  virtual Real diffusivity() override;
  virtual RealVectorValue grad_diffusivity() override;
  virtual Real d_diff_d_uj() override;
  virtual RealVectorValue d_grad_diff_d_uj() override;

  virtual Real strongSink() override;
  virtual Real strongSinkJac() override;

  // Coupled variables
  const VariableValue & _epsilon;

  // Gradients
  const VariableGradient & _grad_epsilon;

  // Variable numberings
  unsigned int _epsilon_var_number;

  // Material properties
  const MaterialProperty<Real> & _diff_k;
  const MaterialProperty<Real> & _d_diff_k_d_k;
  const MaterialProperty<RealVectorValue> & _grad_diff_k;
  const MaterialProperty<RealVectorValue> & _d_grad_diff_k_d_k;
  const MaterialProperty<Real> & _d_grad_diff_k_d_grad_k;
  const MaterialProperty<Real> & _nu_turb;
  const MaterialProperty<Real> & _d_nu_turb_d_k;
};

#endif // INSK_H
