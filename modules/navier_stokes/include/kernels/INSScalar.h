#ifndef INSScalar_H
#define INSScalar_H

#include "INSBase.h"

// Forward Declarations
class INSScalar;

template <>
InputParameters validParams<INSScalar>();

class INSScalar : public INSBase
{
public:
  INSScalar(const InputParameters & parameters);

  virtual ~INSScalar() {}

protected:
  virtual Real diffusivity() override;
  virtual RealVectorValue grad_diffusivity();
  virtual Real d_diff_d_uj();
  virtual RealVectorValue d_grad_diff_d_uj();
  virtual Real d_tau_d_uj();

  virtual Real strongConvection();
  virtual Real strongConvectionJac();

  virtual Real strongDiffusion();
  virtual Real strongDiffusionJac();
  virtual Real weakDiffusion();
  virtual Real weakDiffusionJac();

  virtual Real strongSink();
  virtual Real strongSinkJac();

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  virtual Real computeQpPGResidual();
  virtual Real computeQpPGJacobian();
  virtual Real computeQpPGOffDiagJacobian(unsigned comp);

  const VariableSecond & _second_u;

  const bool _supg;
};

#endif // INSScalar_H
