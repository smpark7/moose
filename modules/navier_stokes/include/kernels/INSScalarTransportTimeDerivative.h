#ifndef INSSCALARTRANSPORTTIMEDERIVATIVE_H
#define INSSCALARTRANSPORTTIMEDERIVATIVE_H

#include "TimeKernel.h"
#include "ScalarTransportBase.h"

// Forward Declaration
class INSScalarTransportTimeDerivative;

template <>
InputParameters validParams<INSScalarTransportTimeDerivative>();

class INSScalarTransportTimeDerivative : public TimeKernel, public ScalarTransportBase
{
public:
  INSScalarTransportTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  bool _lumping;

  const MaterialProperty<Real> & _rho;
};

#endif // INSSCALARTRANSPORTTIMEDERIVATIVE_H
