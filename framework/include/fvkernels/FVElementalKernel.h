//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVKernel.h"
#include "MooseVariableFV.h"
#include "MooseVariableInterface.h"
#include "CoupleableMooseVariableDependencyIntermediateInterface.h"
#include "MaterialPropertyInterface.h"

class FVElementalKernel : public FVKernel,
                          public MooseVariableInterface<Real>,
                          public CoupleableMooseVariableDependencyIntermediateInterface,
                          public MaterialPropertyInterface
{
public:
  static InputParameters validParams();
  FVElementalKernel(const InputParameters & parameters);

  virtual void computeResidual();
  virtual void computeJacobian();
  virtual void computeOffDiagJacobian();
  virtual ADReal computeQpResidual() = 0;

protected:
  MooseVariableFV<Real> & _var;
  const ADVariableValue & _u;
  const unsigned int _qp = 0;
  const Elem * const & _current_elem;
};
