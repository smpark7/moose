/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSMOMENTUMNOBCBCTURBULENTTRACTIONFORM_H
#define INSMOMENTUMNOBCBCTURBULENTTRACTIONFORM_H

#include "INSMomentumNoBCBCTractionForm.h"
#include "ScalarTransportBase.h"

// Forward Declarations
class INSMomentumNoBCBCTurbulentTractionForm;

template <>
InputParameters validParams<INSMomentumNoBCBCTurbulentTractionForm>();

/**
 * This class implements the "No BC" boundary condition for the k-epsilon turbulence model based on
 * the
 * "traction" form of the viscous stress tensor.
 */
class INSMomentumNoBCBCTurbulentTractionForm : public INSMomentumNoBCBCTractionForm,
                                               public ScalarTransportBase
{
public:
  INSMomentumNoBCBCTurbulentTractionForm(const InputParameters & parameters);

  virtual ~INSMomentumNoBCBCTurbulentTractionForm() {}

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned jvar) override;

  const VariableValue & _kin;
  const VariableValue & _epsilon;

  unsigned _kin_var_number;
  unsigned _epsilon_var_number;
};

#endif
