/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef INSTURBULENCEMATERIAL_H
#define INSTURBULENCEMATERIAL_H

#include "Material.h"

class INSTurbulenceMaterial;

template <>
InputParameters validParams<INSTurbulenceMaterial>();

class INSTurbulenceMaterial : public Material
{
public:
  INSTurbulenceMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const VariableValue & _k;
  const VariableValue & _epsilon;
  const VariableGradient & _grad_k;
  const VariableGradient & _grad_epsilon;

  const Real _user_mu;
  const Real _user_rho;
  const Real _Cmu;
  const Real _sigk;
  const Real _sigeps;
  const Real _C1eps;
  const Real _C2eps;

  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _rho;
  MaterialProperty<Real> & _nu;

  MaterialProperty<Real> & _nu_turb;
  MaterialProperty<Real> & _d_nu_turb_d_k;
  MaterialProperty<Real> & _d_nu_turb_d_epsilon;
  MaterialProperty<RealVectorValue> & _grad_nu_turb;
  MaterialProperty<RealVectorValue> & _d_grad_nu_turb_d_k;
  MaterialProperty<Real> & _d_grad_nu_turb_d_grad_k;
  MaterialProperty<RealVectorValue> & _d_grad_nu_turb_d_epsilon;
  MaterialProperty<Real> & _d_grad_nu_turb_d_grad_epsilon;

  MaterialProperty<Real> & _diff_k;
  MaterialProperty<Real> & _d_diff_k_d_k;
  MaterialProperty<Real> & _d_diff_k_d_epsilon;
  MaterialProperty<RealVectorValue> & _grad_diff_k;
  MaterialProperty<RealVectorValue> & _d_grad_diff_k_d_k;
  MaterialProperty<Real> & _d_grad_diff_k_d_grad_k;
  MaterialProperty<RealVectorValue> & _d_grad_diff_k_d_epsilon;
  MaterialProperty<Real> & _d_grad_diff_k_d_grad_epsilon;

  MaterialProperty<Real> & _diff_epsilon;
  MaterialProperty<Real> & _d_diff_epsilon_d_k;
  MaterialProperty<Real> & _d_diff_epsilon_d_epsilon;
  MaterialProperty<RealVectorValue> & _grad_diff_epsilon;
  MaterialProperty<RealVectorValue> & _d_grad_diff_epsilon_d_k;
  MaterialProperty<Real> & _d_grad_diff_epsilon_d_grad_k;
  MaterialProperty<RealVectorValue> & _d_grad_diff_epsilon_d_epsilon;
  MaterialProperty<Real> & _d_grad_diff_epsilon_d_grad_epsilon;
};

#endif // INSTURBULENCEMATERIAL_H
