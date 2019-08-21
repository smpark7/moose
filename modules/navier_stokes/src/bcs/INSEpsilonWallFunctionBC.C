/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "INSEpsilonWallFunctionBC.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<INSEpsilonWallFunctionBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params += validParams<ScalarTransportBase>();

  params.addClassDescription("Epsilon BC for k-epsilon turbulence model.");
  params.addRequiredCoupledVar("kin", "The turbulent kinetic energy");

  // Optional parameters
  params.addParam<MaterialPropertyName>("mu_name", "mu", "The name of the dynamic viscosity");
  params.addParam<MaterialPropertyName>("rho_name", "rho", "The name of the density");

  return params;
}

INSEpsilonWallFunctionBC::INSEpsilonWallFunctionBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    ScalarTransportBase(parameters),

    // Coupled variables
    _kin(coupledValue("kin")),

    // Variable numberings
    _kin_var_number(coupled("kin")),

    // Material properties
    _mu(getMaterialProperty<Real>("mu_name")),
    _rho(getMaterialProperty<Real>("rho_name"))
{
}

Real
INSEpsilonWallFunctionBC::computeQpResidual()
{
  Real Cmu = 0.09;
  // Real yStarPlus = 11.06;
  Real sigmaEpsilon = 1.3;
  Real kappa = 0.41;

  // return -_test[_i][_qp] * std::pow(Cmu, 1.25) * std::pow(computeConcentration(_kin, _qp), 2.5) /
  //        (sigmaEpsilon * yStarPlus * _mu[_qp]);
  Real muT =
      _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) / _u[_qp];
  return -_test[_i][_qp] * (_mu[_qp] + muT / sigmaEpsilon) * kappa * std::pow(Cmu, 0.25) *
         std::sqrt(computeConcentration(_kin, _qp)) / muT * _u[_qp] * _rho[_qp];
}

Real
INSEpsilonWallFunctionBC::computeQpJacobian()
{
  Real Cmu = 0.09;
  // Real yStarPlus = 11.06;
  Real sigmaEpsilon = 1.3;
  Real kappa = 0.41;
  Real muT =
      _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) / _u[_qp];

  // return 0.;
  return -_test[_i][_qp] * (_mu[_qp] + muT / sigmaEpsilon) * kappa * std::pow(Cmu, 0.25) *
         std::sqrt(computeConcentration(_kin, _qp)) / muT * _phi[_j][_qp] * _rho[_qp];
}

Real
INSEpsilonWallFunctionBC::computeQpOffDiagJacobian(unsigned jvar)
{
  if (jvar == _kin_var_number)
  {
    Real Cmu = 0.09;
    // Real yStarPlus = 11.06;
    Real sigmaEpsilon = 1.3;
    Real kappa = 0.41;
    Real muT = _rho[_qp] * Cmu * computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp) /
               _u[_qp];

    // return -_test[_i][_qp] * std::pow(Cmu, 1.25) * 2.5 * std::pow(computeConcentration(_kin,
    // _qp), 1.5) * _phi[_j][_qp]
    // /
    //        (sigmaEpsilon * yStarPlus * _mu[_qp]);
    return -_test[_i][_qp] * (_mu[_qp] + muT / sigmaEpsilon) * kappa * std::pow(Cmu, 0.25) * 0.5 /
           std::sqrt(computeConcentration(_kin, _qp)) *
           computeConcentrationDerivative(_kin, _phi, _j, _qp) / muT * _u[_qp] * _rho[_qp];
  }

  else
    return 0.;
}
