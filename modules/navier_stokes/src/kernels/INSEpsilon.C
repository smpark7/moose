#include "INSEpsilon.h"

template <>
InputParameters
validParams<INSEpsilon>()
{
  InputParameters params = validParams<Kernel>();
  params += validParams<ScalarTransportBase>();

  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", 0, "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", 0, "z-velocity"); // only required in 3D
  params.addRequiredCoupledVar("kin", "turbulent kinetic energy");

  // Optional parameters
  params.addParam<MaterialPropertyName>("mu_name", "mu", "dynamic viscosity");
  params.addParam<MaterialPropertyName>("rho_name", "rho", "density");
  return params;
}

INSEpsilon::INSEpsilon(const InputParameters & parameters)
  : Kernel(parameters),
    ScalarTransportBase(parameters),

    // Coupled variables
    _u_vel(coupledValue("u")),
    _v_vel(coupledValue("v")),
    _w_vel(coupledValue("w")),
    _kin(coupledValue("kin")),

    // Gradients
    _grad_u_vel(coupledGradient("u")),
    _grad_v_vel(coupledGradient("v")),
    _grad_w_vel(coupledGradient("w")),
    _grad_kin(coupledGradient("kin")),

    // Variable numberings
    _u_vel_var_number(coupled("u")),
    _v_vel_var_number(coupled("v")),
    _w_vel_var_number(coupled("w")),
    _kin_var_number(coupled("kin")),

    // Parameters
    _mu(getMaterialProperty<Real>("mu_name")),
    _rho(getMaterialProperty<Real>("rho_name")),
    _Cmu(0.09),
    _sigk(1.00),
    _sigeps(1.30),
    _C1eps(1.44),
    _C2eps(1.92)
{
}

Real
INSEpsilon::computeQpResidual()
{
  // The convection part
  Real convective_part = _rho[_qp] *
                         (_u_vel[_qp] * _grad_u[_qp](0) + _v_vel[_qp] * _grad_u[_qp](1) +
                          _w_vel[_qp] * _grad_u[_qp](2)) *
                         _test[_i][_qp];

  // The diffusive part
  Real eddy_visc = _rho[_qp] * _Cmu * computeConcentration(_kin, _qp) *
                   computeConcentration(_kin, _qp) / _u[_qp];
  Real diffusive_part = _grad_test[_i][_qp] * (_mu[_qp] + eddy_visc / _sigeps) * _grad_u[_qp];

  // Turbulent dissipative sink part
  Real dissipative_part =
      _test[_i][_qp] * (_C2eps * _rho[_qp] * _u[_qp] * _u[_qp] / computeConcentration(_kin, _qp));

  // Source term from velocity gradients
  Real strain_tensor_double_dot_product = 0.;
  for (int m = 0; m < 3; m++)
  {
    for (int n = 0; n < 3; n++)
    {
      Real Em;
      switch (m)
      {
        case 0:
          Em = _grad_u_vel[_qp](n);
          break;
        case 1:
          Em = _grad_v_vel[_qp](n);
          break;
        case 2:
          Em = _grad_w_vel[_qp](n);
          break;
        default:
          mooseError("Unrecognized index");
      }
      Real En;
      switch (n)
      {
        case 0:
          En = _grad_u_vel[_qp](m);
          break;
        case 1:
          En = _grad_v_vel[_qp](m);
          break;
        case 2:
          En = _grad_w_vel[_qp](m);
          break;
        default:
          mooseError("Unrecognized index");
      }
      Real Emn = (Em + En);
      strain_tensor_double_dot_product += Emn * Emn;
    }
  }

  Real source_part = -_test[_i][_qp] * _C1eps * 0.5 * eddy_visc * _u[_qp] /
                     computeConcentration(_kin, _qp) * strain_tensor_double_dot_product;

  return convective_part + diffusive_part + dissipative_part + source_part;
}

Real
INSEpsilon::computeQpJacobian()
{
  RealVectorValue U(_u_vel[_qp], _v_vel[_qp], _w_vel[_qp]);

  // Convective part
  Real convective_part = _rho[_qp] * U * _grad_phi[_j][_qp] * _test[_i][_qp];

  // The diffusive part
  Real eddy_visc = _rho[_qp] * _Cmu * computeConcentration(_kin, _qp) *
                   computeConcentration(_kin, _qp) / _u[_qp];
  Real d_eddy_visc_d_epsilon = -_rho[_qp] * _Cmu * computeConcentration(_kin, _qp) *
                               computeConcentration(_kin, _qp) / (_u[_qp] * _u[_qp]);
  Real diffusive_part =
      _grad_test[_i][_qp] * (_mu[_qp] + eddy_visc / _sigeps) * _grad_phi[_j][_qp] +
      _grad_test[_i][_qp] * d_eddy_visc_d_epsilon * _phi[_j][_qp] / _sigeps * _grad_u[_qp];

  // Turbulent dissipative sink part
  Real dissipative_part = _test[_i][_qp] * (_C2eps * _rho[_qp] * 2. * _u[_qp] * _phi[_j][_qp] /
                                            computeConcentration(_kin, _qp));

  // Source term from velocity gradients
  Real strain_tensor_double_dot_product = 0.;
  for (int m = 0; m < 3; m++)
  {
    for (int n = 0; n < 3; n++)
    {
      Real Em;
      switch (m)
      {
        case 0:
          Em = _grad_u_vel[_qp](n);
          break;
        case 1:
          Em = _grad_v_vel[_qp](n);
          break;
        case 2:
          Em = _grad_w_vel[_qp](n);
          break;
        default:
          mooseError("Unrecognized index");
      }
      Real En;
      switch (n)
      {
        case 0:
          En = _grad_u_vel[_qp](m);
          break;
        case 1:
          En = _grad_v_vel[_qp](m);
          break;
        case 2:
          En = _grad_w_vel[_qp](m);
          break;
        default:
          mooseError("Unrecognized index");
      }
      Real Emn = (Em + En);
      strain_tensor_double_dot_product += Emn * Emn;
    }
  }

  Real source_part = -_test[_i][_qp] * 0.5 * _C1eps *
                     (d_eddy_visc_d_epsilon * _phi[_j][_qp] * _u[_qp] + eddy_visc * _phi[_j][_qp]) /
                     computeConcentration(_kin, _qp) * strain_tensor_double_dot_product;

  return convective_part + diffusive_part + dissipative_part + source_part;
}

Real
INSEpsilon::computeQpOffDiagJacobian(unsigned jvar)
{
  // In Stokes/Laplacian version, off-diag Jacobian entries wrt u,v,w are zero
  if (jvar == _u_vel_var_number)
  {
    int p = 0;
    // Convective part
    Real convective_part = _phi[_j][_qp] * _grad_u[_qp](p) * _test[_i][_qp];

    // Source term from velocity gradients
    Real d_strain_tensor_double_dot_product_d_u_vel = 0.;
    for (int m = 0; m < 3; m++)
    {
      for (int n = 0; n < 3; n++)
      {
        if (p == m || p == n)
        {
          Real Em;
          switch (m)
          {
            case 0:
              Em = _grad_u_vel[_qp](n);
              break;
            case 1:
              Em = _grad_v_vel[_qp](n);
              break;
            case 2:
              Em = _grad_w_vel[_qp](n);
              break;
            default:
              mooseError("Unrecognized index");
          }
          Real En;
          switch (n)
          {
            case 0:
              En = _grad_u_vel[_qp](m);
              break;
            case 1:
              En = _grad_v_vel[_qp](m);
              break;
            case 2:
              En = _grad_w_vel[_qp](m);
              break;
            default:
              mooseError("Unrecognized index");
          }
          Real q1 = 0.;
          Real q2 = 0.;
          if (p == m)
            q1 = _grad_phi[_j][_qp](n);
          if (p == n)
            q2 = _grad_phi[_j][_qp](m);
          d_strain_tensor_double_dot_product_d_u_vel += (Em + En) * (q1 + q2);
        }
      }
    }

    Real eddy_visc = _rho[_qp] * _Cmu * computeConcentration(_kin, _qp) *
                     computeConcentration(_kin, _qp) / _u[_qp];
    Real source_part = -_test[_i][_qp] * 0.5 * _C1eps * eddy_visc * _u[_qp] /
                       computeConcentration(_kin, _qp) * d_strain_tensor_double_dot_product_d_u_vel;

    return convective_part + source_part;
  }

  else if (jvar == _v_vel_var_number)
  {
    int p = 1;
    // Convective part
    Real convective_part = _phi[_j][_qp] * _grad_u[_qp](p) * _test[_i][_qp];

    // Source term from velocity gradients
    Real d_strain_tensor_double_dot_product_d_v_vel = 0.;
    for (int m = 0; m < 3; m++)
    {
      for (int n = 0; n < 3; n++)
      {
        if (p == m || p == n)
        {
          Real Em;
          switch (m)
          {
            case 0:
              Em = _grad_u_vel[_qp](n);
              break;
            case 1:
              Em = _grad_v_vel[_qp](n);
              break;
            case 2:
              Em = _grad_w_vel[_qp](n);
              break;
            default:
              mooseError("Unrecognized index");
          }
          Real En;
          switch (n)
          {
            case 0:
              En = _grad_u_vel[_qp](m);
              break;
            case 1:
              En = _grad_v_vel[_qp](m);
              break;
            case 2:
              En = _grad_w_vel[_qp](m);
              break;
            default:
              mooseError("Unrecognized index");
          }
          Real q1 = 0.;
          Real q2 = 0.;
          if (p == m)
            q1 = _grad_phi[_j][_qp](n);
          if (p == n)
            q2 = _grad_phi[_j][_qp](m);
          d_strain_tensor_double_dot_product_d_v_vel += (Em + En) * (q1 + q2);
        }
      }
    }

    Real eddy_visc = _rho[_qp] * _Cmu * computeConcentration(_kin, _qp) *
                     computeConcentration(_kin, _qp) / _u[_qp];
    Real source_part = -_test[_i][_qp] * 0.5 * _C1eps * eddy_visc * _u[_qp] /
                       computeConcentration(_kin, _qp) * d_strain_tensor_double_dot_product_d_v_vel;

    return convective_part + source_part;
  }

  else if (jvar == _w_vel_var_number)
  {
    int p = 2;
    // Convective part
    Real convective_part = _phi[_j][_qp] * _grad_u[_qp](p) * _test[_i][_qp];

    // Source term from velocity gradients
    Real d_strain_tensor_double_dot_product_d_w_vel = 0.;
    for (int m = 0; m < 3; m++)
    {
      for (int n = 0; n < 3; n++)
      {
        if (p == m || p == n)
        {
          Real Em;
          switch (m)
          {
            case 0:
              Em = _grad_u_vel[_qp](n);
              break;
            case 1:
              Em = _grad_v_vel[_qp](n);
              break;
            case 2:
              Em = _grad_w_vel[_qp](n);
              break;
            default:
              mooseError("Unrecognized index");
          }
          Real En;
          switch (n)
          {
            case 0:
              En = _grad_u_vel[_qp](m);
              break;
            case 1:
              En = _grad_v_vel[_qp](m);
              break;
            case 2:
              En = _grad_w_vel[_qp](m);
              break;
            default:
              mooseError("Unrecognized index");
          }
          Real q1 = 0.;
          Real q2 = 0.;
          if (p == m)
            q1 = _grad_phi[_j][_qp](n);
          if (p == n)
            q2 = _grad_phi[_j][_qp](m);
          d_strain_tensor_double_dot_product_d_w_vel += (Em + En) * (q1 + q2);
        }
      }
    }

    Real eddy_visc = _rho[_qp] * _Cmu * computeConcentration(_kin, _qp) *
                     computeConcentration(_kin, _qp) / _u[_qp];
    Real source_part = -_test[_i][_qp] * 0.5 * _C1eps * eddy_visc * _u[_qp] /
                       computeConcentration(_kin, _qp) * d_strain_tensor_double_dot_product_d_w_vel;

    return convective_part + source_part;
  }

  else if (jvar == _kin_var_number)
  {
    // The diffusive part
    Real eddy_visc = _rho[_qp] * _Cmu * computeConcentration(_kin, _qp) *
                     computeConcentration(_kin, _qp) / _u[_qp];
    Real d_eddy_visc_d_kin = _rho[_qp] * _Cmu * 2. * computeConcentration(_kin, _qp) *
                             computeConcentrationDerivative(_kin, _phi, _j, _qp) / _u[_qp];
    Real diffusive_part = _grad_test[_i][_qp] * d_eddy_visc_d_kin / _sigeps * _grad_u[_qp];

    // Turbulent dissipative sink part
    Real dissipative_part =
        -_test[_i][_qp] * (_C2eps * _rho[_qp] * _u[_qp] * _u[_qp] /
                           (computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp)) *
                           computeConcentrationDerivative(_kin, _phi, _j, _qp));

    // Source term from velocity gradients
    Real strain_tensor_double_dot_product = 0.;
    for (int m = 0; m < 3; m++)
    {
      for (int n = 0; n < 3; n++)
      {
        Real Em;
        switch (m)
        {
          case 0:
            Em = _grad_u_vel[_qp](n);
            break;
          case 1:
            Em = _grad_v_vel[_qp](n);
            break;
          case 2:
            Em = _grad_w_vel[_qp](n);
            break;
          default:
            mooseError("Unrecognized index");
        }
        Real En;
        switch (n)
        {
          case 0:
            En = _grad_u_vel[_qp](m);
            break;
          case 1:
            En = _grad_v_vel[_qp](m);
            break;
          case 2:
            En = _grad_w_vel[_qp](m);
            break;
          default:
            mooseError("Unrecognized index");
        }
        Real Emn = (Em + En);
        strain_tensor_double_dot_product += Emn * Emn;
      }
    }

    Real source_part = -_test[_i][_qp] * 0.5 * _C1eps * _u[_qp] *
                       (d_eddy_visc_d_kin / computeConcentration(_kin, _qp) -
                        eddy_visc * computeConcentrationDerivative(_kin, _phi, _j, _qp) /
                            (computeConcentration(_kin, _qp) * computeConcentration(_kin, _qp))) *
                       strain_tensor_double_dot_product;

    return diffusive_part + dissipative_part + source_part;
  }

  else
    return 0;
}
