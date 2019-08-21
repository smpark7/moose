# This input file tests the k-epsilon turbulence model

# Cmu = 0.9

[GlobalParams]
  gravity = '0 0 0'
  integrate_p_by_parts = false
  epsilon = epsilon
  kin = kin
  use_exp_form = true
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 10
  xmax = 50
  ymin = -15
  ymax = 15
  nx = 25
  ny = 25
  elem_type = QUAD9
  uniform_refine = 1
[]

# [MeshModifiers]
#   [./corner_node]
#     type = AddExtraNodeset
#     new_boundary = top_right
#     coord = '3 1'
#   [../]
# []


[Variables]
  [./vel_x]
    order = SECOND
    family = LAGRANGE
    # scaling = 1e2
    initial_condition = 1e-6
  [../]
  [./vel_y]
    order = SECOND
    family = LAGRANGE
    # scaling = 1e2
    initial_condition = 1e-6
  [../]
  [./p]
    order = FIRST
    family = LAGRANGE
    # scaling = 1e4
    initial_condition = 0
  [../]
  [./kin]
    initial_condition = 1e-6
    order = SECOND
    family = LAGRANGE
  [../]
  [./epsilon]
    initial_condition = 1e-6
    order = SECOND
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
    p = p
  [../]
  [./x_momentum_space]
    type = INSMomentumTractionForm
    variable = vel_x
    u = vel_x
    v = vel_y
    p = p
    component = 0
  [../]
  [./y_momentum_space]
    type = INSMomentumTractionForm
    variable = vel_y
    u = vel_x
    v = vel_y
    p = p
    component = 1
  [../]
  [./x_mom_turb_visc]
    type = INSMomentumTurbulentViscosityTractionForm
    variable = vel_x
    u = vel_x
    v = vel_y
    component = 0
  [../]
  [./y_mom_turb_visc]
    type = INSMomentumTurbulentViscosityTractionForm
    variable = vel_y
    u = vel_x
    v = vel_y
    component = 0
  [../]
  [./kin]
    variable = kin
    type = INSK
    u = vel_x
    v = vel_y
  [../]
  [./epsilon]
    variable = epsilon
    type = INSEpsilon
    u = vel_x
    v = vel_y
  [../]
  [./vel_x_time]
    type = INSMomentumTimeDerivative
    variable = vel_x
  [../]
  [./vel_y_time]
    type = INSMomentumTimeDerivative
    variable = vel_y
  [../]
  [./kin_time]
    type = INSScalarTransportTimeDerivative
    variable = kin
  [../]
  [./epsilon_time]
    type = INSMomentumTimeDerivative
    variable = epsilon
  [../]
  [./vel_x_supg]
    type = INSMomentumSUPG
    variable = vel_x
    u = vel_x
    v = vel_y
    component = 0
  [../]
  [./vel_y_supg]
    type = INSMomentumSUPG
    variable = vel_y
    u = vel_x
    v = vel_y
    component = 1
  [../]
  [./kin_supg]
    type = INSScalarSUPG
    variable = kin
    use_exp_form = true
    u = vel_x
    v = vel_y
  [../]
  [./eps_supg]
    type = INSScalarSUPG
    variable = epsilon
    use_exp_form = false
    u = vel_x
    v = vel_y
  [../]
[]

[BCs]
  [./diri_vel_x]
    type = FunctionDirichletBC
    variable = vel_x
    function = 'vel_x_func'
    boundary = 'left right top bottom'
  [../]
  [./diri_vel_y]
    type = FunctionDirichletBC
    variable = vel_y
    function = 'vel_y_func'
    boundary = 'left right top bottom'
  [../]
  [./diri_kin]
    type = FunctionDirichletBC
    variable = kin
    function = 'kin_func'
    boundary = 'left right top bottom'
  [../]
  [./diri_epsilon]
    type = FunctionDirichletBC
    variable = epsilon
    function = 'epsilon_func'
    boundary = 'left right top bottom'
  [../]
  [./diri_pressure]
    type = DirichletBC
    variable = p
    value = 0
    boundary = 'left right top bottom'
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'rho mu'
    prop_values = '1  1e-6'
  [../]
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = JFNK
    # ksp_norm = none
  [../]
[]

[Executioner]
  # type = Steady
  type = Transient
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor -options_left'
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
  # petsc_options_value = '300                asm	     ilu'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type -sub_pc_factor_levels'
  petsc_options_value = '300                bjacobi  ilu          4'
  # line_search = none
  nl_rel_tol = 1e-6
  nl_max_its = 100
  l_tol = 1e-6
  l_max_its = 300
  end_time = 10
  # num_steps = 30
  dt = .1
  dtmin = .05
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]

[Functions]
  [./vel_x_func]
    type = ParsedFunction
    value = 'max(1. / 2. - 1. / 2. * erf(13.5 * y / x), 1e-8)'
  [../]
  [./vel_y_func]
    type = ParsedFunction
    value = 'max(1. / 2. * 1. / (13.5 * sqrt(pi)) * exp(-(13.5 * y / x)^2), 1e-8)'
  [../]
  [./kin_func]
    type = ParsedFunction
    value = 'log(343. / 75000. * 13.5 / sqrt(pi) * (2e-3 / (343. / 75000. * 13.5 / sqrt(pi)) + exp(-(13.5 * y / x)^2)))'
  [../]
  [./epsilon_func]
    type = ParsedFunction
    value = '343. / 22500. * 0.9 * 13.5^2 / sqrt(pi) / x * (2e-3 / (343. / 75000. * 13.5 / sqrt(pi)) + exp(-(13.5 * y / x)^2))^2'
  [../]
  [./epsilon0_func]
    type = ParsedFunction
    value = '343. / 22500. * 0.9 * 13.5^2 / sqrt(pi)'
  [../]
  [./vel_x0_func]
    type = ParsedFunction
    value = 'if(y<0, 1, 1e-6)'
  [../]
[]

# [ICs]
#   [./eps_ic]
#     type = FunctionIC
#     variable = epsilon
#     function = epsilon_func
#   [../]
#   [./kin_ic]
#     type = FunctionIC
#     variable = kin
#     function = kin_func
#   [../]
#   [./vel_x_ic]
#     type = FunctionIC
#     variable = vel_x
#     function = vel_x_func
#   [../]
#   [./vel_y_ic]
#     type = FunctionIC
#     variable = vel_y
#     function = vel_y_func
#   [../]
# []

[Debug]
  show_var_residual_norms = true
[]

[Adaptivity]
  marker = errorfrac
  max_h_level = 2
  [./Indicators]
    [./error]
      type = GradientJumpIndicator
      variable = kin
      outputs = none
    [../]
  [../]
  [./Markers]
    [./errorfrac]
      type = ErrorFractionMarker
      refine = 0.6
      coarsen = 0.1
      indicator = error
      outputs = none
    [../]
  [../]
[]
