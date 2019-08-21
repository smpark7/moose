# This input file tests the k-epsilon turbulence model with a backwards step

[GlobalParams]
  gravity = '0 0 0'
  integrate_p_by_parts = false
  epsilon = epsilon
  kin = kin
  use_exp_form = true
[]

[Mesh]
  file = backwards_step.msh
  # uniform_refine = 1
[]

[Variables]
  [./vel_x]
    order = SECOND
    family = LAGRANGE
    initial_condition = 29.4
  [../]
  [./vel_y]
    order = SECOND
    family = LAGRANGE
  [../]
  [./p]
    order = FIRST
    family = LAGRANGE
  [../]
  [./kin]
    initial_condition = 2.157
    order = SECOND
    family = LAGRANGE
  [../]
  [./epsilon]
    initial_condition = 327
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
[]

[BCs]
  [./x_walls_slip]
    type = INSMomentumShearStressWallFunctionBC
    variable = vel_x
    boundary = 'norm_y_walls'
    u = vel_x
    v = vel_y
    p = p
    component = 0
  [../]
  [./x_walls_no_slip]
    type = DirichletBC
    variable = vel_x
    boundary = 'norm_x_walls'
    value = 0.0
  [../]
  [./x_inlet]
    type = FunctionDirichletBC
    variable = vel_x
    boundary = 'inlet'
    function = 'inlet_vel_x_func'
  [../]
  [./x_outlet]
    type = INSMomentumNoBCBCTurbulentTractionForm
    variable = vel_x
    u = vel_x
    v = vel_y
    p = p
    component = 0
    boundary = 'outlet'
  [../]
  # [./x_outlet]
  #   type = INSMomentumNoBCBCTractionForm
  #   variable = vel_x
  #   u = vel_x
  #   v = vel_y
  #   p = p
  #   component = 0
  #   boundary = 'right'
  # [../]
  [./y_walls_and_inlet_no_slip]
    type = DirichletBC
    variable = vel_y
    boundary = 'inlet norm_y_walls'
    value = 0.0
  [../]
  [./y_walls_slip]
    type = INSMomentumShearStressWallFunctionBC
    variable = vel_y
    boundary = 'norm_x_walls'
    u = vel_x
    v = vel_y
    p = p
    component = 1
  [../]
  [./y_outlet]
    type = INSMomentumNoBCBCTurbulentTractionForm
    variable = vel_y
    u = vel_x
    v = vel_y
    p = p
    component = 1
    boundary = 'outlet'
  [../]
  # [./y_outlet]
  #   type = INSMomentumNoBCBCTractionForm
  #   variable = vel_y
  #   u = vel_x
  #   v = vel_y
  #   p = p
  #   component = 1
  #   boundary = 'right'
  # [../]
  [./p_corner]
    # Since the pressure is not integrated by parts in this example,
    # it is only specified up to a constant by the natural outflow BC.
    # Therefore, we need to pin its value at a single location.
    type = DirichletBC
    boundary = corner
    value = 0
    variable = p
  [../]
  [./kin_inlet]
    type = FunctionDirichletBC
    variable = kin
    boundary = 'inlet'
    function = 'inlet_kin_func'
  [../]
  [./epsilon_walls]
    type = INSEpsilonWallFunctionBC
    boundary = 'norm_x_walls norm_y_walls'
    variable = epsilon
  [../]
  [./epsilon_inlet]
    type = FunctionDirichletBC
    variable = epsilon
    boundary = 'inlet'
    function = 'inlet_epsilon_func'
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    block = 'domain'
    prop_names = 'rho mu'
    prop_values = '1  1'
  [../]
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = NEWTON
    # ksp_norm = none
  [../]
[]

[Executioner]
  # type = Steady
  type = Transient
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor -options_left'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type -sub_pc_factor_levels -snes_linesearch_minlambda'
  petsc_options_value = '300                bjacobi  ilu          4			1e-3'
  # petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type -sub_ksp_type'
  # petsc_options_value = '300                asm  lu          preonly'
  # petsc_options_iname = '-pc_type -ksp_type -pc_factor_shift_type'
  # petsc_options_value = 'lu preonly positive_definite'
  # petsc_options_iname = '-pc_type -ksp_type -pc_factor_shift_amount'
  # petsc_options_value = 'lu preonly 1e-10'
  # line_search = none
  nl_rel_tol = 1e-6
  nl_max_its = 20
  l_tol = 1e-6
  l_max_its = 50
  end_time = 1
  dt = .001
  dtmin = .001
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    cutback_factor = 0.4
    growth_factor = 1.2
    optimal_iterations = 10
  [../]
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
[]

[Functions]
  # # Turbulent reynolds number
  # [./inlet_vel_x_func]
  #   type = ParsedFunction
  #   value = '-8000 * (y - 0.5)^2 + 2000'
  # [../]
  # [./inlet_kin_func]
  #   type = ParsedFunction
  #   value = '.01 * (-8000 * (y - 0.5)^2 + 2000)^2'
  # [../]
  # [./inlet_epsilon_func]
  #   type = ParsedFunction
  #   value = '${Cmu} / ${channel_width} * (.01 * (-8000 * (y - 0.5)^2 + 2000)^2)^(1.5)'
  # [../]
  # Laminar reynolds number
  # [./inlet_vel_x_func]
  #   type = ParsedFunction
  #   value = '-4 * (y - 0.5)^2 + 1'
  # [../]
  # [./inlet_kin_func]
  #   type = ParsedFunction
  #   value = '.01 * (-4 * (y - 0.5)^2 + 1)^2'
  # [../]
  # [./inlet_epsilon_func]
  #   type = ParsedFunction
  #   value = '${Cmu} / ${channel_width} * (.01 * (-4 * (y - 0.5)^2 + 1)^2)^(1.5)'
  # [../]
  # flat Cammi numbers
  [./inlet_vel_x_func]
    type = ParsedFunction
    value = '2.94e1'
  [../]
  [./inlet_kin_func]
    type = ParsedFunction
    value = '2.157'
  [../]
  [./inlet_epsilon_func]
    type = ParsedFunction
    value = '327'
  [../]
[]

[Adaptivity]
  marker = errorfrac
  max_h_level = 3
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
      refine = 0.5
      coarsen = 0.2
      indicator = error
      outputs = none
    [../]
  [../]
[]
