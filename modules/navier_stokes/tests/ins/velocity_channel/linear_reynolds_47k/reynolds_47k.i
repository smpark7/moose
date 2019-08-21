# This input file tests the k-epsilon turbulence model with a backwards step

[GlobalParams]
  gravity = '0 0 0'
  integrate_p_by_parts = false
  epsilon = epsilon
  kin = kin
  use_exp_form = false
[]

[Mesh]
  file = backwards_step.msh
  uniform_refine = 3
[]

[Variables]
  [./vel_x]
    order = SECOND
    family = LAGRANGE
    initial_condition = 100
    scaling = 1e-5
  [../]
  [./vel_y]
    order = SECOND
    family = LAGRANGE
    scaling = 1e-5
  [../]
  [./p]
    order = FIRST
    family = LAGRANGE
    scaling = 1e-1
  [../]
  [./kin]
    initial_condition = 4.605
    # order = SECOND
    # family = LAGRANGE
    scaling = 1e-4
  [../]
  [./epsilon]
    initial_condition = 12857
    # order = SECOND
    # family = LAGRANGE
    scaling = 1e-7
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
  [./vel_x_art_diff]
    type = INSScalarArtDiff
    variable = vel_x
    use_exp_form = false
    u = vel_x
    v = vel_y
  [../]
  [./vel_y_art_diff]
    type = INSScalarArtDiff
    variable = vel_y
    use_exp_form = false
    u = vel_x
    v = vel_y
  [../]
  [./kin_art_diff]
    type = INSScalarArtDiff
    variable = kin
    u = vel_x
    v = vel_y
  [../]
  [./epsilon_art_diff]
    type = INSScalarArtDiff
    variable = epsilon
    use_exp_form = false
    u = vel_x
    v = vel_y
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
    add_iso_art_diff = true
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
    add_iso_art_diff = true
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
    prop_values = '1  2.1e-3'
  [../]
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = PJFNK
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
  end_time = 10
  # dt = 1e-6
  dtmin = 1e-6
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-4
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
    value = '100'
  [../]
  [./inlet_kin_func]
    type = ParsedFunction
    value = '4.605'
  [../]
  [./inlet_epsilon_func]
    type = ParsedFunction
    value = '12857'
  [../]
[]

[Adaptivity]
  marker = combo
  max_h_level = 5
  [./Indicators]
    [./error_kin]
      type = GradientJumpIndicator
      variable = kin
      outputs = none
    [../]
    [./error_epsilon]
      type = GradientJumpIndicator
      variable = epsilon
      outputs = none
    [../]
    [./error_vel_x]
      type = GradientJumpIndicator
      variable = vel_x
      outputs = none
    [../]
    [./error_vel_y]
      type = GradientJumpIndicator
      variable = vel_y
      outputs = none
    [../]
  [../]
  [./Markers]
    [./errorfrac_kin]
      type = ErrorFractionMarker
      refine = 0.8
      coarsen = 0.1
      indicator = error_kin
      outputs = none
    [../]
    [./errorfrac_epsilon]
      type = ErrorFractionMarker
      refine = 0.8
      coarsen = 0.1
      indicator = error_epsilon
      outputs = none
    [../]
    [./errorfrac_vel_x]
      type = ErrorFractionMarker
      refine = 0.8
      coarsen = 0.1
      indicator = error_vel_x
      outputs = none
    [../]
    [./errorfrac_vel_y]
      type = ErrorFractionMarker
      refine = 0.8
      coarsen = 0.1
      indicator = error_vel_y
      outputs = none
    [../]
    [./combo]
      type = ComboMarker
      markers = 'errorfrac_kin errorfrac_epsilon errorfrac_vel_x errorfrac_vel_y'
    [../]
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]