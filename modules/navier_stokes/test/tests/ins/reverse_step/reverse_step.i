[GlobalParams]
  gravity = '0 0 0'
  laplace = true
  transient_term = true
  supg = true
  pspg = true
  family = LAGRANGE
  order = FIRST
[]

[Mesh]
  file = 'reverse_step.e'
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = NEWTON
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1000

  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-15
  nl_max_its = 20

  trans_ss_check = true
  ss_check_tol = 1e-10
  line_search = 'none'
  [./TimeStepper]
    dt = .1
    type = IterationAdaptiveDT
    cutback_factor = 0.4
    growth_factor = 1.2
    optimal_iterations = 8
  [../]
[]

[Outputs]
  console = true
  [./out]
    type = Exodus
  [../]
[]

[Variables]
  [./vel_x]
  [../]
  [./vel_y]
  [../]
  [./p]
    order = FIRST
  [../]
[]

[BCs]
  [./vel_x_in]
    type = FunctionDirichletBC
    boundary = 'inlet'
    variable = vel_x
    function = '16 * (1 - y) * (y - .5)'
  [../]
  [./vel_y_in]
    type = DirichletBC
    boundary = 'inlet'
    variable = vel_y
    value = 0
  [../]
  [./vel_x_walls]
    type = DirichletBC
    boundary = 'walls'
    variable = vel_x
    value = 0
  [../]
  [./vel_y_walls]
    type = DirichletBC
    boundary = 'walls'
    variable = vel_y
    value = 0
  [../]
[]


[Kernels]
  [./x_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_x
  [../]
  [./y_momentum_time]
    type = INSMomentumTimeDerivative
    variable = vel_y
  [../]
  [./mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
    p = p
  [../]
  [./x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = vel_x
    v = vel_y
    p = p
    component = 0
  [../]
  [./y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_y
    u = vel_x
    v = vel_y
    p = p
    component = 1
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    block = 1
    prop_names = 'rho mu'
    prop_values = '1  1e-2'
  [../]
[]

[Postprocessors]
  [./flow_in]
    type = VolumetricFlowRate
    vel_x = vel_x
    vel_y = vel_y
    boundary = 'inlet'
    outputs = 'console'
    execute_on = 'timestep_end'
  [../]
  [./flow_out]
    type = VolumetricFlowRate
    vel_x = vel_x
    vel_y = vel_y
    boundary = 'outlet'
    outputs = 'console'
    execute_on = 'timestep_end'
  [../]
[]
