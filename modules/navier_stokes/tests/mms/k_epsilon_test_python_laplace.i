mu = 1.5
rho = 2.5

[GlobalParams]
  gravity = '0 0 0'
  integrate_p_by_parts = true
  use_exp_form = false
  u = u
  v = v
  p = p
  kin = kin
  epsilon = epsilon
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  elem_type = QUAD9
[]

[Variables]
  [./u]
    family = LAGRANGE
    order = SECOND
  [../]
  [./v]
    family = LAGRANGE
    order = SECOND
  [../]
  [./p]
  [../]
  [./kin]
    family = LAGRANGE
    order = SECOND
    initial_condition = 1e-6
  [../]
  [./epsilon]
    family = LAGRANGE
    order = SECOND
    initial_condition = 1e-6
  [../]
[]

[Kernels]
  [./vel_x]
    type = INSMomentumLaplaceForm
    variable = u
    component = 0
  [../]
  [./vel_y]
    type = INSMomentumLaplaceForm
    variable = v
    component = 1
  [../]
  [./vel_x_turb]
    type = INSMomentumTurbulentViscosityLaplaceForm
    variable = u
    component = 0
  [../]
  [./vel_y_turb]
    type = INSMomentumTurbulentViscosityLaplaceForm
    variable = v
    component = 1
  [../]
  [./pressure]
    type = INSMass
    variable = p
  [../]
  [./u_source]
    type = UserForcingFunction
    function = u_source_func
    variable = u
  [../]
  [./v_source]
    type = UserForcingFunction
    function = v_source_func
    variable = v
  [../]
  [./p_source]
    function = p_source_func
    variable = p
    type = UserForcingFunction
  [../]
  [./kin]
    type = INSK
    variable = kin
  [../]
  [./kin_source]
    function = k_source_func
    variable = kin
    type = UserForcingFunction
  [../]
  [./epsilon]
    type = INSEpsilon
    variable = epsilon
  [../]
  [./epsilon_source]
    function = eps_source_func
    variable = epsilon
    type = UserForcingFunction
  [../]
[]

[BCs]
  [./u]
    type = FunctionDirichletBC
    function = u_func
    boundary = 'left top right bottom'
    variable = u
  [../]
  [./v]
    type = FunctionDirichletBC
    function = v_func
    boundary = 'left top right bottom'
    variable = v
  [../]
  [./p]
    function = p_func
    variable = p
    type = FunctionDirichletBC
    boundary = 'left top right bottom'
  [../]
  [./k]
    function = k_func
    variable = kin
    type = FunctionDirichletBC
    boundary = 'left top right bottom'
  [../]
  [./eps]
    function = eps_func
    variable = epsilon
    type = FunctionDirichletBC
    boundary = 'left top right bottom'
  [../]
[]

[Functions]
  [./u_source_func]
    type = ParsedFunction
    value = ''
  [../]
  [./v_source_func]
    type = ParsedFunction
    value = ''
  [../]
  [./p_source_func]
    type = ParsedFunction
    value = ''
  [../]
  [./k_source_func]
    type = ParsedFunction
    value = ''
  [../]
  [./eps_source_func]
    type = ParsedFunction
    value = ''
  [../]
  [./u_func]
    type = ParsedFunction
    value = ''
  [../]
  [./v_func]
    type = ParsedFunction
    value = ''
  [../]
  [./p_func]
    type = ParsedFunction
    value = ''
  [../]
  [./k_func]
    type = ParsedFunction
    value = ''
  [../]
  [./eps_func]
    type = ParsedFunction
    value = ''
  [../]
[]

[Materials]
  [./const]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'rho mu'
    prop_values = '${rho}  ${mu}'
  [../]
[]

[Preconditioning]
  [./SMP_PJFNK]
    type = SMP
    full = true
    solve_type = PJFNK
  [../]
[]

[Executioner]
  type = Steady
  petsc_options = '-snes_converged_reason -ksp_converged_reason -snes_linesearch_monitor'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type -sub_pc_factor_levels'
  petsc_options_value = '300                bjacobi  ilu          4'
  # line_search = none
  nl_rel_tol = 1e-12
  nl_max_its = 10
  l_tol = 1e-6
  l_max_its = 50
[]

[Outputs]
  [./exodus]
    type = Exodus
    file_base = ''
  [../]
  [./csv]
    type = CSV
    file_base = ''
  [../]
[]

[Debug]
  show_var_residual_norms = true
[]

[Postprocessors]
  [./L2u]
    type = ElementL2Error
    variable = u
    function = u_func
    outputs = 'console csv'
  [../]
  [./L2v]
    variable = v
    function = v_func
    type = ElementL2Error
    outputs = 'console csv'
  [../]
  [./L2p]
    variable = p
    function = p_func
    type = ElementL2Error
    outputs = 'console csv'
  [../]
  [./L2kin]
    variable = kin
    function = k_func
    type = ElementL2Error
    outputs = 'console csv'
  [../]
  [./L2epsilon]
    variable = epsilon
    function = eps_func
    type = ElementL2Error
    outputs = 'console csv'
  [../]
[]
