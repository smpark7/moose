mu = 1.5
rho = 2.5

[GlobalParams]
  gravity = '0 0 0'
  integrate_p_by_parts = true
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 5
  ny = 5
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
    family = LAGRANGE
    order = FIRST
  [../]
[]

[Kernels]
  [./vel_x]
    type = INSMomentumTractionForm
    variable = u
    component = 0
    u = u
    v = v
    p = p
  [../]
  [./vel_y]
    type = INSMomentumTractionForm
    variable = v
    component = 1
    u = u
    v = v
    p = p
  [../]
  [./pressure]
    type = INSMass
    variable = p
    u = u
    v = v
    p = p
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
    type = UserForcingFunction
    function = p_source_func
    variable = p
  [../]
[]

[BCs]
  [./u]
    type = FunctionDirichletBC
    function = u_func
    boundary = 'left right top bottom'
    variable = u
  [../]
  [./v]
    type = FunctionDirichletBC
    function = v_func
    boundary = 'left right top bottom'
    variable = v
  [../]
  [./p]
    type = FunctionDirichletBC
    function = p_func
    boundary = 'left top right bottom'
    variable = p
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
    solve_type = NEWTON
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
  l_max_its = 300
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
[]