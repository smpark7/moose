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
[]

[Kernels]
  [./vel_x]
    type = Diffusion
    variable = u
  [../]
  [./u_source]
    type = UserForcingFunction
    function = u_source_func
    variable = u
  [../]
[]

[BCs]
  [./u]
    type = FunctionDirichletBC
    function = u_func
    boundary = 'left right bottom'
    variable = u
  [../]
  [./u_vacuum]
    type = VacuumBC
    variable = u
    boundary = 'top'
  [../]
  [./u_fn_neumann]
    type = FunctionNeumannBC
    function = u_bc_source_func
    variable = u
    boundary = 'top'
  [../]
[]

[Functions]
  [./u_source_func]
    type = ParsedFunction
    value = '0.028*pi^2*x^2*sin(0.2*pi*x*y) + 0.028*pi^2*y^2*sin(0.2*pi*x*y) + 0.1*pi^2*sin(0.5*pi*x) + 0.4*pi^2*sin(pi*y)'
  [../]
  [./u_func]
    type = ParsedFunction
    value = '0.4*sin(0.5*pi*x) + 0.4*sin(pi*y) + 0.7*sin(0.2*pi*x*y) + 0.5'
  [../]
  [./u_bc_source_func]
    type = ParsedFunction
    value = '0.14*pi*x*cos(0.2*pi*x*y) + 0.2*sin(0.5*pi*x) + 0.2*sin(pi*y) + 0.35*sin(0.2*pi*x*y) + 0.4*pi*cos(pi*y) + 0.25'
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
  line_search = none
  nl_rel_tol = 1e-12
  nl_max_its = 6
  l_tol = 1e-6
  l_max_its = 300
[]

[Outputs]
  [./out]
    type = Exodus
  [../]
  [./csv]
    file_base = ''
    type = CSV
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
[]