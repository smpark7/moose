CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes         num_elem      
num_el_blk        num_node_sets      	   num_side_sets         num_el_in_blk1        num_nod_per_el1       num_side_ss1      num_side_ss2      num_side_ss3      num_side_ss4      num_side_ss5      num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_ns4       num_nod_ns5       num_nod_ns6       num_nod_ns7       num_nod_ns8       num_nod_ns9       num_nod_var       num_info  @   num_glo_var             api_version       @��H   version       @��H   floating_point_word_size            	file_size               title         	rz_out.e       maximum_name_length                 (   
time_whole                            t�   	eb_status                             
�   eb_prop1               name      ID              
�   	ns_status         	              $      
�   ns_prop1      	         name      ID        $      
�   	ss_status         
                       ss_prop1      
         name      ID              ,   coordx                      @      @   coordy                      @      �   eb_names                       $      �   ns_names      	                ,      �   ss_names      
                 �         
coor_names                         D      �   connect1                  	elem_type         QUAD4         P      �   elem_num_map                          L   elem_ss1                          `   side_ss1                          d   elem_ss2                          h   side_ss2                          x   elem_ss3                          �   side_ss3                          �   elem_ss4                          �   side_ss4                          �   elem_ss5                          �   side_ss5                          �   node_ns1                          �   node_ns2                          �   node_ns3                          �   node_ns4                          �   node_ns5                          �   node_ns6                          �   node_ns7                          �   node_ns8                          �   node_ns9                          �   vals_nod_var1                          @      t�   vals_nod_var2                          @      t�   name_nod_var                       D      �   info_records                      e@      $   name_glo_var                       $      td   vals_glo_var                             u                                    e      
      f      h      g                     
         ?�      ?��
=p��?��G�z�?��
=p�?�G�z�H?�      ?�\(�?��
=p��                ?��Q��?�z�G�{?�z�G�{?��Q��?�z�G�{?��Q��PATCH                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  disp_x                           disp_y                             ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                   -i                                                                               rz.i                                                                                                                                                              ### Input File ###                                                                                                                                                [BCs]                                                                                                                                                               [./101x]                                                                           type                         = FunctionDirichletBC                               boundary                     = 101                                               function                     = x101                                              time_periods                 =                                                   variable                     = disp_x                                          [../]                                                                                                                                                             [./101y]                                                                           type                         = FunctionDirichletBC                               boundary                     = 101                                               function                     = y101                                              time_periods                 =                                                   variable                     = disp_y                                          [../]                                                                                                                                                             [./102x]                                                                           type                         = FunctionDirichletBC                               boundary                     = 102                                               function                     = x102                                              time_periods                 =                                                   variable                     = disp_x                                          [../]                                                                                                                                                             [./102y]                                                                           type                         = FunctionDirichletBC                               boundary                     = 102                                               function                     = y102                                              time_periods                 =                                                   variable                     = disp_y                                          [../]                                                                                                                                                             [./103x]                                                                           type                         = FunctionDirichletBC                               boundary                     = 103                                               function                     = x103                                              time_periods                 =                                                   variable                     = disp_x                                          [../]                                                                                                                                                             [./103y]                                                                           type                         = FunctionDirichletBC                               boundary                     = 103                                               function                     = y103                                              time_periods                 =                                                   variable                     = disp_y                                          [../]                                                                                                                                                             [./104x]                                                                           type                         = FunctionDirichletBC                               boundary                     = 104                                               function                     = x104                                              time_periods                 =                                                   variable                     = disp_x                                          [../]                                                                                                                                                             [./104y]                                                                           type                         = FunctionDirichletBC                               boundary                     = 104                                               function                     = y104                                              time_periods                 =                                                   variable                     = disp_y                                          [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      l_abs_step_tol                 = -1                                              l_max_its                      = 20                                              l_tol                          = 1e-05                                           nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-10                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-08                                           no_fe_reinit                   = 0                                               petsc_options                  = '-snes_mf_operator -ksp_monitor'                petsc_options_iname            = '-pc_type -snes_type -snes_ls -ksp_gmres_r... estart'                                                                            petsc_options_value            = 'lu ls basic 101'                               scheme                         = backward-euler                                  type                           = Transient                                       _mesh                          = 0x10215d9e0                                     dt                             = 1                                               dtmax                          = 1e+30                                           dtmin                          = 0                                               end_time                       = 6                                               growth_factor                  = 2                                               n_startup_steps                = 0                                               num_steps                      = 6                                               predictor_scale                = 0                                               restart_file_base              =                                                 ss_check_tol                   = 1e-08                                           ss_tmin                        = 0                                               start_time                     = 0                                               sync_times                     = -1                                              time_dt                        =                                                 time_period_ends               =                                                 time_period_starts             =                                                 time_periods                   =                                                 time_t                         =                                                 trans_ss_check                 = 0                                             []                                                                                                                                                                [Functions]                                                                                                                                                         [./x101]                                                                           type                         = PiecewiseLinear                                   axis                         = 1                                                 scale_factor                 = 1                                                 x                            = '0 5 6'                                           y                            = '0 0 0.24'                                      [../]                                                                                                                                                             [./x102]                                                                           type                         = PiecewiseLinear                                   axis                         = 1                                                 scale_factor                 = 1                                                 x                            = '0 4 5'                                           y                            = '0 0 0.24'                                      [../]                                                                                                                                                             [./x103]                                                                           type                         = PiecewiseLinear                                   axis                         = 34690998                                          scale_factor                 = 1                                                 x                            = '0 4 5'                                           y                            = '0 0 0.24'                                      [../]                                                                                                                                                             [./x104]                                                                           type                         = PiecewiseLinear                                   axis                         = 34690998                                          scale_factor                 = 1                                                 x                            = '0 5 6'                                           y                            = '0 0 0.24'                                      [../]                                                                                                                                                             [./y101]                                                                           type                         = PiecewiseLinear                                   axis                         = 1                                                 scale_factor                 = 1                                                 x                            = '0 6'                                             y                            = '0 0'                                           [../]                                                                                                                                                             [./y102]                                                                           type                         = PiecewiseLinear                                   axis                         = 1                                                 scale_factor                 = 1                                                 x                            = '0 1 2 3'                                         y                            = '0 0 0.12 0'                                    [../]                                                                                                                                                             [./y103]                                                                           type                         = PiecewiseLinear                                   axis                         = 34690998                                          scale_factor                 = 1                                                 x                            = '0 1 3 4'                                         y                            = '0 0.12 0.12 0'                                 [../]                                                                                                                                                             [./y104]                                                                           type                         = PiecewiseLinear                                   axis                         = 34782669                                          scale_factor                 = 1                                                 x                            = '0 2 3 4'                                         y                            = '0 0 0.12 0'                                    [../]                                                                          []                                                                                                                                                                [Materials]                                                                                                                                                         [./density]                                                                        type                         = Density                                           block                        = PATCH                                             density                      = 12.3353                                           disp_r                       = disp_x                                            disp_x                       =                                                   disp_y                       =                                                   disp_z                       = disp_y                                          [../]                                                                                                                                                             [./stiffStuff1]                                                                    type                         = Elastic                                           active_crack_planes          =                                                   appended_property_name       =                                                   block                        = PATCH                                             bulk_modulus                 = 4.94066e-324                                      cracking_release             = abrupt                                            cracking_stress              = 0                                                 disp_r                       = disp_x                                            disp_x                       =                                                   disp_y                       =                                                   disp_z                       = disp_y                                            formulation                  =                                                   increment_calculation        = RashidApprox                                      lambda                       = 4.94066e-324                                      large_strain                 = 0                                                 max_cracks                   = 3                                                 poissons_ratio               = 0                                                 shear_modulus                = 4.94066e-324                                      temp                         =                                                   thermal_expansion            = 0                                                 youngs_modulus               = 1e+06                                           [../]                                                                          []                                                                                                                                                                [Mesh]                                                                             displacements                  = 'disp_x disp_y'                                 uniform_refine                 = 0                                               displacements                  = 'disp_x disp_y'                                 file                           = elastic_patch_rz.e                              ghosted_boundaries             =                                                 ghosted_boundaries_inflation   =                                                 nemesis                        = 0                                               patch_size                     = 40                                              skip_partitioning              = 0                                               type                           = MooseMesh                                       block_id                       =                                                 block_name                     =                                                 boundary_id                    =                                                 boundary_name                  =                                                 centroid_partitioner_direction =                                                 construct_side_list_from_node_list = 0                                           partitioner                    =                                                 second_order                   = 0                                               _dimension                     = 1                                             []                                                                                                                                                                [Output]                                                                           elemental_as_nodal             = 1                                               exodus                         = 1                                               exodus_inputfile_output        = 1                                               file_base                      = rz_out                                          gmv                            = 0                                               gnuplot_format                 = ps                                              interval                       = 1                                               iteration_plot_start_time      = 1.79769e+308                                    max_pps_rows_screen            = 0                                               nemesis                        = 0                                               num_restart_files              = 0                                               output_displaced               = 0                                               output_initial                 = 1                                               output_solution_history        = 0                                               output_variables               =                                                 perf_log                       = 1                                               postprocessor_csv              = 0                                               postprocessor_gnuplot          = 0                                               postprocessor_screen           = 1                                               print_linear_residuals         = 0                                               screen_interval                = 1                                               show_setup_log_early           = 0                                               tecplot                        = 0                                               tecplot_binary                 = 0                                               xda                            = 0                                             []                                                                                                                                                                [Postprocessors]                                                                                                                                                    [./mass]                                                                           type                         = Mass                                              block                        = ANY_BLOCK_ID                                      execute_on                   = timestep                                          output                       = both                                              variable                     = disp_x                                          [../]                                                                          []                                                                                                                                                                [Problem]                                                                          coord_type                     = RZ                                              fe_cache                       = 0                                             []                                                                                                                                                                [SolidMechanics]                                                                                                                                                    [./solid]                                                                          appended_property_name       =                                                   disp_r                       = disp_x                                            disp_x                       =                                                   disp_y                       =                                                   disp_z                       = disp_y                                            temp                         =                                                   use_displaced_mesh           = 1                                               [../]                                                                          []                                                                                                                                                                [Variables]                                                                                                                                                         [./disp_x]                                                                         block                        =                                                   family                       = LAGRANGE                                          initial_condition            = 0                                                 order                        = FIRST                                             scaling                      = 1                                                 initial_from_file_timestep   = 2                                                 initial_from_file_var        =                                                 [../]                                                                                                                                                             [./disp_y]                                                                         block                        =                                                   family                       = LAGRANGE                                          initial_condition            = 0                                                 order                        = FIRST                                             scaling                      = 1                                                 initial_from_file_timestep   = 2                                                 initial_from_file_var        =                                                 [../]                                                                          []                                                                               mass                                                                                                                                                                        @������?�                      ?k%�R��?c�@�W?v
"���        ?q���Y��                        ?��Ψ&�?aը���?�"d[⎇        ?�֦�.� ?��Q��@������@                       �ؙ��c?=f����\Y�        >� Hb��                ?��Q��?��uM�?���-H?���C3@        ?��3j,�?��Q��@      @                      �Q�;�I  <9�70  �d��ǀ         �T��@                 <��_�� ?��Q��?�z�G�{?�z�G�{?��Q��?�z�G�|?��Q��@������@                      �K�`��  ;��A   �p�'s          �rFw��                 �'��1   <K��   �A��+�  ���G�   ;���   <:�i   ;���   @������@              ?θQ��?�8f��H=?�+~��Z?��k�k        ?ĸX3��?θQ��        2`T`�@  ��;�N?9ǹ�o��F�Ɂ�6r=���  ?�-6r=���  @������@      ?θQ��?θQ��?ί`8k��?ΰ���k?Ϋm�%h�?θQ��?ά�k�?θQ��        ��N�U`  >���J�h���EKT>�Xe���P0���!0  ��:80���!0  @������