CDF      
   "   
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes         num_elem      
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1        num_nod_per_el1       num_side_ss1      
num_df_ss1        num_side_ss2      
num_df_ss2        num_side_ss3      
num_df_ss3        num_side_ss4      
num_df_ss4        num_side_ss5      
num_df_ss5        num_side_ss6      
num_df_ss6        num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_ns4       num_nod_ns5       num_nod_ns6       num_nod_var    
   num_elem_var      num_info  U         api_version       @��H   version       @��H   floating_point_word_size            	file_size               title         out.e      maximum_name_length                 <   
time_whole                            �P   	eb_status                             �   eb_prop1               name      ID              �   	ns_status         	                    �   ns_prop1      	         name      ID              �   	ss_status         
                    �   ss_prop1      
         name      ID              �   coordx                      @      �   coordy                      @      (   coordz                      @      h   eb_names                       $      �   ns_names      	                 �      �   ss_names      
                 �      �   
coor_names                         d      \   connect1                  	elem_type         HEX8                 �   elem_num_map                          �   elem_ss1                          �   side_ss1                          �   dist_fact_ss1                              �   elem_ss2                             side_ss2                             dist_fact_ss2                                 elem_ss3                          4   side_ss3                          8   dist_fact_ss3                              <   elem_ss4                          \   side_ss4                          `   dist_fact_ss4                              d   elem_ss5                          �   side_ss5                          �   dist_fact_ss5                              �   elem_ss6                          �   side_ss6                          �   dist_fact_ss6                              �   node_ns1                          �   node_ns2                          �   node_ns3                          �   node_ns4                             node_ns5                             node_ns6                          $   vals_nod_var1                          @      �X   vals_nod_var2                          @      ��   vals_nod_var3                          @      ��   vals_nod_var4                          @      �   vals_nod_var5                          @      �X   vals_nod_var6                          @      ��   vals_nod_var7                          @      ��   vals_nod_var8                          @      �   vals_nod_var9                          @      �X   vals_nod_var10                         @      ��   name_nod_var                      L      4   name_elem_var                           �      �   vals_elem_var1eb1                                ��   vals_elem_var2eb1                                ��   vals_elem_var3eb1                                ��   vals_elem_var4eb1                                ��   vals_elem_var5eb1                                ��   vals_elem_var6eb1                                �    vals_elem_var7eb1                                �   info_records      !                k�      h                                                                              ��      ��      ��      ��      ?�      ?�      ?�      ?�      ��      ��      ?�      ?�      ��      ��      ?�      ?�      ?�      ��      ��      ?�      ?�      ��      ��      ?�                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     disp_x                           disp_y                           disp_z                           strain_xx                        stress_xx                        stress_yy                        stress_zz                        stress_xy                        stress_yz                        stress_zx                          strain_xx                        stress_xx                        stress_yy                        stress_zz                        stress_xy                        stress_yz                        stress_zx                         ####################                                                             # Created by MOOSE #                                                             ####################                                                             ### Command Line Arguments ###                                                   -i                                                                               cracking_test.i                                                                                                                                                   ### Input File ###                                                                                                                                                                                                                                 [Mesh]                                                                             second_order                   = 0                                               displacements                  = 'disp_x disp_y disp_z'                          file                           = cracking_test.e                                 nemesis                        = 0                                               patch_size                     = 40                                              skip_partitioning              = 0                                               uniform_refine                 = 0                                             []                                                                                                                                                                [Functions]                                                                        [./displ]                                                                          type                         = PiecewiseLinear                                   execute_on                   = residual                                          scale_factor                 = 1                                                 x                            = '0 1 2 3 4'                                       y                            = '0 1 0 -1 0'                                    [../]                                                                          []                                                                                                                                                                [Variables]                                                                        [./disp_z]                                                                         family                       = LAGRANGE                                          initial_condition            = 0                                                 order                        = FIRST                                             scaling                      = 1                                               [../]                                                                                                                                                             [./disp_y]                                                                         family                       = LAGRANGE                                          initial_condition            = 0                                                 order                        = FIRST                                             scaling                      = 1                                               [../]                                                                                                                                                             [./disp_x]                                                                         family                       = LAGRANGE                                          initial_condition            = 0                                                 order                        = FIRST                                             scaling                      = 1                                               [../]                                                                                                                                                             [./disp_z]                                                                         initial_from_file_timestep   = 2                                               [../]                                                                                                                                                             [./disp_y]                                                                         initial_from_file_timestep   = 2                                               [../]                                                                                                                                                             [./disp_x]                                                                         initial_from_file_timestep   = 2                                               [../]                                                                          []                                                                                                                                                                [AuxVariables]                                                                     [./strain_xx]                                                                      initial_from_file_timestep   = 2                                               [../]                                                                                                                                                             [./stress_xx]                                                                      initial_from_file_timestep   = 2                                               [../]                                                                                                                                                             [./stress_zx]                                                                      initial_from_file_timestep   = 2                                               [../]                                                                                                                                                             [./stress_yy]                                                                      initial_from_file_timestep   = 2                                               [../]                                                                                                                                                             [./stress_zz]                                                                      initial_from_file_timestep   = 2                                               [../]                                                                                                                                                             [./stress_xy]                                                                      initial_from_file_timestep   = 2                                               [../]                                                                                                                                                             [./stress_yz]                                                                      initial_from_file_timestep   = 2                                               [../]                                                                                                                                                             [./stress_zx]                                                                      family                       = MONOMIAL                                          initial_condition            = 0                                                 order                        = CONSTANT                                          scaling                      = 1                                               [../]                                                                                                                                                             [./stress_yz]                                                                      family                       = MONOMIAL                                          initial_condition            = 0                                                 order                        = CONSTANT                                          scaling                      = 1                                               [../]                                                                                                                                                             [./stress_xy]                                                                      family                       = MONOMIAL                                          initial_condition            = 0                                                 order                        = CONSTANT                                          scaling                      = 1                                               [../]                                                                                                                                                             [./stress_zz]                                                                      family                       = MONOMIAL                                          initial_condition            = 0                                                 order                        = CONSTANT                                          scaling                      = 1                                               [../]                                                                                                                                                             [./stress_yy]                                                                      family                       = MONOMIAL                                          initial_condition            = 0                                                 order                        = CONSTANT                                          scaling                      = 1                                               [../]                                                                                                                                                             [./stress_xx]                                                                      family                       = MONOMIAL                                          initial_condition            = 0                                                 order                        = CONSTANT                                          scaling                      = 1                                               [../]                                                                                                                                                             [./strain_xx]                                                                      family                       = MONOMIAL                                          initial_condition            = 0                                                 order                        = CONSTANT                                          scaling                      = 1                                               [../]                                                                          []                                                                                                                                                                [AuxKernels]                                                                       [./stress_zx]                                                                      type                         = MaterialTensorAux                                 execute_on                   = residual                                          index                        = 5                                                 quantity                     =                                                   tensor                       = stress                                            variable                     = stress_zx                                       [../]                                                                                                                                                             [./stress_yz]                                                                      type                         = MaterialTensorAux                                 execute_on                   = residual                                          index                        = 4                                                 quantity                     =                                                   tensor                       = stress                                            variable                     = stress_yz                                       [../]                                                                                                                                                             [./stress_xy]                                                                      type                         = MaterialTensorAux                                 execute_on                   = residual                                          index                        = 3                                                 quantity                     =                                                   tensor                       = stress                                            variable                     = stress_xy                                       [../]                                                                                                                                                             [./stress_zz]                                                                      type                         = MaterialTensorAux                                 execute_on                   = residual                                          index                        = 2                                                 quantity                     =                                                   tensor                       = stress                                            variable                     = stress_zz                                       [../]                                                                                                                                                             [./strain_xx]                                                                      type                         = MaterialTensorAux                                 execute_on                   = residual                                          index                        = 0                                                 quantity                     =                                                   tensor                       = total_strain                                      variable                     = strain_xx                                       [../]                                                                                                                                                             [./stress_yy]                                                                      type                         = MaterialTensorAux                                 execute_on                   = residual                                          index                        = 1                                                 quantity                     =                                                   tensor                       = stress                                            variable                     = stress_yy                                       [../]                                                                                                                                                             [./stress_xx]                                                                      type                         = MaterialTensorAux                                 execute_on                   = residual                                          index                        = 0                                                 quantity                     =                                                   tensor                       = stress                                            variable                     = stress_xx                                       [../]                                                                          []                                                                                                                                                                [BCs]                                                                              [./back]                                                                           type                         = PresetBC                                          boundary                     = 3                                                 execute_on                   = residual                                          value                        = 0                                                 variable                     = disp_z                                          [../]                                                                                                                                                             [./bottom]                                                                         type                         = PresetBC                                          boundary                     = 2                                                 execute_on                   = residual                                          value                        = 0                                                 variable                     = disp_y                                          [../]                                                                                                                                                             [./left]                                                                           type                         = PresetBC                                          boundary                     = 1                                                 execute_on                   = residual                                          value                        = 0                                                 variable                     = disp_x                                          [../]                                                                                                                                                             [./pull]                                                                           type                         = FunctionPresetBC                                  boundary                     = 4                                                 execute_on                   = residual                                          function                     = displ                                             variable                     = disp_x                                          [../]                                                                          []                                                                                                                                                                [Materials]                                                                        [./fred]                                                                           type                         = Elastic                                           block                        = 1                                                 cracking_stress              = 1.68e+06                                          disp_x                       = disp_x                                            disp_y                       = disp_y                                            disp_z                       = disp_z                                            execute_on                   = residual                                          increment_calculation        = RashidApprox                                      max_cracks                   = 3                                                 poissons_ratio               = 0                                                 thermal_expansion            = 0                                                 youngs_modulus               = 2.8e+07                                         [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      l_abs_step_tol                 = -1                                              l_max_its                      = 100                                             l_tol                          = 1e-05                                           nl_abs_step_tol                = 1e-50                                           nl_abs_tol                     = 1e-08                                           nl_max_funcs                   = 10000                                           nl_max_its                     = 100                                             nl_rel_step_tol                = 1e-50                                           nl_rel_tol                     = 1e-08                                           no_fe_reinit                   = 0                                               petsc_options                  = '-snes_mf_operator -ksp_monitor'                petsc_options_iname            = '-snes_type -snes_ls -ksp_gmres_restart -p... c_type -sub_pc_type'                                                               petsc_options_value            = 'ls basic 101 asm lu'                           scheme                         = backward-euler                                  type                           = Transient                                       dt                             = 0.025                                           dtmax                          = 1e+30                                           dtmin                          = 0                                               end_time                       = 0.1                                             execute_on                     = residual                                        n_startup_steps                = 0                                               num_steps                      = 1.79769e+308                                    ss_check_tol                   = 1e-08                                           ss_tmin                        = 0                                               start_time                     = 0                                               sync_times                     = -1                                              trans_ss_check                 = 0                                             []                                                                                                                                                                [Output]                                                                           exodus                         = 1                                               file_base                      = out                                             gmv                            = 0                                               gnuplot_format                 = ps                                              interval                       = 1                                               iteration_plot_start_time      = 1.79769e+308                                    nemesis                        = 0                                               output_displaced               = 0                                               output_initial                 = 1                                               perf_log                       = 1                                               postprocessor_csv              = 0                                               postprocessor_ensight          = 0                                               postprocessor_gnuplot          = 0                                               postprocessor_screen           = 1                                               print_linear_residuals         = 0                                               show_setup_log_early           = 0                                               tecplot                        = 0                                               tecplot_binary                 = 0                                               xda                            = 0                                             []                                                                                                                                                                [setup_quadrature]                                                                 order                          = AUTO                                            type                           = GAUSS                                         []                                                                                                                                                                [setup_dampers]                                                                  []                                                                                                                                                                [no_action]                                                                      []                                                                                                                                                                [init_problem]                                                                   []                                                                                                                                                                [SolidMechanics]                                                                   [./solid]                                                                          disp_r                       =                                                   disp_x                       = disp_x                                            disp_y                       = disp_y                                            disp_z                       = disp_z                                            temp                         =                                                 [../]                                                                          []                                                                                                                                                                [check_integrity]                                                                []                                                                                                                                                                [no_action]                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         ?�������                                ?�������?�������?�������?�������                                                                                                                                ?�C��Mj?�C��Mj?�C��Mj?�C��Mj?�C��Mj?�C��Mj?�C��Mj?�C��MjA%����A%����A%����A%����A%����A%����A%����A%�����
˱��u�
˱��u�
˱��u�
˱��u�
˱��u�
˱��u�
˱��u�
˱��u                                                                =���]OH=���]OH=���]OH=���]OH=���]OH=���]OH=���]OH=���]OH                                                                                                                                ?�C��MjA%�����
˱��u        =���]OH                ?�������                                ?�������?�������?�������?�������                                                                                                                                ?������?������?������?������?������?������?������?������A4�z��A4�z��A4�z��A4�z��A4�z��A4�z��A4�z��A4�z�藺!I�A���!I�A���!I�A���!I�A���!I�A���!I�A���!I�A���!I�A��                                                                =�f����}=�f����}=�f����}=�f����}=�f����}=�f����}=�f����}=�f����}                                                                                                                                ?������A4�z�藺!I�A��        =�f����}                ?�333334                                ?�333334?�333334?�333334?�333334                                                                                                                                ?��?X�r?��?X�r?��?X�r?��?X�r?��?X�r?��?X�r?��?X�r?��?X�r=��,�~=��,�~=��,�~=��,�~=��,�~=��,�~=��,�~=��,�~�0�I|C�0�I|C�0�I|C�0�I|C�0�I|C�0�I|C�0�I|C�0�I|C                                                                <��8�=�<��8�=�<��8�=�<��8�=�<��8�=�<��8�=�<��8�=�<��8�=�                                                                                                                                ?��?X�r=��,�~�0�I|C        <��8�=�                ?�������                                ?�������?�������?�������?�������                                                                                                                                ?�aU��A�?�aU��A�?�aU��A�?�aU��A�?�aU��A�?�aU��A�?�aU��A�?�aU��A�:��_����:��_����:��_����:��_����:��_����:��_����:��_����:��_���ϺAwEэF�AwEэF�AwEэF�AwEэF�AwEэF�AwEэF�AwEэF�AwEэF                                                                :]�J��:]�J��:]�J��:]�J��:]�J��:]�J��:]�J��:]�J��                                                                                                                                ?�aU��A�:��_���ϺAwEэF        :]�J��                