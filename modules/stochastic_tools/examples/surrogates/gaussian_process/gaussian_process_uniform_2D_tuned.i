[StochasticTools]
[]

[Distributions]
  [k_dist]
    type = Uniform
    lower_bound = 1
    upper_bound = 10
  []
  [q_dist]
    type = Uniform
    lower_bound = 9000
    upper_bound = 11000
  []
[]

[Samplers]
  [train_sample]
    type = MonteCarlo
    num_rows = 50
    distributions = 'k_dist q_dist'
    execute_on = PRE_MULTIAPP_SETUP
  []
  [cart_sample]
    type = CartesianProduct
    linear_space_items = '1 0.09 100
                          9000 20 100 '
    execute_on = initial
  []
[]

[MultiApps]
  [sub]
    type = SamplerFullSolveMultiApp
    input_files = sub.i
    sampler = train_sample
  []
[]

[Controls]
  [cmdline]
    type = MultiAppCommandLineControl
    multi_app = sub
    sampler = train_sample
    param_names = 'Materials/conductivity/prop_values Kernels/source/value'
  []
[]

[Transfers]
  [data]
    type = SamplerPostprocessorTransfer
    multi_app = sub
    sampler = train_sample
    to_vector_postprocessor = results
    from_postprocessor = 'avg'
  []
[]

[VectorPostprocessors]
  [results]
    type = StochasticResults
  []
[]

[Trainers]
  [GP_avg_trainer]
    type = GaussianProcessTrainer
    execute_on = timestep_end
    covariance_function = 'rbf'
    standardize_params = 'true'               #Center and scale the training params
    standardize_data = 'true'                 #Center and scale the training data
    tao_options = '-tao_bncg_type kd'
    distributions = 'k_dist q_dist'
    sampler = train_sample
    results_vpp = results
    results_vector = data:avg
    tune_parameters = ' signal_variance length_factor'
    tuning_min = ' 1e-9 1e-9'
    tuning_max = ' 1e16  1e16'
  []
[]


[Covariance]
  [rbf]
    type=SquaredExponentialCovariance
    signal_variance = 1                       #Use a signal variance of 1 in the kernel
    noise_variance = 1e-3                     #A small amount of noise can help with numerical stability
    length_factor = '0.38971 0.38971'         #Select a length factor for each parameter (k and q)
  []
[]

[Surrogates]
  [GP_avg]
    type = GaussianProcess
    trainer = 'GP_avg_trainer'
  []
[]

[VectorPostprocessors]
  [train_avg]
    type = EvaluateGaussianProcess
    model = GP_avg
    sampler = train_sample
    output_samples = true
    execute_on = final
  []
  [cart_avg]
    type = EvaluateGaussianProcess
    model = GP_avg
    sampler = cart_sample
    output_samples = true
    execute_on = final
  []
  [hyperparams]
    type = GaussianProcessData
    gp_name = 'GP_avg'
    execute_on = final
  []
[]

[Outputs]
  [out]
    type = CSV
    execute_on = FINAL
  []
[]
