module HetaSimulator

  using DiffEqBase: isempty
  using Reexport
  # heta compiler support
  using NodeJS
  # diffeq-related pkgs
  using LabelledArrays
  using DataStructures
  using DiffEqBase
  @reexport using OrdinaryDiffEq
  @reexport using DiffEqBase.EnsembleAnalysis
  using Sundials
  @reexport using NaNMath

  # fitting
  using NLopt

  # utils
  @reexport using DataFrames
  @reexport using Distributions
  using LinearAlgebra: pinv, diag
  using Distributed
  using RecursiveArrayTools: vecvec_to_mat, VectorOfArray
  using ProgressMeter
  #ProgressMeter.ijulia_behavior(:clear)
  # measurements 
  using CSV
  using XLSX
  #plots
  using RecipesBase

  const HetaSimulatorDir = dirname(Base.@__DIR__)

  include("types.jl")
  include("heta_connect.jl")
  include("utils.jl")
  include("events.jl")
  include("measurements.jl")
  include("ode_problem.jl")
  include("scenario.jl")
  include("parameters.jl")
  include("simulate.jl")
  include("saving.jl")
  include("solution_interface.jl")
  include("plots.jl")
  include("loss.jl")
  include("fit.jl")
  include("estimator.jl")
  include("monte_carlo.jl")
  include("import_platform.jl")
  include("gsa.jl")
  include("save_as_heta.jl")

  
  heta_update()

  export heta_update, heta_update_dev, heta_build, load_platform, load_jlplatform, load_jlmodel
  export Platform, Model, Scenario, Params
  export read_scenarios, add_scenarios!
  export read_measurements, add_measurements!, measurements_as_table
  export read_parameters
  export models, scenarios, scenario, constants, records, switchers, events, parameters, events_active, events_save, observables  # variables, dynamic, static
  export measurements, tspan 
  export CVODE_BDF, CVODE_Adams
  export optim, obj
  export sim, mc, mc!
  export fit, loss, estimator
  export HetaSimulatorDir
  export update
  export times, vals, status, status_summary
  export save_results, read_mcvecs
  export gsa, pearson, partial, standard
  export save_as_heta
  export scale_params, unscale_params
end
