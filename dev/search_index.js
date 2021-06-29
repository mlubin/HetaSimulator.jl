var documenterSearchIndex = {"docs":
[{"location":"api/#API-references","page":"API","title":"API references","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [HetaSimulator]\nOrder   = [:type, :function]","category":"page"},{"location":"api/#HetaSimulator.fit-Tuple{Platform, Vector{Pair{Symbol, Float64}}}","page":"API","title":"HetaSimulator.fit","text":"fit(platform::Platform,\n  params::Vector{Pair{Symbol,Float64}};\n  conditions::Union{AbstractVector{Symbol}, Nothing} = nothing,\n  kwargs...\n) where C<:AbstractCond\n\nFit parameters to experimental measurements. Returns FitResults type. Example: fit(platform, [:k1=>0.1,:k2=>0.2,:k3=>0.3];conditions=[:cond2,:cond3])\n\nArguments:\n\nplatform : platform of Platform type\nparams : optimization parameters and their initial values\nconditions : vector of conditions of type Cond or nothing to fit all conditions. Default is nothing\nkwargs : other solver related arguments supported by fit(condition_pairs::Vector{<:Pair}, params::Vector{<:Pair}\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.fit-Union{Tuple{C}, Tuple{AbstractArray{Pair{Symbol, C}, 1}, Vector{Pair{Symbol, Float64}}}} where C<:HetaSimulator.AbstractCond","page":"API","title":"HetaSimulator.fit","text":"fit(condition_pairs::AbstractVector{Pair{Symbol, C}},\n  params::Vector{Pair{Symbol,Float64}};\n  alg=DEFAULT_ALG,\n  reltol=DEFAULT_FITTING_RELTOL,\n  abstol=DEFAULT_FITTING_ABSTOL,\n  parallel_type=EnsembleSerial(),\n  ftol_abs = 0.0,\n  ftol_rel = 1e-4, \n  xtol_rel = 0.0,\n  xtol_abs = 0.0, \n  fit_alg = :LN_NELDERMEAD,\n  maxeval = 10000,\n  maxtime = 0.0,\n  lbounds = fill(0.0, length(params)),\n  ubounds = fill(Inf, length(params)),\n  kwargs... \n) where C<:AbstractCond\n\nFit parameters to experimental measurements. Returns FitResults type. Example: fit([:x=>cond2, :y=>cond3, :z=>cond4], [:k1=>0.1,:k2=>0.2,:k3=>0.3])\n\nArguments:\n\ncondition_pairs : vector of pairs containing names and conditions of type Cond\nparams : optimization parameters and their initial values\nalg : ODE solver. See SciML docs for details. Default is AutoTsit5(Rosenbrock23())\nreltol : relative tolerance. Default is 1e-6\nabstol : relative tolerance. Default is 1e-8\nparallel_type : parallel setup. See SciML docs for details. Default is no parallelism: EnsembleSerial()\nftol_abs : absolute tolerance on function value. See NLopt.jl docs for details. Default is 0.0\nftol_rel : relative tolerance on function value. See NLopt.jl docs for details. Default is 1e-4\nxtol_rel : relative tolerance on optimization parameters. See NLopt.jl docs for details. Default is 0.0\nxtol_rel : absolute tolerance on optimization parameters. See NLopt.jl docs for details. Default is 0.0\nfit_alg : fitting algorithm. See NLopt.jl docs for details. Default is :LN_NELDERMEAD\nmaxeval : maximum number of function evaluations. See NLopt.jl docs for details. Default is 1e4\nmaxtime : maximum optimization time (in seconds). See NLopt.jl docs for details. Default is 0\nlbounds : lower parameters bounds. See NLopt.jl docs for details. Default is fill(0.0, length(params))\nubounds : upper parameters bounds. See NLopt.jl docs for details. Default is fill(Inf, length(params))\nkwargs : other solver related arguments supported by DiffEqBase.solve. See SciML docs for details\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.fit-Union{Tuple{C}, Tuple{AbstractVector{C}, Vector{Pair{Symbol, Float64}}}} where C<:HetaSimulator.AbstractCond","page":"API","title":"HetaSimulator.fit","text":"fit(conditions::AbstractVector{C},\n  params::Vector{Pair{Symbol,Float64}};\n  kwargs...\n) where C<:AbstractCond\n\nFit parameters to experimental measurements. Returns FitResults type. Example: fit([cond2, cond3, cond4], [:k1=>0.1,:k2=>0.2,:k3=>0.3])\n\nArguments:\n\nconditions : vector of conditions of type Cond\nparams : optimization parameters and their initial values\nkwargs : other solver related arguments supported by fit(condition_pairs::Vector{<:Pair}, params::Vector{<:Pair}\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.heta_build-Tuple{AbstractString}","page":"API","title":"HetaSimulator.heta_build","text":"heta_build(  \n  heta_index::AbstractString;\n  declaration::String = \"platform\",\n  skip_export::Bool = false,\n  log_mode::String = \"error\",\n  debug::Bool = false,\n  julia_only::Bool = false,\n  dist_dir::String = \"dist\",\n  meta_dir::String = \"meta\",\n  source::String = \"index.heta\",\n  type::String = \"heta\"\n)\n\nBuilds the model from Heta-based reactions\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.heta_update-Tuple{}","page":"API","title":"HetaSimulator.heta_update","text":"heta_update(version::String)\n\nInstalls heta-compiler from NPM.\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.load_jlplatform-Tuple{AbstractString}","page":"API","title":"HetaSimulator.load_jlplatform","text":"load_jlplatform(  \n  model_jl::AbstractString; \n  rm_out::Bool = false\n)\n\nLoads prebuild julia model as part of platform\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.load_platform-Tuple{AbstractString}","page":"API","title":"HetaSimulator.load_platform","text":"load_platform(  \n  heta_index::AbstractString;\n  rm_out::Bool = true,\n  julia_only::Bool = true, \n  dist_dir::String = \".\",\n  source::String = \"index.heta\",\n  type::String = \"heta\",\n  kwargs...\n)\n\nConverts heta model to Julia and outputs platform type\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.mc-Tuple{Cond, DataFrame}","page":"API","title":"HetaSimulator.mc","text":"mc(cond::Cond,\n  params::DataFrame,\n  num_iter::Int64;\n  kwargs...\n)\n\nRun Monte-Carlo simulations with single condition cond. Returns MCResults type. Example: mc(cond1, DataFrame(k2=rand(3),k3=rand(3)), 1000)\n\nArguments:\n\ncond : simulation condition of type Cond\nparams : DataFrame with pre-generated parameters.\nnum_iter : number of Monte-Carlo iterations \nkwargs : other solver related arguments supported by mc(cond::Cond, params::Vector, num_iter::Int64)\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.mc-Union{Tuple{PP}, Tuple{CP}, Tuple{Vector{CP}, Vector{PP}, Int64}} where {CP<:Pair, PP<:Pair}","page":"API","title":"HetaSimulator.mc","text":"mc(cond_pairs::Vector{<:Pair},\n  params::Vector{<:Pair},\n  num_iter::Int64;\n  verbose=false,\n  alg=DEFAULT_ALG,\n  reltol=DEFAULT_SIMULATION_RELTOL,\n  abstol=DEFAULT_SIMULATION_ABSTOL,\n  parallel_type=EnsembleSerial(),\n  kwargs...\n)\n\nRun Monte-Carlo simulations with single condition cond. Returns Vector{MCResults} type. Example: mc([:c1=>cond1,:c2=>cond2], [:k2=>Normal(1e-3,1e-4), :k3=>Uniform(1e-4,1e-2)], 1000)\n\nArguments:\n\ncond_pairs : vector of pairs containing names and conditions of type Cond\nparams : parameters variation setup\nnum_iter : number of Monte-Carlo iterations\nverbose : print iteration progress. Default is false\nalg : ODE solver. See SciML docs for details. Default is AutoTsit5(Rosenbrock23())\nreltol : relative tolerance. Default is 1e-3\nabstol : relative tolerance. Default is 1e-6\nparallel_type : parallel setup. See SciML docs for details. Default is no parallelism: EnsembleSerial()\nkwargs : other solver related arguments supported by DiffEqBase.solve. See SciML docs for details\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.mc-Union{Tuple{P}, Tuple{Cond, Vector{P}, Int64}} where P<:Pair","page":"API","title":"HetaSimulator.mc","text":"mc(cond::Cond,\n  params::Vector{<:Pair},\n  num_iter::Int64;\n  verbose=false,\n  alg=DEFAULT_ALG,\n  reltol=DEFAULT_SIMULATION_RELTOL,\n  abstol=DEFAULT_SIMULATION_ABSTOL,\n  parallel_type=EnsembleSerial(),\n  kwargs...\n)\n\nRun Monte-Carlo simulations with single condition cond. Returns MCResults type. Example: mc(cond, [:k2=>Normal(1e-3,1e-4), :k3=>Uniform(1e-4,1e-2)], 1000)\n\nArguments:\n\ncond : simulation condition of type Cond\nparams : parameters variation setup\nnum_iter : number of Monte-Carlo iterations\nverbose : print iteration progress. Default is false\nalg : ODE solver. See SciML docs for details. Default is AutoTsit5(Rosenbrock23())\nreltol : relative tolerance. Default is 1e-3\nabstol : relative tolerance. Default is 1e-6\nparallel_type : parallel setup. See SciML docs for details. Default is no parallelism: EnsembleSerial()\nkwargs : other solver related arguments supported by DiffEqBase.solve. See SciML docs for details\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.mc-Union{Tuple{P}, Tuple{C}, Tuple{Vector{C}, Vector{P}, Int64}} where {C<:HetaSimulator.AbstractCond, P<:Pair}","page":"API","title":"HetaSimulator.mc","text":"mc(conds::Vector{<:AbstractCond},\n  params::Vector{<:Pair},\n  num_iter::Int64;\n  kwargs...\n)\n\nRun Monte-Carlo simulations with single condition cond. Returns Vector{MCResults} type. Example: mc([cond1,cond2], [:k2=>Normal(1e-3,1e-4), :k3=>Uniform(1e-4,1e-2)], 1000)\n\nArguments:\n\ncond_pairs : vector of conditions of type Cond\nparams : parameters variation setup\nnum_iter : number of Monte-Carlo iterations\nkwargs : other solver related arguments supported by mc(cond_pairs::Vector{<:Pair}, params::Vector, num_iter::Int64)\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.mc-Union{Tuple{P}, Tuple{Model, Vector{P}, Int64}} where P<:Pair","page":"API","title":"HetaSimulator.mc","text":"mc(model::Model,\n  params::Vector{<:Pair},\n  num_iter::Int64;\n  measurements::Vector{AbstractMeasurementPoint} = AbstractMeasurementPoint[],\n  events_active::Union{Nothing, Vector{Pair{Symbol,Bool}}} = Pair{Symbol,Bool}[],\n  events_save::Union{Tuple,Vector{Pair{Symbol, Tuple{Bool, Bool}}}}=(true,true), \n  observables::Union{Nothing,Vector{Symbol}} = nothing,\n  saveat::Union{Nothing,AbstractVector} = nothing,\n  tspan::Union{Nothing,Tuple} = nothing,\n  save_scope::Bool=false,\n  time_type::DataType=Float64,\n  kwargs...\n)\n\nRun Monte-Carlo simulations with Model. Returns MCResults type. Example: mc(model, [:k2=>Normal(1e-3,1e-4), :k3=>Uniform(1e-4,1e-2)], 1000)\n\nArguments:\n\nmodel : model of type Model\nparams : parameters variation setup\nnum_iter : number of Monte-Carlo iterations\nmeasurements : Vector of measurements. Default is empty vector \nevents_active : Vector of Pairs containing events' names and true/false values. Overwrites default model's values. Default is empty vector \nevents_save : Tuple or Vector{Tuple} marking whether to save solution before and after event. Default is (true,true) for all events\nobservables : names of output observables. Overwrites default model's values. Default is empty vector\nsaveat : time points, where solution should be saved. Default nothing values stands for saving solution at timepoints reached by the solver \ntspan : time span for the ODE problem\nsave_scope : should scope be saved together with solution. Default is false\nkwargs : other solver related arguments supported by mc(cond::Cond, params::Vector, num_iter::Int64)\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.mc-Union{Tuple{P}, Tuple{Platform, Vector{P}, Int64}} where P<:Pair","page":"API","title":"HetaSimulator.mc","text":"mc(platform::Platform, \n  params::Vector{<:Pair},\n  num_iter::Int64;\n  kwargs...\n)\n\nRun Monte-Carlo simulations with single condition cond. Returns Vector{MCResults} type. Example: mc(platform, [:k2=>Normal(1e-3,1e-4), :k3=>Uniform(1e-4,1e-2)], 1000)\n\nArguments:\n\nplatform : platform of Platform type\nparams : parameters variation setup\nnum_iter : number of Monte-Carlo iterations\nkwargs : other solver related arguments supported by mc(cond_pairs::Vector{<:Pair}, params::Vector, num_iter::Int64)\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.sim-Tuple{Model}","page":"API","title":"HetaSimulator.sim","text":"sim(model::Model; \n  parameters::Vector{Pair{Symbol,Float64}} = Pair{Symbol,Float64}[],\n  measurements::Vector{AbstractMeasurementPoint} = AbstractMeasurementPoint[],\n  events_active::Union{Nothing, Vector{Pair{Symbol,Bool}}} = Pair{Symbol,Bool}[],\n  events_save::Union{Tuple,Vector{Pair{Symbol, Tuple{Bool, Bool}}}}=(true,true), \n  observables::Union{Nothing,Vector{Symbol}} = nothing,\n  saveat::Union{Nothing,AbstractVector} = nothing,\n  tspan::Union{Nothing,Tuple} = nothing,\n  save_scope::Bool=true,\n  time_type::DataType=Float64,\n  kwargs...)\n\nSimulate model of type Model. Returns SimResults type. Example: sim(model; tspan = (0., 200.), parameters_upd = [:k1=>0.01])\n\nArguments:\n\nmodel : model of type Model\nparameters : Vector of Pairs containing constants' names and values. Overwrites default model's values. Default is empty vector \nmeasurements : Vector of measurements. Default is empty vector \nevents_active : Vector of Pairs containing events' names and true/false values. Overwrites default model's values. Default is empty vector \nevents_save : Tuple or Vector{Tuple} marking whether to save solution before and after event. Default is (true,true) for all events\nobservables : names of output observables. Overwrites default model's values. Default is empty vector\nsaveat : time points, where solution should be saved. Default nothing values stands for saving solution at timepoints reached by the solver \ntspan : time span for the ODE problem\nsave_scope : should scope be saved together with solution. Default is true\nkwargs : other solver related arguments supported by sim(cond::Cond)\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.sim-Tuple{Platform}","page":"API","title":"HetaSimulator.sim","text":"sim(platform::Platform; \n  conditions::Union{AbstractVector{Symbol}, Nothing} = nothing,\n  kwargs...) where {C<:AbstractCond}\n\nSimulate conditions included in platform. Returns Vector{Pair}. Example: sim(platform)\n\nArguments:\n\nplatform : platform of Platform type\nconditions : Vector containing names of conditions included in platform. Default value nothing stands for all conditions in the platform \nkwargs : other kwargs supported by sim(condition_pairs::Vector{Pair})\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.sim-Union{Tuple{AbstractVector{C}}, Tuple{C}} where C<:HetaSimulator.AbstractCond","page":"API","title":"HetaSimulator.sim","text":"sim(conditions::AbstractVector{C}; kwargs...) where {C<:AbstractCond}\n\nSimulate multiple conditions. Returns Vector{Pair}. Example: sim([cond1, cond2, cond3])\n\nArguments:\n\nconditions : Vector containing names and conditions of type Cond\nkwargs : other kwargs supported by sim(condition_pairs::Vector{Pair})\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.sim-Union{Tuple{Cond}, Tuple{P}} where P<:Pair","page":"API","title":"HetaSimulator.sim","text":"sim(cond::Cond; \n  parameters_upd::Union{Nothing, Vector{P}}=nothing,\n  alg=DEFAULT_ALG, \n  reltol=DEFAULT_SIMULATION_RELTOL, \n  abstol=DEFAULT_SIMULATION_ABSTOL,\n  kwargs...)\n\nSimulate single condition cond. Returns SimResults type. Example: Cond(model; tspan = (0., 200.), saveat = [0.0, 150., 250.]) |> sim\n\nArguments:\n\ncond : simulation condition of type Cond\nparameters_upd : constants, which overwrite both Model and Cond constants. Default is nothing\nalg : ODE solver. See SciML docs for details. Default is AutoTsit5(Rosenbrock23())\nreltol : relative tolerance. Default is 1e-3\nabstol : relative tolerance. Default is 1e-6\nkwargs : other solver related arguments supported by DiffEqBase.solve. See SciML docs for details\n\n\n\n\n\n","category":"method"},{"location":"api/#HetaSimulator.sim-Union{Tuple{Vector{P}}, Tuple{P}} where P<:Pair","page":"API","title":"HetaSimulator.sim","text":"sim(condition_pairs::Vector{P}; \n  parameters_upd::Union{Nothing, Vector}=nothing,\n  alg=DEFAULT_ALG, \n  reltol=DEFAULT_SIMULATION_RELTOL, \n  abstol=DEFAULT_SIMULATION_ABSTOL,\n  parallel_type=EnsembleSerial(),\n  kwargs...) where P<:Pair\n\nSimulate multiple conditions. Returns Vector{Pair}. Example: sim([:x => cond1, :y=>cond2, :z=>cond3])\n\nArguments:\n\ncondition_pairs : vector of pairs containing names and conditions of type Cond\nparameters_upd : constants, which overwrite both Model and Cond constants. Default is nothing\nalg : ODE solver. See SciML docs for details. Default is AutoTsit5(Rosenbrock23())\nreltol : relative tolerance. Default is 1e-3\nabstol : relative tolerance. Default is 1e-6\nparallel_type : type of multiple simulations parallelism. Default is no parallelism. See SciML docs for details\nkwargs : other solver related arguments supported by DiffEqBase.solve. See SciML docs for details\n\n\n\n\n\n","category":"method"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"HetaSimulator is a simulation and parameters estimation (fitting) platform for Heta modeling language. The main purpose of the platform is to establish the linkage between emerging QSP frameworks and fast computational methods (parallel simulations, automatic differentiation, etc.). HetaSimulator is inspired by the user experience of the software packages like SBMLToolbox, mrgsolve, DBSolve, dMod. From the computational point of view, it utilizes the unique features of Julia and SciML ecosystem.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"It is assumed that you have Julia v1.6 installed. Latest Julia release can be downloaded from julialang.org","category":"page"},{"location":"","page":"Home","title":"Home","text":"To install or update HetaSimulator and Heta compiler run:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ]\n(@v1.6) pkg> add https://github.com/hetalang/HetaSimulator.jl.git\njulia> using HetaSimulator\njulia> heta_update() # installs \"Heta compiler\" in NodeJS","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The user of HetaSimulator typically deals with the following three types:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Model - an ODE model, containing rhs, rules, initial parameters and vector of events.\nCond - condition representing a special model's setup for simulations or fitting. This setup can include initial parameters and events settings, output variables etc. In case of fitting Cond should also include experimental data. A common usage of Cond can be model's simulation with different drugs (parameters and events setup). Different Cond's can be united to run multi-conditional simulations and fitting.\nPlatform - container for different Models and Conds.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The user can perform the following three operations with both Model, Cond and Platform","category":"page"},{"location":"","page":"Home","title":"Home","text":"sim - run a single simulation or multi-conditional simulations. \nfit - fit a model to experimental data. \nmc - run Monte-Carlo or virtual patients simulations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"See documentation for detailed overview of HetaSimulator types and functions' arguments.","category":"page"},{"location":"#Basic-usage","page":"Home","title":"Basic usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A basic use-case example is provided in /cases/story_1 folder","category":"page"},{"location":"","page":"Home","title":"Home","text":"using HetaSimulator, Plots\n\nplatform = load_platform(\"$HetaSimulatorDir/cases/story_1\", rm_out=false);\nmodel = platform.models[:nameless]\n\n## single simulation\n\nsim(model; tspan = (0., 200.)) |> plot #1\n\n## condition simulation\n\ncond1 = Cond(model; tspan = (0., 200.), events_active=[:ss1 => false], saveat = [0.0, 150., 250.]);\nsim(cond1) |> plot\ncond2 = Cond(model; tspan = (0., 200.), events_active=[:sw1=>false, :ss1 => false], parameters = [:k2 => 0.001, :k3 => 0.02]);\nsim(cond2) |> plot\ncond3 = Cond(model; tspan = (0., 200.), events_active=[:ss1 => false], parameters = [:k1=>0.01]);\nsim(cond3) |> plot \n\nsim([:x => cond1, :y=>cond2, :z=>cond3]) |> plot\n\n## fitting\n\nmeasurements_csv = read_measurements(\"$HetaSimulatorDir/cases/story_1/measurements.csv\")\ncond4 = Cond(model; parameters = [:k2=>0.001, :k3=>0.04], events_active=[:ss1 => false], saveat = [0.0, 50., 150., 250.]);\nadd_measurements!(cond4, measurements_csv; subset = Dict(:condition => :dataone))\nres2 = fit([cond2, cond3, cond4], [:k1=>0.1,:k2=>0.2,:k3=>0.3])\n\n## Monte-Carlo simulations\n\nmccond1 = Cond(model; tspan = (0., 200.), parameters = [:k1=>0.01], saveat = [50., 80., 150.], events_active=[:ss1 => false]);\nmccond2 = Cond(model; tspan = (0., 200.), parameters = [:k1=>0.02], saveat = [50., 100., 200.], events_active=[:ss1 => false]);\nmccond3 = Cond(model; tspan = (0., 200.), parameters = [:k1=>0.03], saveat = [50., 100., 180.], events_active=[:ss1 => false]);\n\nmc(mccond1, [:k2=>Normal(1e-3,1e-4), :k3=>Normal(1e-4,1e-5)], 1000) |> plot\nmc([:mc1=>mccond1,:mc2=>mccond2,:mc3=>mccond3], [:k1=>0.01, :k2=>Normal(1e-3,1e-4), :k3=>Uniform(1e-4,1e-2)], 1000) |> plot\n\n## Simulations and fitting with Platform interface\n\n# load conditions\nconditions_csv = read_conditions(\"$HetaSimulatorDir/cases/story_1/conditions.csv\")\nadd_conditions!(platform, conditions_csv)\n\n# load measurements\nmeasurements = read_measurements(\"$HetaSimulatorDir/cases/story_1/measurements.csv\");\nadd_measurements!(platform, measurements)\n\nsim(platform, conditions = [:three]) |> plot\nfit1 = fit(platform, [:k1=>0.1,:k2=>0.2,:k3=>0.3], conditions = [:dataone, :withdata2])","category":"page"}]
}
