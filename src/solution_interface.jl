############################ Accessing Solution Values #####################

observables(sr::SimResults) = observables(sr.sim)
observables(sim::Simulation) = collect(keys(sim.vals.u[1]))
observables(sim::MCResults) = collect(keys(sim.u[1][1]))
constants(sim::SimResults) = sim.scenario.prob.p.constants

@inline Base.getindex(sim::Simulation, I...) = sim.vals[I...]
@inline Base.getindex(sim::Simulation, i::Symbol,::Colon) = [sim.vals[j][i] for j in 1:length(sim.vals)]
@inline Base.getindex(sr::SimResults, I...) = sr.sim[I...]
@inline Base.getindex(mc::MCResults, I...) = mc.sim[I...]

solat(sr::SimResults, args...) = solat(sr.sim, args...)

solat(sim::Simulation, t, idx, scope) = solat(sim, t, scope)[idx]
solat(sim::Simulation, t, idx::Colon, scope) = solat(sim, t, scope)

function solat(sim::Simulation, t, scope)
  if scope == :ode_ || scope == :start_
      _id = findfirst(x->x==t, sim.vals.t) # change to searchsortedfirst ?
  else
      t_id = findall(x->x==t, sim.vals.t)
      _id = t_id[findfirst(x->x==scope, @view(sim.scope[t_id]))]
  end
  return sim[_id]
end

(s::Simulation)(t, idx=:, scope=:ode_) = solat(s, t, idx, scope)
(sr::SimResults)(t, idx=:, scope=:ode_) = solat(sr.sim, t, idx, scope)

############################ Plots ########################################

function layout_choice(n)
  #=
  n == 1 && return (1,1)
  n == 2 && return (2,1)
  n == 3 && return (3,1)
  n == 4 && return (2,2)
  n > 4  && return n
  =#
  return (n,1)
end

@recipe function plot(sim::Simulation; vars = observables(sim))
  @assert !isempty(sim.vals) "Results don't contain output. You should probably add output observables to your model"

  time = sim.vals.t
  vals = [sim[id,:] for id in vars]
 
  xguide --> "time"
  yguide --> "output"
  label --> permutedims(string.(vars))
  xlims --> (time[1],time[end])
  linewidth --> 3
  (time, vals)
end


@recipe function plot(sr::SimResults; vars=observables(sr), measurements=true)

  @assert !isempty(sr.sim.vals) "Results don't contain output. You should probably add output observables to your model"

  time = sr.sim.vals.t
  vals = [sr[id,:] for id in vars]
 
  @series begin
    xguide --> "time"
    yguide --> "output"
    label --> permutedims(string.(vars))
    xlims --> (time[1],time[end])
    linewidth --> 3
    (time, vals)
  end

  if measurements == true && !isempty(sr.scenario.measurements)
    t_meas = NamedTuple{Tuple(vars)}([Float64[] for i in eachindex(vars)])
    vals_meas = NamedTuple{Tuple(vars)}([Float64[] for i in eachindex(vars)])
    for meas in sr.scenario.measurements
      μ = meas.μ 
      if isa(μ,Symbol) && μ ∈ vars 
        push!(t_meas[μ], meas.t)
        push!(vals_meas[μ], meas.val)
      end
    end
    for v in vars
      if !isempty(t_meas[v])
        @series begin
          seriestype --> :scatter
          xguide --> "time"
          yguide --> "output"
          label --> "$(v)"
          (t_meas[v], vals_meas[v])
        end
      end
    end
  end
  nothing
end
#= XXX: do we need it?
@recipe function plot(sim::Vector{S}) where S <:AbstractResults
  [Symbol("_$i")=>s for (i,s) in enumerate(sim)]
end
=#
@recipe function plot(s::Pair{Symbol,S}) where S<:AbstractResults
  @series begin
    title := "$(first(s))"
    last(s)
  end
end
@recipe function plot(sim::Vector{Pair{Symbol,S}}) where S<:AbstractResults
  (m,n) = layout_choice(length(sim))
  layout := (m,n)
  size := (400*n,300*m)
  for (i, s) in enumerate(sim)
    @series begin
      subplot := i
      s
    end
  end
end

#https://github.com/SciML/SciMLBase.jl/blob/7151bbe784df70cc572073d76d3a818aa8d1f4d0/src/ensemble/ensemble_solutions.jl#L102
@recipe function plot(sol::MCResults)
  for i in 1:length(sol)
    @series begin
      legend := false
      sol[i]
    end
  end
end

############################ DataFrames ########################################

function DataFrame(s::Simulation; vars=observables(s))
  df = DataFrame(t=s.vals.t)

  [df[!, v] = s[v,:] for v in vars[in.(vars, Ref(observables(s)))]]
  !isnothing(s.scope) && (df[!,:scope]=s.scope)
  
  return df
end

DataFrame(sr::SimResults; kwargs...) = DataFrame(sr.sim; kwargs...)

function DataFrame(sr::Pair{Symbol,S}; kwargs...) where S<:SimResults
  df = DataFrame(last(sr); kwargs...)
  df[!, :scenario] .= first(sr) # add new column

  return df
end

function DataFrame(res::Vector{Pair{Symbol,S}}; kwargs...) where S<:SimResults
  df_vectors = DataFrame.(res; kwargs...)

  return vcat(df_vectors...; cols=:union)
end

function DataFrame(mcr::MCResults; kwargs...)
  # df performance
  df = DataFrame()

  for (i,s) in enumerate(mcr.sim)
      dfs = DataFrame(s; kwargs...)
      insertcols!(dfs, 1, :iter => fill(i, length(s)))
      df = vcat(df,dfs)
  end
  
  return df
end

function DataFrame(mcr::Pair{Symbol,S}; kwargs...) where S<:MCResults
  df = DataFrame(last(mcr); kwargs...)
  df[!, :scenario] .= first(mcr) # add new column

  return df
end

function DataFrame(res::Vector{Pair{Symbol,S}}; kwargs...) where S<:MCResults
  df_vectors = DataFrame.(res; kwargs...)

  return vcat(df_vectors...; cols=:union)
end

############################ Save Results ########################################
"""
    save_results(filepath::String, sim::AbstractResults) 

Save results as csv file

Arguments:

- `filepath`: path and name of the file to write to
- `sim`: simulation results of `AbstractResults` type
"""
save_results(filepath::String, sim::AbstractResults) = save_results(filepath, DataFrame(sim))

save_results(filepath::String, df::DataFrame) = CSV.write(filepath, df, delim=";")

function save_optim(filepath::String, fr::FitResults)
  optim_params = optim(fr)
  open(filepath,"a") do io
    for op in optim_params
      println(io,"$(first(op)) = $(last(op));")
    end
 end
end

#=FIXME
function save_results(path::String, mcsim::MCResults; groupby::Symbol=:observables) 
  if groupby == :simulations
    for i in 1:length(mcsim)
      save_results("$path/$i.csv", mcsim[i])
    end
  elseif groupby == :observables
    obs = observables(mcsim[1])
    lobs = length(obs)
    for ob in obs
      df = DataFrame(t = mcsim[1].t)
      [df[!,string(i)] = mcsim[i][ob,:] for i in 1:length(mcsim)]
      CSV.write("$path/$ob.csv", df, delim=";")
    end
  end
end
=#