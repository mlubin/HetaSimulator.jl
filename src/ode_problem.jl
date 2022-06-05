function build_ode_problem(
  model::Model,
  tspan;
  parameters::Vector{Pair{Symbol,Float64}} = Pair{Symbol,Float64}[],
  events_active::Union{Nothing, Vector{Pair{Symbol,Bool}}} = Pair{Symbol,Bool}[],
  events_save::Union{Tuple,Vector{Pair{Symbol, Tuple{Bool, Bool}}}} = (true,true)
)

  # initial values and params
  init_func = model.init_func
  parameters_nt = NamedTuple(parameters)
  u0, params = init_values(init_func, merge(model.constants, parameters_nt))

  # events
  ev_on_nt = !isnothing(events_active) ? NamedTuple(events_active) : NamedTuple()
  events = active_events(model.events, merge(model.events_active, ev_on_nt), events_save)
  cbs = CallbackSet(events...)
  
  # ODEFunction 
  ode_func = ODEFunction(model.ode_func; observed=model.saving_generator, syms=collect(keys(u0)))
  
  # problem setup
  return ODEProblem(
    ode_func, # ODEFunction
    u0, # u0
    tspan, # tspan
    params; # constants and static
    callback = cbs # callback
  )
end

#=
function saveat_tspan(saveat, tspan, time_type)

  if !isnothing(saveat) && !isempty(saveat)
    _saveat = collect_saveat(saveat)
    _tspan = (zero(time_type), time_type(maximum(_saveat)))
  elseif !isnothing(tspan)
    _saveat = time_type[]
    _tspan = (zero(time_type), time_type(last(tspan))) # tspan should start from zero?
  else
    error("Please, provide either `saveat` or `tspan` value.")
  end
  return (_saveat, _tspan)
end


collect_saveat(saveat::Tuple) = Float64[]
collect_saveat(saveat::Vector{S}) where S<:Real = Float64.(saveat)
collect_saveat(saveat::AbstractRange{S}) where S<:Real = Float64.(saveat)
=#

function init_values(init_func, constants)
  u0, p0 = init_func(constants)
  lvcons = LVector(constants)
  lvstatics = LVector(p0)
  lvu0 = LVector(u0)
  
  return (lvu0, Params(lvcons, lvstatics, StaticCache([1], [copy(lvstatics)])))
end

# NamedTuple{}(fill(typeof(3)[], length(p0)))
#=
function  update_init_values(prob, init_func, x)
  constants = merge(NamedTuple(prob.p.constants),x)
  u0, p0 = init_values(init_func, constants)

  prob_upd = remake(prob; u0=u0, p=p0)
  
  return prob_upd
end
=#