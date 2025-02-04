function build_ode_problem(
  model::Model,
  tspan;
  parameters::Vector{Pair{Symbol,Float64}} = Pair{Symbol,Float64}[],
  events_active::Union{Nothing, Vector{Pair{Symbol,Bool}}} = Pair{Symbol,Bool}[],
  events_save::Union{Tuple,Vector{Pair{Symbol, Tuple{Bool, Bool}}}} = (true,true), 
  observables_::Union{Nothing,Vector{Symbol}} = nothing,
  saveat::Union{Nothing,AbstractVector} = nothing,
  save_scope::Bool = true,
  time_type::DataType = Float64
)
  
  _saveat = isnothing(saveat) ? time_type[] : saveat 
  # initial values and params
  init_func = model.init_func
  parameters_nt = NamedTuple(parameters)
  u0, params = init_values(init_func, merge(model.constants, parameters_nt))

  # check observables
  if !isnothing(observables_)
    records_ind = indexin(observables_, records(model))
    if any((x)-> x === nothing, records_ind)
      lost_observables_ind = findall((x)-> x===nothing, records_ind)
      lost_observables = observables_[lost_observables_ind]
      error("The following observables have not been found in the model: $lost_observables")
    end
  end

  # saving setup
  utype = eltype(u0)
  merged_observables = isnothing(observables_) ? observables(model) : observables_ # use default if not set
  saved_values = SavedValues(
    LArray{utype,1,Array{utype,1},Tuple(merged_observables)}[],
    time_type[],
    save_scope ? Symbol[] : nothing
    )
  saving_func = model.saving_generator(merged_observables)
  scb = saving_wrapper(saving_func, saved_values; saveat=_saveat, save_scope)

  # events
  ev_on_nt = !isnothing(events_active) ? NamedTuple(events_active) : NamedTuple()
  events = active_events(model.events, merge(model.events_active, ev_on_nt), events_save)
  cbs = CallbackSet(scb, events...)

  # problem setup
  return ODEProblem(
    model.ode_func, # ODE function
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
  return (u0,Params{typeof(lvcons),typeof(p0)}(lvcons, p0))
end

function  update_init_values(prob, init_func, x)
  constants = merge(NamedTuple(prob.p.constants),x)
  u0, p0 = init_values(init_func, constants)

  prob_upd = remake(prob; u0=u0, p=p0)
  
  return prob_upd
end

function remake_saveat(prob, saveat; tspan = prob.tspan)
  
  scb_orig = prob.kwargs[:callback].discrete_callbacks[1].affect!
  utype = eltype(prob.u0)
  save_scope = scb_orig.save_scope

  saved_values = SavedValues(
    LArray{utype,1,Array{utype,1},observables(prob)}[],
    eltype(tspan)[],
    save_scope ? Symbol[] : nothing
    )
  save_func = scb_orig.save_func
  scb_new = saving_wrapper(save_func, saved_values; saveat, save_scope=save_scope)

  cbs = list_callbacks(prob)
  cb_set = CallbackSet(scb_new,cbs[1]...,cbs[2]...)
  remake(prob; callback=cb_set, tspan=tspan)
end

function remake_saveat(prob, measurements::Vector{M}) where M<:AbstractMeasurementPoint
  saveat = unique([dp.t for dp in measurements])
  tspan = (prob.tspan[1], eltype(prob.tspan)(maximum(saveat)))
  remake_saveat(prob, saveat; tspan)
end

function list_callbacks(prob)
  discr_cbs = prob.kwargs[:callback].discrete_callbacks
  cont_cbs = prob.kwargs[:callback].continuous_callbacks
  (discr_cbs[1:end .!=1], cont_cbs)
end

observables(prob::SciMLBase.AbstractODEProblem) = observables(
  eltype(prob.kwargs[:callback].discrete_callbacks[1].affect!.saved_values.u)
)

observables(::Type{LArray{T,N,D,Syms}}) where {T,N,D,Syms} = Syms
