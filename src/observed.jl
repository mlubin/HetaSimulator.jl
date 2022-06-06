#=
struct StaticWithCache{T,C} <: AbstractVector{T}
  cur::AbstractVector{T}
  cache::C
end


Base.size(s::StaticWithCache) = size(s.cur)
Base.getindex(s::StaticWithCache, i) = s.cur[i]
Base.setindex!(s::StaticWithCache, val, i) = (s.cur[i] = val)
=#

struct StaticCache{I,T} <: AbstractVector{T}
  idxs::I
  cache::AbstractVector{T}
end


Base.size(s::StaticCache) = !isempty(s.idxs) ? (s.idxs[end],) : (0,)
function Base.getindex(s::StaticCache, j::Int)
  (j > length(s)) && throw(BoundsError())
  
  i = searchsortedlast(s.idxs, j)
  return s.cache[i]
end

const ParamsWithCache{C,S,CA} = Params{C,S,CA} where CA<:StaticCache
Base.getindex(p::ParamsWithCache, i::Int) = Params(p.constants, p.static_cache[i], nothing)
Base.getindex(p::ParamsWithCache, i::AbstractVector{<:Integer}) = Params.((p.constants,), p.static_cache[i], (nothing,))
Base.getindex(p::ParamsWithCache, ::Colon) = Params.((p.constants,), p.static_cache, (nothing,))

const ODEProblemWithStaticCache{uType,tType,isinplace,P,F,K,PT} = ODEProblem{uType,tType,isinplace,P,F,K,PT} where P<:ParamsWithCache
const ODESolutionWithStaticCache{T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE} = Union{
  OrdinaryDiffEq.ODECompositeSolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE},
  SciMLBase.ODESolution{T,N,uType,uType2,EType,tType,rateType,P,A,IType,DE}} where P <: ODEProblemWithStaticCache
const DiffEqArrayWithStaticCache{T, N, A, B, C, D, E, F} = RecursiveArrayTools.DiffEqArray{T, N, A, B, C, D, E, F} where F<: Union{ParamsWithCache, Vector{P}} where P <: Params

function RecursiveArrayTools.observed(A::DiffEqArrayWithStaticCache,sym,i::Int)
  A.observed(sym,A.u[i],A.p[i],A.t[i])
end

function RecursiveArrayTools.observed(A::DiffEqArrayWithStaticCache,sym,i::AbstractArray{Int})
  A.observed.((sym,),A.u[i],A.p[i],A.t[i])
end

function RecursiveArrayTools.observed(A::DiffEqArrayWithStaticCache,sym,::Colon)
  A.observed.((sym,),A.u,A.p[:],A.t)
end

function SciMLBase.observed(A::ODESolutionWithStaticCache,sym,i::Int)
  SciMLBase.getobserved(A)(sym,A[i],A.prob.p[i],A.t[i])
end

function SciMLBase.observed(A::ODESolutionWithStaticCache,sym,i::AbstractArray{Int})
  SciMLBase.getobserved(A).((sym,),A.u[i],A.prob.p[i],A.t[i])
end

function SciMLBase.observed(A::ODESolutionWithStaticCache,sym,i::Colon)
  SciMLBase.getobserved(A).((sym,),A.u,A.prob.p[:],A.t)
end

#static(p::Params) = p.static
#cache(p::Params) = p.cache

#=
(p::Params)(t) = _interpolate(p.static_cache,t)

# https://github.com/PumasAI/DataInterpolations.jl/blob/447bab144d3a27b90edd8de3b5df18ab5215b309/src/interpolation_methods.jl#L119
function _interpolate(static_cache,t)
  i = searchsortedlast(static_cache.t, t)
  return static_cache.p[max(1, i)]
end
=#

# https://github.com/SciML/SciMLBase.jl/blob/0274b5f9f42acbc98cb4fc57bab997bddcb3050c/src/solutions/ode_solutions.jl#L43
(sol::ODESolutionWithStaticCache)(t, ::Type{deriv}=Val{0}; idxs=nothing, continuity=:left) where {deriv} = sol(t, deriv, idxs, continuity)

function (sol::ODESolutionWithStaticCache)(t::Number, ::Type{deriv}, idxs, continuity) where {deriv}
  SciMLBase.issymbollike(idxs) || error("Incorrect specification of `idxs`")
  interp_p = static_interpolation(sol,t,continuity)
  _augment(sol.interp([t], nothing, deriv, interp_p, continuity), sol, interp_p)[idxs][1]
end

function (sol::ODESolutionWithStaticCache)(t::Number, ::Type{deriv}, idxs::AbstractVector, continuity) where {deriv}
  all(SciMLBase.issymbollike.(idxs)) || error("Incorrect specification of `idxs`")
  interp_p = static_interpolation(sol,t,continuity)
  interp_sol = _augment(sol.interp([t], nothing, deriv, interp_p, continuity), sol, interp_p)
  [first(interp_sol[idx]) for idx in idxs]
end

function (sol::ODESolutionWithStaticCache)(t::AbstractVector{<:Number}, ::Type{deriv}, idxs, continuity) where {deriv}
  SciMLBase.issymbollike(idxs) || error("Incorrect specification of `idxs`")
  interp_p = static_interpolation(sol,t,continuity)
  interp_sol = _augment(sol.interp(t, nothing, deriv, interp_p, continuity), sol, interp_p)
  observed = SciMLBase.has_observed(sol.prob.f) ? sol.prob.f.observed : SciMLBase.DEFAULT_OBSERVED
  RecursiveArrayTools.DiffEqArray(interp_sol[idxs], t, [idxs], SciMLBase.getindepsym(sol), observed, interp_p)
end

function (sol::ODESolutionWithStaticCache)(t::AbstractVector{<:Number}, ::Type{deriv}, idxs::AbstractVector, continuity) where {deriv}
  all(SciMLBase.issymbollike.(idxs)) || error("Incorrect specification of `idxs`")
  interp_p = static_interpolation(sol,t,continuity)
  interp_sol = _augment(sol.interp(t, nothing, deriv, interp_p, continuity), sol, interp_p)
  observed = SciMLBase.has_observed(sol.prob.f) ? sol.prob.f.observed : SciMLBase.DEFAULT_OBSERVED
  RecursiveArrayTools.DiffEqArray([[interp_sol[idx][i] for idx in idxs] for i in 1:length(t)], t, idxs, SciMLBase.getindepsym(sol), observed, interp_p)
end

function _augment(A::RecursiveArrayTools.DiffEqArray, sol::ODESolutionWithStaticCache, p)
  syms = hasproperty(sol.prob.f, :syms) ? sol.prob.f.syms : nothing
  observed = SciMLBase.has_observed(sol.prob.f) ? sol.prob.f.observed : SciMLBase.DEFAULT_OBSERVED
  RecursiveArrayTools.DiffEqArray(A.u, A.t, syms, SciMLBase.getindepsym(sol),observed,p)
end

function static_interpolation(sol::ODESolutionWithStaticCache, t::Number, continuity)
  if continuity == :left
    i = searchsortedfirst(sol.t, t)
  else
    i = searchsortedlast(sol.t, t)
  end
  return sol.prob.p[i]
end

function static_interpolation(sol::ODESolutionWithStaticCache, t::AbstractVector{<:Number}, continuity)
  if continuity == :left
    i = searchsortedfirst.((sol.t,), t)
  else
    i = searchsortedlast.((sol.t,), t)
  end
  return sol.prob.p[i]
end