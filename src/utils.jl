
collect_saveat(saveat::Tuple{S1,S2}) where {S1<:Real,S2<:Real} = Float64[]
collect_saveat(saveat::Vector{S}) where S<:Real = Float64.(saveat)
collect_saveat(saveat::AbstractRange{S}) where S<:Real = Float64.(saveat)

dictkeys(d::Dict) = (collect(keys(d))...,)
dictvalues(d::Dict) = (collect(values(d))...,)

#=
to_nt(d::Dict{Symbol,T}) where {T} =
    NamedTuple{dictkeys(d)}(dictvalues(d))
=#
typed(v::Number) = Float64(v)
function typed(v::String) 
  v_float = tryparse(Float64,v)
  ret = isnothing(v_float) ? Symbol(v) : v_float
  ret
end

typed(v::Symbol) = v

function _subset(vector::AbstractVector, subset::Union{Dict{Symbol}, Nothing} = nothing)
  if subset === nothing
    # use the whole vector
    return vector
  else
    # this function performs matching
    # returs true if all values in subset equal to values in row
    subset_fun = (row) -> begin
      res = [row[key] == value for (key, value) in subset]
      return all(res)
    end

    return filter(subset_fun, vector)
  end
end

has_saveat(c::Cond) = !isnothing(c.saveat) && length(c.saveat) != 0

# auxilary method to update values in LVector
# TODO: not optimal method
# TODO: write for any number of pairs update(la::LArray, argv...)
function update(la::LArray, pairs::AbstractVector{Pair{Symbol, C}}) where C<:Real
  clone = copy(la)
  keys_la = keys(la)

  for (key, value) in pairs
    key ∉ keys_la && throw("key :$key is not found in LVector.")
    clone[key] = value
  end

  return clone
end

update(p1::AbstractVector{P}, p2::AbstractVector{P}) where P<:Pair = update(Dict(p1),Dict(p2))

function update(d1::Dict, d2::Dict)
  d = copy(d1)
  keys_d = keys(d)

  for (k,v) in d2
    k ∉ keys_d && throw("key :$k is not found in Dict.")
    d[k] = v
  end

  return d
end