const NORMAL = :normal
const LOGNORMAL = :lognormal
const SIGMA = :sigma
const MEAN = :mean

# CSV methods

function add_measurements!(
  scenario::Scenario,
  vector::AbstractVector;
  subset::AbstractVector{P} = Pair{Symbol, Symbol}[]
) where P <: Pair{Symbol, Symbol}
  selected_rows = _subset(vector, subset)

  for row in selected_rows
    _add_measurement!(scenario, row)
  end
end

function add_measurements!(
  platform::Platform,
  vector::AbstractVector;
  subset::AbstractVector{P} =  Pair{Symbol, Symbol}[]
) where P <: Pair{Symbol, Symbol}
  selected_rows = _subset(vector, subset)

  # we will store here lost scenario names
  lost_refs = Symbol[]

  for row in selected_rows
    scenario_ref = row[:scenario]

    if !haskey(platform.scenarios, scenario_ref)
      push!(lost_refs, scenario_ref)
    else
      scenario = platform.scenarios[scenario_ref]
      _add_measurement!(scenario, row)
    end
  end

  if length(lost_refs) > 0
    @warn "Lost scenario names: $(unique(lost_refs)). Some measurement points haven't been added."
  end

  return nothing
end

# DataFrame methods
"""
    add_measurements!(
      scenario::Scenario,
      df::DataFrame;
      kwargs...
    )

Adds measurements to `Scenario`

Arguments:

- `scenario` : simulation scenario of type [`Scenario`](@ref)
- `df` : `DataFrame` with measurements, typically obtained with [`read_measurements`](@ref) function
- `subset` : subset of measurements which will be added to the `Scenario`. Default `Pair{Symbol, Symbol}[]` adds all measurements from the `df`
"""
function add_measurements!(
  scenario::Scenario,
  df::DataFrame;
  kwargs...
)
  add_measurements!(scenario, eachrow(df); kwargs...)
end

"""
    add_measurements!(
      platform::Platform,
      df::DataFrame;
      kwargs...
    )

Adds measurements to `Scenario`

Arguments:

- `platform` : platform of [`Platform`](@ref) type
- `df` : `DataFrame` with measurements, typically obtained with [`read_measurements`](@ref) function
- `subset` : subset of measurements which will be added to the `Scenario`. Default `Pair{Symbol, Symbol}[]` adds all measurements from the `df`
"""
function add_measurements!(
  platform::Platform,
  df::DataFrame;
  kwargs...
)
  add_measurements!(platform, eachrow(df); kwargs...)
end

# private function to add one measurement row into scenario 

function _add_measurement!(scenario::Scenario, row::Any) # maybe not any
  _t = row[:t]
  _val = row[:measurement]
  _scope = haskey(row, :scope) ? row[:scope] : missing

  _type = haskey(row, Symbol("prob.type")) ? row[Symbol("prob.type")] : missing
  type = ismissing(_type) ? NORMAL : _type

  if type in [NORMAL, LOGNORMAL]
    _mean = typed(row[Symbol("prob.$MEAN")])
    _sigma = typed(row[Symbol("prob.$SIGMA")])

    point = (type == LOGNORMAL) ? LogNormalMeasurementPoint(_t, _val, _scope, _mean, _sigma) : NormalMeasurementPoint(_t, _val, _scope, _mean, _sigma)
  else 
    error("Distribution value $type is wrong or not supported. Supported distributions are: $Normal, $Lognormal")
  end

  push!(scenario.measurements, point)
end

# helper to read from csv and xlsx

function read_measurements_csv(filepath::String; kwargs...)
  df = DataFrame(CSV.File(
    filepath;
    typemap = Dict(Int64=>Float64, Int32=>Float64),
    types = Dict(:t=>Float64, :measurement=>Float64, :scope=>Symbol, :scenario=>Symbol, Symbol("prob.type")=>Symbol),
    kwargs...)
  )
  assert_measurements(df)
  
  return df
end

function read_measurements_xlsx(filepath::String, sheet=1; kwargs...)
  df = DataFrame(XLSX.readtable(filepath, sheet,infer_eltypes=true)...)
  assert_measurements(df)

  df[!,:t] .= Float64.(df[!,:t])
  df[!,:measurement] .= Float64.(df[!,:measurement])
  hasproperty(df, :scope) && (df[!,:scope] .= Symbol.(df[!,:scope]))
  df[!,:scenario] .= Symbol.(df[!,:scenario])
  hasproperty(df, Symbol("prob.type")) && (df[!,Symbol("prob.type")] .= Symbol.(df[!,Symbol("prob.type")]))

  return df
end

"""
    read_measurements(filepath::String, sheet=1; kwargs...)

Reads table file with measurements to `DataFrame`

Arguments:

- `filepath` : path to table file. Supports ".csv" and ".xlsx" files
- `sheet` : number of sheet in case of ".xlsx" file. Default is `1`
- kwargs : other arguments supported by `CSV.File`
"""
function read_measurements(filepath::String, sheet=1; kwargs...)
  ext = splitext(filepath)[end]

  if ext == ".csv"
    df = read_measurements_csv(filepath; kwargs...)
  elseif ext == ".xlsx"
    df = read_measurements_xlsx(filepath, sheet)
  else  
    error("Extension $ext is not supported.")
  end
  return df
end

function assert_measurements(df)
  names_df = names(df)
  for f in ["t", "measurement", "scenario"]
    @assert f ∈ names_df "Required column name $f is not found in measurements table."
  end
  return nothing
end

#################################### GET MEASUREMENTS ###########################

# ! Current assumption is that error dist is normal or lognormal

measurements_as_table(sr::SimResult) = measurements_as_table(sr.scenario)

function measurements_as_table(scn::Scenario)
  meas = measurements(scn)
  @assert !isempty(meas) "Scenario doesn't contain measurements."

  lm = length(meas)
  _t = Vector{Float64}(undef,lm)
  _mea = []
  _scope = Vector{Symbol}(undef,lm)
  _mu = []
  _sigma = []
  _type = Vector{Symbol}(undef,lm)
  
  for (i,m) in enumerate(meas)
    _t[i] = m.t
    push!(_mea, m.val)
    _scope[i] = m.scope
    push!(_mu, m.μ)
    push!(_sigma, m.σ)
    _type[i] = typeof(m) <: NormalMeasurementPoint ? NORMAL : LOGNORMAL
  end
  df = DataFrame(t =_t, measurement = _mea, scope = _scope)
  df[!,"prob.mean"] = _mu
  df[!,"prob.sigma"] = _sigma 
  df[!,"prob.type"] = _type
  return df
end