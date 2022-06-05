using OrdinaryDiffEq, LabelledArrays

function lorenz_f(du,u,p,t)
  du.x = p.σ*(u.y-u.x)
  du.y = u.x*(p.ρ-u.z) - u.y
  du.z = u.x*u.y - p.β*u.z
end

u0 = @LArray [1.0,0.0,0.0] (:x,:y,:z)
p = @LArray [10.0, 28.0, 8/3]  (:σ,:ρ,:β)
tspan = (0.0,10.0)

evt1 = DiscreteCallback((u,t,integrator)->t==5.0, (integrator)->integrator.u.x = 1.0)
evt2 = DiscreteCallback((u,t,integrator)->t==5.0, (integrator)->integrator.u.x = 1.2)

function myobserved(sym, u, p, t)
  if sym == :v
    return u.x * p.σ
  else
    throw("key $sym not found")
  end
end

f = ODEFunction(lorenz_f; observed=myobserved, syms=collect(keys(u0)))

prob = ODEProblem(f,u0,tspan,p,callback=CallbackSet(evt1,evt2), tstops=[5.0])
sol = solve(prob, Tsit5(), saveat=[1.0])
