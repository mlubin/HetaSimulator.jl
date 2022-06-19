using HetaSimulator

m = load_jlmodel("./_julia/model.jl")

s = Scenario(m, (0,100))

sol = sim(s, alg=CVODE_BDF())

mcsol = mc(s, [:k2=>Normal(1e-3,1e-4), :k3=>Normal(1e-4,1e-5)], 100)