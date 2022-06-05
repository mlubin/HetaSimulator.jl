using HetaSimulator

m = load_jlmodel("./_julia/model.jl")

s = Scenario(m, (0,100))

sol = sim(s)
