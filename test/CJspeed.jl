import SDtoolbox: CJspeed
using Plots

P₁ = 1e5
T₁ = 300
ϕ = 1

X₁ = "H2:$(42*ϕ), O2:21,N2:79"
try
    global mech = "gri30.yaml"
catch LoadError
    global mech = "gri30.xml"
end

# reference value from online NASA CEA
V_CJ_ref = 1964.6

# compare results for all CJ algorithms
algorithms = "CEA", "umin", "aeq"
V_CJ = Float64[]

for method in algorithms
    push!(V_CJ,CJspeed(P₁,T₁,X₁,mech,method=method))
end

rel_err = maximum(abs.(1 .- V_CJ / (sum(V_CJ)/3) ))
ref_err = maximum(abs.(V_CJ_ref .- V_CJ))/V_CJ_ref
