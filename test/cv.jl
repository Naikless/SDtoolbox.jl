import SDtoolbox: ct, cvsolve
using Plots

P₁ = 2e5
T₁ = 1500

# for ϕ = 0.5:0.05:1.1
    ϕ = 1

    X₁ = "H2:$(42*ϕ), O2:21,N2:79"
    mech = "gri30.xml"

    gas = ct.Solution(mech)
    gas.TPX = T₁,P₁,X₁

    out = cvsolve(gas,t_end=1e-3)

    plot!(out["time"],out["T"])
# end
