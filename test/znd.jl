import SDtoolbox: ct, PostShock_fr, CJspeed, zndsolve
using Plots

P₁ = 1e5
T₁ = 300

# for ϕ = 0.5:0.05:1.1
    ϕ = 1

    X₁ = "H2:$(42*ϕ), O2:21,N2:79"
    mech = "gri30.xml"

    gas₁ = ct.Solution(mech)::PyObject
    gas₁.TPX = T₁,P₁,X₁

    U₁ = CJspeed(P₁,T₁,X₁,mech)::Float64

    gas = PostShock_fr(U₁, P₁, T₁, X₁, mech)::PyObject

    out = zndsolve(gas,gas₁,U₁,advanced_output=true)

    plot!(out["time"],out["thermicity"])
# end
