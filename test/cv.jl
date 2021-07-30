import SDtoolbox: ct, cvsolve
using Plots
using PyCall
using StatsBase
using Interpolations

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

# compare with original sdtoolbox
pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__)
cvsolve_py = pyimport("sdtoolbox.cv").cvsolve
gas.TPX = T₁,P₁,X₁
out_py = cvsolve_py(gas,t_end=out["time"][end])
plot!(out_py["time"],out_py["T"])

# calculate cross correlation to quantify error
itp = LinearInterpolation(out["time"][1:(end-1)], out["T"][1:(end-1)])
sol_itp = itp(out_py["time"])

err = maximum(crosscor(sol_itp,out_py["T"]))
