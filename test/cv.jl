import SDtoolbox: ct, cvsolve
using Plots
using PyCall
using StatsBase
using Interpolations

P₁ = 2e5
T₁ = 1500


ϕ = 1

X₁ = "H2:$(42*ϕ), O2:21,N2:79"
try
    global mech = "gri30.yaml"
catch LoadError
    global mech = "gri30.xml"
end

gas = ct.Solution(mech)
gas.TPX = T₁,P₁,X₁

out = cvsolve(gas,t_end=1e-3)

plot!(out["time"],out["T"])


# compare with original sdtoolbox
pyimport_conda("scipy","scipy")
pyimport_conda("matplotlib","matplotlib")
pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__)
cvsolve_py = pyimport("sdtoolbox.cv").cvsolve
gas.TPX = T₁,P₁,X₁
out_py = cvsolve_py(gas,t_end=out["time"][end])
plot!(out_py["time"],out_py["T"])

# calculate cross correlation to quantify error
itp = LinearInterpolation(out["time"][1:(end-1)], out["T"][1:(end-1)])
sol_itp = itp(out_py["time"])

err = maximum(crosscor(sol_itp,out_py["T"]))
