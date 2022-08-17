import SDtoolbox: ct, PostShock_fr, CJspeed, zndsolve, CANTERA_MECH_FILETYPE
using PyCall
using StatsBase
using Interpolations

P₁ = 1e5
T₁ = 300


ϕ = 1

X₁ = "H2:$(42*ϕ), O2:21,N2:79"
mech = "gri30."*CANTERA_MECH_FILETYPE

gas₁ = ct.Solution(mech)
gas₁.TPX = T₁,P₁,X₁

U₁ = CJspeed(P₁,T₁,X₁,mech)

gas = PostShock_fr(U₁, P₁, T₁, X₁, mech)

out = zndsolve(gas,gas₁,U₁,advanced_output=true)

# Tplot = plot(out["time"],out["T"])
# therm_plot = plot(out["time"],out["thermicity"])


# compare with original sdtoolbox
pyimport_conda("scipy","scipy<=1.8")
pyimport_conda("matplotlib","matplotlib")
pushfirst!(PyVector(pyimport("sys")."path"), @__DIR__)
zndsolve_py = pyimport("sdtoolbox.znd").zndsolve
gas = pyimport("sdtoolbox.postshock").PostShock_fr(U₁,P₁,T₁,X₁,mech)
out_py = zndsolve_py(gas,gas₁,U₁,t_end=out["time"][end])
# plot!(Tplot, out_py["time"],out_py["T"])
# plot!(therm_plot, out_py["time"],out_py["thermicity"])

# calculate cross correlation to quantify error
itp_T = LinearInterpolation(out["time"][1:(end-1)], out["T"][1:(end-1)])
sol_itp_T = itp_T(out_py["time"])
itp_therm = LinearInterpolation(out["time"][1:(end-1)], out["thermicity"][1:(end-1)])
sol_itp_therm = itp_therm(out_py["time"])

err_T = maximum(crosscor(sol_itp_T,out_py["T"]))
err_therm = maximum(crosscor(sol_itp_therm,out_py["thermicity"]))
