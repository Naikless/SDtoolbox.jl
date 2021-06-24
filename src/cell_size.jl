using SDtoolbox
using PyCall
ct = pyimport_conda("cantera","cantera","cantera")
pushfirst!(PyVector(pyimport("sys")."path"), "") # adds local dir to python path
postshock = pyimport("sdtoolbox.postshock")
CJspeed = postshock.CJspeed
PostShock_fr = postshock.PostShock_fr
PostShock_eq = postshock.PostShock_eq

include("correlations.jl")

function cell_size(T₁::Real,P₁::Real,X₁::Union{AbstractString,AbstractDict},mech::AbstractString)

    # Find CJ speed
    cj_speed = CJspeed(P₁,T₁,X₁,mech)::Float64

    # Set up gas object
    gas₁ = ct.Solution(mech)::PyObject
    gas₁.TPX = T₁,P₁,X₁

    # Find equilibrium post shock state for given speed
    gas = PostShock_eq(cj_speed, P₁, T₁, X₁, mech)::PyObject
    u_cj = cj_speed*gas₁.density/gas.density

    # Find frozen post shock state for given speed
    gas = PostShock_fr(cj_speed, P₁, T₁, X₁, mech)::PyObject
    Tₛ,Pₛ = gas.TP::Tuple{Float64, Float64}

    # Solve ZND ODEs
    out = zndsolve(gas,gas₁,cj_speed,advanced_output=true)::AbstractDict

    # Find CV parameters including effective activation energy
    T₊ = Tₛ*1.02
    gas.TPX = T₊,Pₛ,X₁
    CV_out₊ = cvsolve(gas)
    T₋ = Tₛ*0.98
    gas.TPX = T₋,Pₛ,X₁
    CV_out₋ = cvsolve(gas)
    # Approximate effective activation energy for CV explosion
    τ₊ = CV_out₊["ind_time"]
    τ₋ = CV_out₋["ind_time"]
    if τ₊ == 0 && τ₋ == 0
        θₑ_CV = 0
    else
        θₑ_CV = 1/Tₛ*((log(τ₊)-log(τ₋))/((1/T₊)-(1/T₋)))
    end

    Δᵣ = u_cj / out["max_thermicity_ZND"]
    Δᵢ = out["ind_len_ZND"]
    χ = θₑ_CV * Δᵢ / Δᵣ

    # obtain cell size through correlation
    λ = ng(Δᵢ,χ)
end
