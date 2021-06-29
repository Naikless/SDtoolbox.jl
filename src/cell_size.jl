# to be able to run this standalone via include()
using PyCall
import SDtoolbox: zndsolve, cvsolve, CJspeed, PostShock_fr, PostShock_eq, ct

include("correlations.jl")

function cell_size(T₁,P₁,X₁::Union{AbstractString,AbstractDict},mech;kwargs...)

    # Find CJ speed
    cj_speed = CJspeed(P₁,T₁,X₁,mech)::Float64
    println("CJ velocity: $(round(cj_speed,digits=3)) m/s")

    # Set up gas object
    gas₁ = ct.Solution(mech)::PyObject
    gas₁.TPX = T₁,P₁,X₁

    # Find equilibrium post shock state for given speed
    gas = PostShock_eq(cj_speed, P₁, T₁, X₁, mech)::PyObject
    u_cj = cj_speed*gas₁.density/gas.density
    println("u_CJ: $(round(u_cj,digits=3)) m/s")

    # Find frozen post shock state for given speed
    gas = PostShock_fr(cj_speed, P₁, T₁, X₁, mech)::PyObject
    Tₛ,Pₛ = gas.TP::Tuple{Float64, Float64}
    println("Tₛ: $(round(Tₛ,digits=3)), Pₛ: $(round(Pₛ,digits=3))")

    # Solve ZND ODEs
    out = zndsolve(gas,gas₁,cj_speed,advanced_output=true;kwargs...)::AbstractDict

    # Find CV parameters including effective activation energy
    T₊ = Tₛ*1.02
    gas.TPX = T₊,Pₛ,X₁
    CV_out₊ = cvsolve(gas;kwargs...)::AbstractDict
    T₋ = Tₛ*0.98
    gas.TPX = T₋,Pₛ,X₁
    CV_out₋ = cvsolve(gas;kwargs...)::AbstractDict
    # Approximate effective activation energy for CV explosion
    τ₊ = CV_out₊["ind_time"]
    τ₋ = CV_out₋["ind_time"]
    if τ₊ == zero(τ₊) && τ₋ == zero(τ₋)
        θₑ_CV = 0.0
    else
        θₑ_CV = 1.0/Tₛ*((log(τ₊)-log(τ₋))/((1.0/T₊)-(1.0/T₋)))
    end

    Δᵣ = u_cj / out["max_thermicity_ZND"]
    Δᵢ = out["ind_len_ZND"]
    χ = θₑ_CV * Δᵢ / Δᵣ

    # obtain cell size through correlation
    λ = ng(Δᵢ,χ)::Float64
end
