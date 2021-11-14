"""
    gavrikov(δ,θ,Tvn,T₀)

Correlation function for detonation cell width proposed by Gavrikov et al
COMBUSTION AND FLAME 120:19�33 (2000) based on using a reaction zone length based on time
to limiting reactant consumption in constant volume explosion approximation using vn
postshock velocity to convert time to distance.
Tested against a range of fuel-oxidizer diluent mixtures

INPUT:
    δ = reaction zone length based on time to 50consumption of limiting
    reactant from CV computation and delta = time * w_VN
    θ = Eₐ/RT_VN,  effective reduced activation energy based on CV
    computation
    Tvn = von Neumann (postshock temperature behind CJ shock wave)
    T₀ = initial temperature
"""
function gavrikov(δ,θ,Tvn,T₀)

    # Constants
    a = -0.007843787493
    b = 0.1777662961
    c = 0.02371845901
    d = 1.477047968
    e = 0.1545112957
    f = 0.01547021569
    g = -1.446582357
    h = 8.730494354
    i = 4.599907939
    j = 7.443410379
    k = 0.4058325462
    m = 1.453392165
    # Define nondimensional parameters
    X = θ
    Y = Tvn/T₀
    z = Y*(a*Y-b) + X*(c*X-d + (e-f*Y)*Y) + g*log(Y) + h*log(X) + Y*(i/X - k*Y/X^m) - j
    λ = δ*10^z
end

"""
     ng(Δᵢ,χ)

Correlation function for detonation cell size from
Ng, Hoi Dick, Yiguang Ju, and John H. S. Lee. 2007. Assessment of Detonation Hazards in
High-Pressure Hydrogen Storage from Chemical Sensitivity Analysis.
INTERNATIONAL JOURNAL OF HYDROGEN ENERGY 32 (1):93-99.
Tested only against low pressure H2-air data.

INPUT:
    Δᵢ = induction zone length based on peak thermicity in ZND simulation
    χ = 𝜺ᵢ*Δᵢ/Δᵣ where
          𝜺ᵢ = reduced effective activation energy from CV computation
          Δᵢ = distance to peak thermicity from ZND computation
          Δᵣ = w_vN/σ̇_max from ZND computation

See Ng et al. Combustion Theory and Modeling 2005 for a discussion of
the χ parameter.
"""
function ng(Δᵢ,χ)

    # Constants
    A₀ = 30.465860763763
    a₁ = 89.55438805808153
    a₂ = -130.792822369483
    a₃ = 42.02450507117405
    b₁ = -0.02929128383850
    b₂ = 1.0263250730647101E-5
    b₃ = -1.031921244571857E-9
    λ = Δᵢ*(A₀ + ((a₃/χ + a₂)/χ + a₁)/χ + ((b₃*χ + b₂)*χ + b₁)*χ)
end
