import SDtoolbox: θₑ, Λ, χ, cell_size, CANTERA_MECH_FILETYPE

P₁ = 1e5
T₁ = 300
ϕ = 1

X₁ = "H2:$(42*ϕ), O2:21,N2:79"
mech = "gri30."*CANTERA_MECH_FILETYPE

λ = cell_size(T₁,P₁,X₁,mech)

println("P₁ = $(P₁*1e-5) bar, T₁ = $T₁ K, ϕ = $ϕ : λ = $(λ*1e3) mm")

# case from "The propagation and failure mechanism of gaseous 
# detonations: experiments in porous-walled tubes" by M. Radulescu
P₁ = 23.5e3
T₁ = 298.15
X₁ = "H2:2,O2:1"

θₑ_ = θₑ(T₁,P₁,X₁,mech,t_end=1e-6)
Λ_ = Λ(T₁,P₁,X₁,mech,t_end=1e-6)
χ_ = χ(T₁,P₁,X₁,mech,t_end=1e-6)
