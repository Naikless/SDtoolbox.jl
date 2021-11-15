import SDtoolbox: ct, cell_size
using Plots

P₁ = 1e5
T₁ = 300
ϕ = 1

X₁ = "H2:$(42*ϕ), O2:21,N2:79"
try
    global mech = "gri30.yaml"
catch PyError
    global mech = "gri30.xml"
end

λ = cell_size(T₁,P₁,X₁,mech)

println("P₁ = $(P₁*1e-5) bar, T₁ = $T₁ K, ϕ = $ϕ : λ = $(λ*1e3) mm")
