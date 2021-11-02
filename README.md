# SDtoolbox

This package aims to (at least partly) reimplement [CalTech's Shock and Detonation Toolbox](https://shepherd.caltech.edu/EDL/PublicResources/sdt/) ("sdtoolbox") in Julia.

Focus is put on performance critical parts such as ODE integration. With this in mind, the original SciPy libraries are replaced with Julia's [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). Besides calculation speed improvements, this should also provide more stability for more complex kinetic mechanisms.

# Usage
## Define initial conditions and select kinetic mechanism:
```julia
P₁ = 1e5
T₁ = 300
X₁ = "H2:42, O2:21,N2:79"

# uses the GRI 3.0 kinetic mechanism
mech = "gri30.xml"
```

## Construct Cantera gas objects
SDtoolbox exports direct access to Cantera Python Package via [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) as `ct`.
```julia
gas₁ = ct.Solution(mech)
gas₁.TPX = T₁,P₁,X₁
```
## Calculate CJ detonation velocity
```julia
U₁ = CJspeed(P₁,T₁,X₁,mech)
```
`CJspeed` accepts an optional keyword `method` which defines the algorithm to calculate the Chapman-Jouguet velocity for the given mixture. The two methods proposed by the original sdtoolbox (`umin`,`aeq`) are complemented by the algorithm used by the well-established [NASA CEARUN](https://cearun.grc.nasa.gov/). So far, tests reveiled the latter  to be the most stable and efficient, which is why `CEA` is currently the default setting.

## Calculate post shock states:
Uses shock relations to calculate a frozen (i.e. no change in composition) state. The equilibrium state is then reached by forcing chemical equilibrium.
```julia
gas_fr = PostShock_fr(U₁, P₁, T₁, X₁, mech)
gas_eq = PostShock_eq(U₁, P₁, T₁, X₁, mech)
```


```
out = zndsolve(gas,gas₁,U₁,advanced_output=true)
```
## `cvsolve(gas)`
## `cell_size(T₁,P₁,X₁,mech)`

# Future improvements

Since code transfer to Julia is done gradually and there is currently no direct Julia interface for [Cantera](https://cantera.org/), this package heavily relies on [PyCall](https://github.com/JuliaPy/PyCall.jl). This results in all of the cantera calls to be the most obvious remaining performance bottlenecks, although speedups compared to the pure Python interpretation of the sdtoolbox can still be quite significant. 

Possible solutions for this could be:
 - Wait until there is a native Julia interface for Cantera
 - Use the [Arrhenius.jl package](https://juliahub.com/ui/Packages/Arrhenius/ymlFQ/0.1.1)
 - Call Cantera C++ subroutines directly

## Additional to dos
 - use multiple dispatch on cell_size to calculate CJ speed only when necessary.
 - move try-catch from cell_size to CJspeed.
