# SDtoolbox

This package aims to (at least partly) reimplement [CalTech's Shock and Detonation Toolbox](https://shepherd.caltech.edu/EDL/PublicResources/sdt/) in Julia. 

Focus is put on performance critical parts such as ODE integration. For this the original SciPy libraries are replaced with Julia's [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). 

Since code transfer to Julia is done gradually and there is currently no direct Julia interface for [Cantera](https://cantera.org/), this package heavily relies on [PyCall](https://github.com/JuliaPy/PyCall.jl).
