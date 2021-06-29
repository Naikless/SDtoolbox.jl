module SDtoolbox

export zndsolve, cvsolve, cell_size, CJspeed, PostShock_eq, PostShock_fr

include("postshock.jl")
include("cv.jl")
include("znd.jl")

# submodules that are essentially standalone
using .Postshock
using .CV
using .ZND

# ensure all python packages are loaded at run time to make the package
# precompile safe.
using PyCall

const ct = PyNULL()

function __init__()
    # load cantera python package
    copy!(ct, pyimport_conda("cantera","cantera","cantera"))
end

include("cell_size.jl")

end
