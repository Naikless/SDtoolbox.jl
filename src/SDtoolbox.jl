module SDtoolbox

export zndsolve, cvsolve, cell_size

include("cv.jl")
include("znd.jl")

# submodules that are essentially standalone
using .CV
using .ZND

# ensure all python packages are loaded at run time to make the package
# precompile safe.
using PyCall

const ct = PyNULL()
const postshock = PyNULL()
const CJspeed = PyNULL()
const PostShock_fr = PyNULL()
const PostShock_eq = PyNULL()

function __init__()
    copy!(ct, pyimport_conda("cantera","cantera","cantera"))

    # adds package dir to python path
    pushfirst!(PyVector(pyimport("sys")."path"), pkgdir(SDtoolbox))
    copy!(postshock, pyimport("sdtoolbox.postshock"))
    copy!(CJspeed, postshock.CJspeed)
    copy!(PostShock_fr, postshock.PostShock_fr)
    copy!(PostShock_eq, postshock.PostShock_eq)
end

include("cell_size.jl")

end
