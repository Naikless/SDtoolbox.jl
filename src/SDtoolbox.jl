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

    # use xml mech files in tests for cantera before v2.5, otherwise yaml
    if VersionNumber(ct.__version__) < v"2.5"
        global CANTERA_MECH_FILETYPE = "xml"
    else
        global CANTERA_MECH_FILETYPE = "yaml"
    end
end

include("cell_size.jl")

end
