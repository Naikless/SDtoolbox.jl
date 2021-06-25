using SDtoolbox
using Test

@testset "SDtoolbox.jl" begin
    include("cv.jl")
    include("znd.jl")
    incldue("cell_size.jl")
end
