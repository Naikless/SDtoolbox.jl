using SDtoolbox, Test

@testset "SDtoolbox.jl" begin
    include("cv.jl")
    @test err > 0.99

    include("znd.jl")
    @test err_T > 0.99 && err_therm > 0.99

    include("cell_size.jl")
    @test 9e-3 <= λ <= 10e-3
end
