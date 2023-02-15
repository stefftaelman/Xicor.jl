using Xicor
using Test

@testset "Xicor.jl" begin
    X = [1, 2, 3, 5]
    Y = [1.0, 3.0, 8.0, 10.0]

    @test ξ(X, Y) > 0 
    @test ξ(X, Y) ≈ ξ(X, Y, noties=true)
end
