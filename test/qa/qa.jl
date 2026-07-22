using SciMLTesting, GlobalDiffEq, Test
using JET
using OrdinaryDiffEq, OrdinaryDiffEqSSPRK

run_qa(GlobalDiffEq)

@testset "GlobalRichardson static analysis" begin
    @test_opt GlobalRichardson(SSPRK33())
    @test GlobalRichardson{typeof(SSPRK33())} <: GlobalDiffEq.GlobalDiffEqAlgorithm
    @test GlobalRichardson(Tsit5()) isa GlobalDiffEq.GlobalDiffEqAlgorithm
end
