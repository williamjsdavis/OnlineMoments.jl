# Testing original algorithms for turbulence calculations

## Histogram Based Regression, 1D (normal methods)

@testset "HBR moments, turbulence" begin
    @testset "Size" begin
        @test size(M1_tu_est_A) == (N_tau1, N_x)
        @test size(M2_tu_est_A) == (N_tau1, N_x)
        @test size(M1_tu_est_C) == (N_tau1, N_x)
        @test size(M2_tu_est_C) == (N_tau1, N_x)
    end

    @testset "Values" begin
        # Algorithms A and C give almost the same results
        @test all(M1_tu_est_A .≈ M1_tu_est_C)
        @test all(M2_tu_est_A .≈ M2_tu_est_C)
    end
end