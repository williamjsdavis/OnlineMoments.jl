# Testing original algorithms

## Histogram Based Regression, 1D (normal methods), uncorrected

@testset "HBR moments, uncorrected" begin
    @testset "Size" begin
        @test size(M1_u_ref_A) == (N_tau, N_x)
        @test size(M2_u_ref_A) == (N_tau, N_x)
        @test size(M1_u_ref_B) == (N_tau, N_x)
        @test size(M2_u_ref_B) == (N_tau, N_x)
        #=
        @test size(M1_ref_B) == (N_tau, N_x)
        @test size(M2_ref_B) == (N_tau, N_x)
        @test size(M1_ref_C) == (N_tau, N_x)
        @test size(M2_ref_C) == (N_tau, N_x)
        @test size(M1_ref_C2) == (N_tau, N_x)
        @test size(M2_ref_C2) == (N_tau, N_x)
        @test size(M1_ref_mod) == (N_tau, N_x)
        @test size(M2_ref_mod) == (N_tau, N_x)
        @test size(M2_ref_mod_shift) == (N_tau, N_x)
        @test size(M2_ref_mod_shift) == (N_tau, N_x)
        =#
    end

    @testset "Values" begin
        # Algorithms A and B give different results
        @test !all(M1_u_ref_A .≈ M1_u_ref_B)
        @test !all(M2_u_ref_A .≈ M2_u_ref_B)

        #=
        # Algorithms A and C give almost the same results
        @test all(M1_ref_A .≈ M1_ref_C)
        @test all(M2_ref_A .≈ M2_ref_C)

        # Algorithms C and C2 give identical results
        @test all(M1_ref_C .== M1_ref_C2)
        @test all(M2_ref_C .== M2_ref_C2)

        # Modulo moments
        # Algorithms C and "mod" give identical results, for large mod period
        @test all(M1_ref_C .== M1_ref_mod)
        @test all(M2_ref_C .== M2_ref_mod)

        # Test for translational invariance,
        # for translation = k*period, k>>1
        @test all(M1_ref_mod .≈ M1_ref_mod_shift)
        @test all(M2_ref_mod .≈ M2_ref_mod_shift)
        =#
    end
end