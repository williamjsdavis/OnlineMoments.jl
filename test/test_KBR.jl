# Testing original algorithms

## Kernel Based Regression, 1D (normal methods)

@testset "KBR moments" begin
    @testset "Size" begin
        @test size(M1_K_ref_A) == (N_tau, N_x)
        @test size(M2_K_ref_A) == (N_tau, N_x)
        @test size(M1_K_ref_A2) == (N_tau, N_x)
        @test size(M2_K_ref_A2) == (N_tau, N_x)
        @test size(M1_KE_ref_A) == (N_tau, N_x)
        @test size(M2_KE_ref_A) == (N_tau, N_x)
        @test size(M1_KE_ref_A2) == (N_tau, N_x)
        @test size(M2_KE_ref_A2) == (N_tau, N_x)
        @test size(M1_K_ref_mod) == (N_tau, N_x)
        @test size(M2_K_ref_mod) == (N_tau, N_x)
        @test size(M1_K_ref_mod_shift) == (N_tau, N_x)
        @test size(M2_K_ref_mod_shift) == (N_tau, N_x)
    end

    @testset "Values" begin
        # Algorithms A and A2 give identical results
        @test all(M1_K_ref_A .== M1_K_ref_A2)
        @test all(M2_K_ref_A .== M2_K_ref_A2)
        @test all(M1_KE_ref_A .== M1_KE_ref_A2)
        @test all(M2_KE_ref_A .== M2_KE_ref_A2)
        
        # Algorithm A gives almost the same results as (HBR) A
        #NOTE: I should highlight that this is a regression test...
        #NOTE: Make it its own section?
        @test all(M1_K_ref_A .≈ M1_ref_A)
        @test all(M2_K_ref_A .≈ M2_ref_A)

        # Algorithm A gives almost the same results as (HBR) C
        #NOTE: Same comment as last...
        @test all(M1_K_ref_A .≈ M1_ref_C)
        @test all(M2_K_ref_A .≈ M2_ref_C)

        # Modulo moments
        # Algorithms A and "mod" give identical results, for large mod period
        @test all(M1_K_ref_A .== M1_K_ref_mod)
        @test all(M2_K_ref_A .== M2_K_ref_mod)

        # Algorithm "mod" gives almost the same results as (HBR) C
        #NOTE: Same comment as above
        @test all(M1_K_ref_mod .≈ M1_ref_C)
        @test all(M2_K_ref_mod .≈ M2_ref_C)

        # Test for translational invariance,
        # for translation = k*period, k>>1
        @test all(M1_K_ref_mod .≈ M1_K_ref_mod_shift)
        @test all(M2_K_ref_mod .≈ M2_K_ref_mod_shift)
    end
end