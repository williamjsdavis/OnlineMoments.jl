# Testing streaming algorithms

## Online Kernel Based Regression, 1D (online methods), uncorrected

@testset "OKBR (single, uncorrected)" begin
    X_stream = stream_data(X_small)

    okbr_u_single = OKBRu_single(x_centers, kernel_boxcar, h)
    @testset "Structs" begin
        @test okbr_u_single.x_eval_points == x_centers
        @test size(okbr_u_single.w) == (N_x,)
        @test size(okbr_u_single.M1) == (N_x,)
        @test size(okbr_u_single.M2) == (N_x,)
        @test size(okbr_u_single.mem) == ()
    end

    for _ in 1:N_data
        add_data!(okbr_u_single, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm and algorithm A should be almost the same
        @test all(okbr_u_single.M1 .≈ M1_Ku_ref_A[1,:])
        @test all(okbr_u_single.M2 .≈ M2_Ku_ref_A[1,:])

    end
end

@testset "OKBR (single, mod, uncorrected)" begin
    X_stream = stream_data(X_small)

    okbr_u_mod_single = OKBRu_mod_single(x_centers, modulo_period_large, kernel_boxcar, h)
    @testset "Structs" begin
        @test okbr_u_mod_single.x_eval_points == x_centers
        @test size(okbr_u_mod_single.w) == (N_x,)
        @test size(okbr_u_mod_single.M1) == (N_x,)
        @test size(okbr_u_mod_single.M2) == (N_x,)
        @test size(okbr_u_mod_single.mem) == ()
        @test okbr_u_mod_single.hinv == inv(h)
    end

    for _ in 1:N_data
        add_data!(okbr_u_mod_single, X_stream())
    end
    @testset "Moments" begin
        # Online algorithm almost matches non-online version
        @test all(okbr_u_mod_single.M1 .≈ M1_Ku_ref_mod[1,:])
        @test all(okbr_u_mod_single.M2 .≈ M2_Ku_ref_mod[1,:])

        # Algorithms A and streaming "mod" give identical results, for large mod period
        @test all(okbr_u_mod_single.M1 .≈ M1_Ku_ref_A[1,:])
        @test all(okbr_u_mod_single.M2 .≈ M2_Ku_ref_A[1,:])

        # Test for translational invariance,
        # for translation = k*period, k>>1
        @test all(okbr_u_mod_single.M1 .≈ M1_Ku_ref_mod_shift[1,:])
        @test all(okbr_u_mod_single.M2 .≈ M2_Ku_ref_mod_shift[1,:])
    end
end

@testset "OKBR (multiple, uncorrected)" begin
    X_stream = stream_data(X_small)

    okbr_u_multiple = OKBRu_multiple(x_centers, tau_i_range, kernel_boxcar, h)
    @testset "Structs" begin
        @test okbr_u_multiple.x_eval_points == x_centers
        @test size(okbr_u_multiple.w) == (N_tau, N_x)
        @test size(okbr_u_multiple.M1) == (N_tau, N_x)
        @test size(okbr_u_multiple.M2) == (N_tau, N_x)
        @test size(okbr_u_multiple.mem) == (N_tau,)
    end

    for _ in 1:N_data
        add_data!(okbr_u_multiple, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm and algorithm A should be almost the same
        @test all(okbr_u_multiple.M1 .≈ M1_Ku_ref_A)
        @test all(okbr_u_multiple.M2 .≈ M2_Ku_ref_A)
    end
end

@testset "OKBR (multiple, mod, uncorrected)" begin
    X_stream = stream_data(X_small)

    okbr_u_mod_multiple = OKBRu_mod_multiple(x_centers, tau_i_range, modulo_period_large, kernel_boxcar, h)
    @testset "Structs" begin
        @test okbr_u_mod_multiple.x_eval_points == x_centers
        @test okbr_u_mod_multiple.tau_i == tau_i_range
        @test okbr_u_mod_multiple.period == modulo_period_large
        @test size(okbr_u_mod_multiple.w) == (N_tau, N_x)
        @test size(okbr_u_mod_multiple.M1) == (N_tau, N_x)
        @test size(okbr_u_mod_multiple.M2) == (N_tau, N_x)
        @test size(okbr_u_mod_multiple.mem) == (N_tau,)
        @test okbr_u_mod_multiple.hinv == inv(h)
    end

    for _ in 1:N_data
        add_data!(okbr_u_mod_multiple, X_stream())
    end
    @testset "Moments" begin
        # Online algorithm almost matches non-online version
        @test all(okbr_u_mod_multiple.M1 .≈ M1_Ku_ref_mod)
        @test all(okbr_u_mod_multiple.M2 .≈ M2_Ku_ref_mod)

        # Algorithms A and streaming "mod" give identical results, for large mod period
        @test all(okbr_u_mod_multiple.M1 .≈ M1_Ku_ref_A)
        @test all(okbr_u_mod_multiple.M2 .≈ M2_Ku_ref_A)

        # Test for translational invariance,
        # for translation = k*period, k>>1
        @test all(okbr_u_mod_multiple.M1 .≈ M1_Ku_ref_mod_shift)
        @test all(okbr_u_mod_multiple.M2 .≈ M2_Ku_ref_mod_shift)
    end
end
