# Testing streaming algorithms

## Online Kernel Based Regression, 1D (online methods)

@testset "OKBR (single)" begin
    X_stream = stream_data(X_small)

    hinv = inv(h)
    kernel_scaled(x) = hinv*kernel_boxcar(hinv*x)

    #okbr_single = OKBR_single(x_centers, kernel_scaled)
    okbr_single = OKBR_single(x_centers, kernel_boxcar, h)
    @testset "Structs" begin
        @test okbr_single.x_eval_points == x_centers
        @test size(okbr_single.w) == (N_x,)
        @test size(okbr_single.M1) == (N_x,)
        @test size(okbr_single.M2) == (N_x,)
        @test size(okbr_single.mem) == ()
    end

    for _ in 1:N_data
        add_data!(okbr_single, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm and algorithm A should be almost the same
        @test all(okbr_single.M1 .≈ M1_K_ref_A[1,:])
        @test all(okbr_single.M2 .≈ M2_K_ref_A[1,:])

    end
end

@testset "OKBR (multiple)" begin
    X_stream = stream_data(X_small)

    hinv = inv(h)
    kernel_scaled(x) = hinv*kernel_boxcar(hinv*x)

    #kbr_multiple = OKBR_multiple(x_centers, tau_i_range, kernel_scaled)
    kbr_multiple = OKBR_multiple(x_centers, tau_i_range, kernel_boxcar, h)
    @testset "Structs" begin
        @test kbr_multiple.x_eval_points == x_centers
        @test size(kbr_multiple.w) == (N_tau, N_x)
        @test size(kbr_multiple.M1) == (N_tau, N_x)
        @test size(kbr_multiple.M2) == (N_tau, N_x)
        @test size(kbr_multiple.mem) == (N_tau,)
    end

    for _ in 1:N_data
        add_data!(kbr_multiple, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm and algorithm A should be almost the same
        @test all(kbr_multiple.M1 .≈ M1_K_ref_A)
        @test all(kbr_multiple.M2 .≈ M2_K_ref_A)
    end
end