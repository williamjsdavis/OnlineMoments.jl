# Testing streaming algorithms for turbulence calculations

## Online Kernel Based Regression, 1D (online methods), uncorrected

@testset "OKBR (multiple, turbulence)" begin
    X_stream = stream_data(X_small)

    # Time-difference sampling
    tau2_range = 1:tau2
    N_tau2 = length(tau2_range)

    okbr_u_multiple = OKBRu_turb_multiple(x_centers, tau2_range, kernel_boxcar, h)
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

    # Comparing with offline methods
    tau_i = tau2 .- tau1_samples

    @testset "Moments" begin
        # This streaming algorithm and HBR algorithm A/C should be almost the same
        @test all(okbr_u_multiple.M1[tau_i,:] .≈ M1_tu_est_A)
        @test all(okbr_u_multiple.M2[tau_i,:] .≈ M2_tu_est_A)
        @test all(okbr_u_multiple.M1[tau_i,:] .≈ M1_tu_est_C)
        @test all(okbr_u_multiple.M2[tau_i,:] .≈ M2_tu_est_C)
    end
end
