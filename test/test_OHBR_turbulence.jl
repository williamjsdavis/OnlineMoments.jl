# Testing streaming algorithms for turbulence calculations

## Online Histogram Based Regression, 1D (online methods)

@testset "OHBR (multiple, turbulence)" begin
    X_stream = stream_data(X_small)

    tau2_range = 1:tau2
    N_tau2 = length(tau2_range)

    ohbr_u_multiple = OHBRu_turb_multiple(x_edges, tau2_range)
    @testset "Structs" begin
        @test ohbr_u_multiple.edges == x_edges
        @test ohbr_u_multiple.tau_i == tau_i_range
        @test size(ohbr_u_multiple.N) == (N_tau2, N_x)
        @test size(ohbr_u_multiple.M1) == (N_tau2, N_x)
        @test size(ohbr_u_multiple.M2) == (N_tau2, N_x)
        @test size(ohbr_u_multiple.mem) == (N_tau2,)
    end

    for _ in 1:N_data
        add_data!(ohbr_u_multiple, X_stream())
    end

    # Comparing with offline methods
    tau_i = tau2 .- tau1_samples

    @testset "Moments" begin
        # This streaming algorithm should be almost identical to algorithm A/C
        @test all(ohbr_u_multiple.M1[tau_i,:] .≈ M1_tu_est_A)
        @test all(ohbr_u_multiple.M2[tau_i,:] .≈ M2_tu_est_A)
        @test all(ohbr_u_multiple.M1[tau_i,:] .≈ M1_tu_est_C)
        @test all(ohbr_u_multiple.M2[tau_i,:] .≈ M2_tu_est_C)
    end
end
