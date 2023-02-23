# Testing autocorrelation

@testset "Autocorrelation" begin
    X_stream = stream_data(X_small)

    acf = AutoCov(maximum(tau_i_range))

    @testset "Structs" begin
        @test size(acf.m1) == (N_lag,)
    end

    for _ in 1:N_data
        add_data!(acf, X_stream())
    end

    acf_online = online_autocorr(acf)

    @testset "ACF values" begin
        # This streaming algorithm should be similar to the offline algorithm
        @test all(isapprox.(acf_offline, acf_online; rtol=1e-1))
    end
end