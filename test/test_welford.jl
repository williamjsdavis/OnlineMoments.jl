# Testing Welford streaming algorithms

## Online Histogram Based Regression, 1D (online methods)

@testset "OHBR (single, Welford)" begin
    X_stream = stream_data(X_small)

    ohbr_single = OHBR_single(x_edges)
    ohbr_welford_single = OHBR_welford_single(x_edges)
    @testset "Structs" begin
        @test ohbr_welford_single.edges == x_edges
        @test size(ohbr_welford_single.N) == (N_x,)
        @test size(ohbr_welford_single.M1) == (N_x,)
        @test size(ohbr_welford_single.S) == (N_x,)
        @test size(ohbr_welford_single.mem) == ()
    end

    for _ in 1:N_data
        X_data = X_stream()
        add_data!(ohbr_single, X_data)
        add_data!(ohbr_welford_single, X_data)
    end

    @testset "Moments" begin
        # The original streaming algorithm should be identical to algorithm C
        @test all(ohbr_single.M1 .== M1_ref_C[1,:])
        @test all(ohbr_single.M2 .== M2_ref_C[1,:])

        # Welford's algorithm streaming algorithm should be identical to algorithm C
        @test all(ohbr_welford_single.M1 .== M1_ref_C[1,:])
        @test all(M2(ohbr_welford_single) .== M2_ref_C[1,:])

        #  Welford's algorithm should be identical to the other streaming algorithm
        @test all(ohbr_welford_single.M1 .== ohbr_single.M1)
        @test all(M2(ohbr_welford_single) .== ohbr_single.M2)
    end
end