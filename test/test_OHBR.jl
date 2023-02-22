# Testing streaming algorithms

## Online Histogram Based Regression, 1D (online methods)

@testset "OHBR (single)" begin
    X_stream = stream_data(X_small)

    ohbr_single = OHBR_single(x_edges)
    @testset "Structs" begin
        @test ohbr_single.edges == x_edges
        @test size(ohbr_single.N) == (N_x,)
        @test size(ohbr_single.M1) == (N_x,)
        @test size(ohbr_single.M2) == (N_x,)
        @test size(ohbr_single.mem) == ()
    end

    for _ in 1:N_data
        add_data!(ohbr_single, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm should be identical to algorithm C
        @test all(ohbr_single.M1 .== M1_ref_C[1,:])
        @test all(ohbr_single.M2 .== M2_ref_C[1,:])
    end
end

@testset "OHBR (single, mod)" begin
    X_stream = stream_data(X_small)

    ohbr_mod_single = OHBR_mod_single(x_edges, modulo_period_large)
    @testset "Structs" begin
        @test ohbr_mod_single.edges == x_edges
        @test ohbr_mod_single.period == modulo_period_large
        @test size(ohbr_mod_single.N) == (N_x,)
        @test size(ohbr_mod_single.M1) == (N_x,)
        @test size(ohbr_mod_single.M2) == (N_x,)
        @test size(ohbr_mod_single.mem) == ()
    end

    for _ in 1:N_data
        add_data!(ohbr_mod_single, X_stream())
    end
    @testset "Moments" begin
        # Algorithms C and streaming "mod" give identical results, for large mod period
        @test all(ohbr_mod_single.M1 .== M1_ref_C[1,:])
        @test all(ohbr_mod_single.M2 .== M2_ref_C[1,:])

        # Test for translational invariance,
        # for translation = k*period, k>>1
        @test all(ohbr_mod_single.M1 .≈ M1_ref_mod_shift[1,:])
        @test all(ohbr_mod_single.M2 .≈ M2_ref_mod_shift[1,:])
    end
end

@testset "OHBR (multiple)" begin
    X_stream = stream_data(X_small)

    ohbr_multiple = OHBR_multiple(x_edges, tau_i_range)
    @testset "Structs" begin
        @test ohbr_multiple.edges == x_edges
        @test ohbr_multiple.tau_i == tau_i_range
        @test size(ohbr_multiple.N) == (N_tau, N_x)
        @test size(ohbr_multiple.M1) == (N_tau, N_x)
        @test size(ohbr_multiple.M2) == (N_tau, N_x)
        @test size(ohbr_multiple.mem) == (N_tau,)
    end

    for _ in 1:N_data
        add_data!(ohbr_multiple, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm should be identical to algorithm C
        @test all(ohbr_multiple.M1 .== M1_ref_C)
        @test all(ohbr_multiple.M2 .== M2_ref_C)
    end
end

@testset "OHBR (multiple, mod)" begin
    X_stream = stream_data(X_small)

    ohbr_mod_multiple = OHBR_mod_multiple(x_edges, tau_i_range, modulo_period_large)
    @testset "Structs" begin
        @test ohbr_mod_multiple.edges == x_edges
        @test ohbr_mod_multiple.tau_i == tau_i_range
        @test ohbr_mod_multiple.period == modulo_period_large
        @test size(ohbr_mod_multiple.N) == (N_tau, N_x)
        @test size(ohbr_mod_multiple.M1) == (N_tau, N_x)
        @test size(ohbr_mod_multiple.M2) == (N_tau, N_x)
        @test size(ohbr_mod_multiple.mem) == (N_tau,)
    end

    for _ in 1:N_data
        add_data!(ohbr_mod_multiple, X_stream())
    end
    @testset "Moments" begin
        # Algorithms C and streaming "mod" give identical results, for large mod period
        @test all(ohbr_mod_multiple.M1 .== M1_ref_C)
        @test all(ohbr_mod_multiple.M2 .== M2_ref_C)

        # Test for translational invariance,
        # for translation = k*period, k>>1
        @test all(ohbr_mod_multiple.M1 .≈ M1_ref_mod_shift)
        @test all(ohbr_mod_multiple.M2 .≈ M2_ref_mod_shift)
    end
end