# Compare online algorithms

@testset "Comparing online algorithms" begin
    X_stream = stream_data(X_small)

    ohbr_single = OHBR_single(x_edges)
    ohbr_mod_single = OHBR_mod_single(x_edges, modulo_period_large)
    okbr_single = OKBR_single(x_centers, kernel_boxcar, h)
    okbr_mod_single = OKBR_mod_single(x_centers, modulo_period_large, kernel_boxcar, h)

    ohbr_multiple = OHBR_multiple(x_edges, tau_i_range)
    ohbr_mod_multiple = OHBR_mod_multiple(x_edges, tau_i_range, modulo_period_large)
    okbr_multiple = OKBR_multiple(x_centers, tau_i_range, kernel_boxcar, h)
    okbr_mod_multiple = OKBR_mod_multiple(x_centers, tau_i_range, modulo_period_large, kernel_boxcar, h)

    for _ in 1:N_data
        X_data = X_stream()
        add_data!(ohbr_single, X_data)
        add_data!(ohbr_mod_single, X_data)
        add_data!(okbr_single, X_data)
        add_data!(okbr_mod_single, X_data)

        add_data!(ohbr_multiple, X_data)
        add_data!(ohbr_mod_multiple, X_data)
        add_data!(okbr_multiple, X_data)
        add_data!(okbr_mod_multiple, X_data)
    end
    @testset "Moments" begin
        # Boxcar kernel is almost the same as OHBR, single
        @test all(ohbr_single.M1 .≈ okbr_single.M1)
        @test all(ohbr_single.M2 .≈ okbr_single.M2)

        # Boxcar kernel is almost the same as OHBR, multiple
        @test all(ohbr_multiple.M1 .≈ okbr_multiple.M1)
        @test all(ohbr_multiple.M2 .≈ okbr_multiple.M2)

        # Boxcar kernel is almost the same as OHBR, single (mod)
        @test all(ohbr_mod_single.M1 .≈ okbr_mod_single.M1)
        @test all(ohbr_mod_single.M2 .≈ okbr_mod_single.M2)

        # Boxcar kernel is almost the same as OHBR, multiple (mod)
        @test all(ohbr_mod_multiple.M1 .≈ okbr_mod_multiple.M1)
        @test all(ohbr_mod_multiple.M2 .≈ okbr_mod_multiple.M2)

        # First slice of OHBR multiple regresses to OHBR single
        @test all(ohbr_single.M1 .== ohbr_multiple.M1[1,:])
        @test all(ohbr_single.M2 .== ohbr_multiple.M2[1,:])

        # First slice of OKBR multiple regresses to OKBR single
        @test all(okbr_single.M1 .== okbr_multiple.M1[1,:])
        @test all(okbr_single.M2 .== okbr_multiple.M2[1,:])

        # First slice of OHBR multiple regresses to OHBR single (mod)
        @test all(ohbr_mod_single.M1 .== ohbr_mod_multiple.M1[1,:])
        @test all(ohbr_mod_single.M2 .== ohbr_mod_multiple.M2[1,:])

        # First slice of OKBR multiple regresses to OKBR single (mod)
        @test all(okbr_mod_single.M1 .== okbr_mod_multiple.M1[1,:])
        @test all(okbr_mod_single.M2 .== okbr_mod_multiple.M2[1,:])
    end
end