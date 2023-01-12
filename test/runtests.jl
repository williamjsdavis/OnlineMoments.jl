using OnlineMoments
using FileIO
using Test
include("./testutils.jl")

#TODO: Write some tests for skipping tau (e.g. [2,4,6])

# Load test data
const X_small = load_small_data()
const N_data = length(X_small)

# Test settings
const tau_i_range = 1:6
const x_edges = range(0.01, 0.41, length=5)

# Derived test variables
X_range = maximum(X_small) - minimum(X_small)
x_centers = 0.5*(x_edges[1:end-1] + x_edges[2:end])
N_tau = length(tau_i_range)
N_x = length(x_centers)

# Shifted test data (for modulo moments)
const modulo_period_large = 100*X_range
X_offset = 1000*modulo_period_large
const X_shift = X_small .- X_offset

# Kernel variables
h = 0.5*step(x_centers)
kernel_boxcar = Boxcar()
kernel_epan = Epaneknikov()

## Histogram Based Regression, 1D (normal methods)

M1_ref_A, M2_ref_A = HBR_moments_A(X_small, tau_i_range, x_edges)
#M1_ref_A2, M2_ref_A2 = HBR_moments_A2(X_small, tau_i_range, x_edges) # Bad method
M1_ref_B, M2_ref_B = HBR_moments_B(X_small, tau_i_range, x_edges)
M1_ref_C, M2_ref_C = HBR_moments_C(X_small, tau_i_range, x_edges)
M1_ref_C2, M2_ref_C2 = HBR_moments_C2(X_small, tau_i_range, x_edges)

# Modulo moments
M1_ref_mod, M2_ref_mod =
    HBR_moments_mod(X_small, tau_i_range, x_edges, modulo_period_large)
M1_ref_mod_shift, M2_ref_mod_shift =
    HBR_moments_mod(X_shift, tau_i_range, x_edges, modulo_period_large)

@testset "HBR moments" begin
    @testset "Size" begin
        @test size(M1_ref_A) == (N_tau, N_x)
        @test size(M2_ref_A) == (N_tau, N_x)
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
    end

    @testset "Values" begin
        # Algorithms A and B give different results
        @test !all(M1_ref_A .≈ M1_ref_B)
        @test !all(M2_ref_A .≈ M2_ref_B)

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
    end
end

## Online Histogram Based Regression, 1D (online methods)

@testset "Streaming data" begin
    X_stream = stream_data(X_small)

    @test X_small[1] == X_stream()
    @test X_small[2] == X_stream()
    @test X_small[3] == X_stream()
end

@testset "OHBR (single)" begin
    X_stream = stream_data(X_small)

    hbr_single = OHBR_single(x_edges)
    @testset "Structs" begin
        @test hbr_single.edges == x_edges
        @test size(hbr_single.N) == (N_x,)
        @test size(hbr_single.M1) == (N_x,)
        @test size(hbr_single.M2) == (N_x,)
        @test size(hbr_single.mem) == ()
    end

    for _ in 1:N_data
        add_data!(hbr_single, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm should be identical to algorithm C
        @test all(hbr_single.M1 .== M1_ref_C[1,:])
        @test all(hbr_single.M2 .== M2_ref_C[1,:])
    end
end

@testset "OHBR (multiple)" begin
    X_stream = stream_data(X_small)

    hbr_multiple = OHBR_multiple(x_edges, tau_i_range)
    @testset "Structs" begin
        @test hbr_multiple.edges == x_edges
        @test size(hbr_multiple.N) == (N_tau, N_x)
        @test size(hbr_multiple.M1) == (N_tau, N_x)
        @test size(hbr_multiple.M2) == (N_tau, N_x)
        @test size(hbr_multiple.mem) == (N_tau,)
    end

    for _ in 1:N_data
        add_data!(hbr_multiple, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm should be identical to algorithm C
        @test all(hbr_multiple.M1[1,:] .== M1_ref_C[1,:])
        @test all(hbr_multiple.M2[1,:] .== M2_ref_C[1,:])
        @test all(hbr_multiple.M1 .== M1_ref_C)
        @test all(hbr_multiple.M2 .== M2_ref_C)
    end
end

#TODO: Was I doing something here???
#M1_ref_A, M2_ref_A = HBR_moments_A(X_small, tau_i_range, x_edges)

## Kernels
include("./testkernels.jl")

## Kernel Based Regression, 1D (normal methods)

# Boxcar kernel (can compare)
M1_K_ref_A, M2_K_ref_A = KBR_moments_A(X_small, tau_i_range, x_centers, h, kernel_boxcar)
M1_K_ref_A2, M2_K_ref_A2 = KBR_moments_A2(X_small, tau_i_range, x_centers, h, kernel_boxcar)

# Boxcar kernel: modulo moments
M1_K_ref_mod, M2_K_ref_mod =
    KBR_moments_mod(X_small, tau_i_range, x_centers, h, kernel_boxcar, modulo_period_large)
M1_K_ref_mod_shift, M2_K_ref_mod_shift =
    KBR_moments_mod(X_shift, tau_i_range, x_centers, h, kernel_boxcar, modulo_period_large)

# Epaneknikov kernel (cannot validate)
M1_KE_ref_A, M2_KE_ref_A = KBR_moments_A(X_small, tau_i_range, x_centers, h, kernel_epan)
M1_KE_ref_A2, M2_KE_ref_A2 = KBR_moments_A2(X_small, tau_i_range, x_centers, h, kernel_epan)

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

        # Algorithm A gives almost the same results as (HBR) C
        #NOTE: I should highlight that this is a regression test...
        #NOTE: Make it its own section?
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

## Online Kernel Based Regression, 1D (online methods)

@testset "OKBR (single)" begin
    X_stream = stream_data(X_small)

    hinv = inv(h)
    kernel_scaled(x) = hinv*kernel_boxcar(hinv*x)

    kbr_single = OKBR_single(x_centers, kernel_scaled)
    @testset "Structs" begin
        @test kbr_single.x_eval_points == x_centers
        @test size(kbr_single.w) == (N_x,)
        @test size(kbr_single.M1) == (N_x,)
        @test size(kbr_single.M2) == (N_x,)
        @test size(kbr_single.mem) == ()
    end

    for _ in 1:N_data
        add_data!(kbr_single, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm and algorithm A should be almost the same
        @test all(kbr_single.M1 .≈ M1_K_ref_A[1,:])
        @test all(kbr_single.M2 .≈ M2_K_ref_A[1,:])

    end
end

@testset "OKBR (multiple)" begin
    X_stream = stream_data(X_small)

    hinv = inv(h)
    kernel_scaled(x) = hinv*kernel_boxcar(hinv*x)

    kbr_multiple = OKBR_multiple(x_centers, tau_i_range, kernel_scaled)
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
        # Consistancy with single step algorithm (move these to next section)
        #@test all(kbr_multiple.M1[1,:] .== kbr_single.M1)
        #@test all(kbr_multiple.M2[1,:] .== kbr_single.M2)

        # This streaming algorithm and algorithm A should be almost the same
        @test all(kbr_multiple.M1 .≈ M1_K_ref_A)
        @test all(kbr_multiple.M2 .≈ M2_K_ref_A)
    end
end

## Comparing online algorithms

@testset "Comparing online algorithms" begin
    X_stream = stream_data(X_small)

    hinv = inv(h)
    kernel_scaled(x) = hinv*kernel_boxcar(hinv*x)

    hbr_single = OHBR_single(x_edges)
    kbr_single = OKBR_single(x_centers, kernel_scaled)

    hbr_multiple = OHBR_multiple(x_edges, tau_i_range)

    for _ in 1:N_data
        X_data = X_stream()
        add_data!(hbr_single, X_data)
        add_data!(kbr_single, X_data)
        add_data!(hbr_multiple, X_data)
    end
    @testset "Moments" begin
        #
        @test all(hbr_single.M1 .≈ kbr_single.M1)
        @test all(hbr_single.M2 .≈ kbr_single.M2)

        #
        @test all(hbr_single.M1 .== hbr_multiple.M1[1,:])
        @test all(hbr_single.M2 .== hbr_multiple.M2[1,:])

        #
        @test all(kbr_single.M1 .≈ hbr_multiple.M1[1,:])
        @test all(kbr_single.M2 .≈ hbr_multiple.M2[1,:])
    end
end
