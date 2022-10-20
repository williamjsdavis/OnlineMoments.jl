using OnlineMoments
using FileIO
using Test
include("./load_test_data.jl")

# Test settings
const tau_i_range = 1:1:6
const x_edges = LinRange(-0.0,0.4,5)

# Load test data
const X_small = load_small_data()
#X_small = X_small[1:10]

# Derived test variables
x_centers = 0.5*(x_edges[1:end-1] + x_edges[2:end])
N_tau = length(tau_i_range)
N_x = length(x_centers)

h = 0.5*step(x_centers)
kernel_boxcar = Boxcar()
kernel_epan = Epaneknikov()

## Histogram Based Regression, 1D (normal methods)

M1_ref_A, M2_ref_A = HBR_moments_A(X_small, tau_i_range, x_edges)
#M1_ref_A2, M2_ref_A2 = HBR_moments_A2(X_small, tau_i_range, x_edges) # Bad method
M1_ref_B, M2_ref_B = HBR_moments_B(X_small, tau_i_range, x_edges)
M1_ref_C, M2_ref_C = HBR_moments_C(X_small, tau_i_range, x_edges)
M1_ref_C2, M2_ref_C2 = HBR_moments_C2(X_small, tau_i_range, x_edges)

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
    end
end

## Online Histogram Based Regression, 1D (online methods)

@testset "Streaming data" begin
    X_sample = deepcopy(X_small)
    X_stream() = popfirst!(X_sample)

    @test X_stream() == X_small[1]
    @test X_stream() == X_small[2]
    @test X_stream() == X_small[3]
end

@testset "OHBR (single)" begin
    X_sample = deepcopy(X_small)
    X_stream() = popfirst!(X_sample)

    hbr_single = OHBR_single(x_edges)
    @testset "Structs" begin
        @test hbr_single.edges == x_edges
        @test size(hbr_single.N) == (N_x,)
        @test size(hbr_single.M1) == (N_x,)
        @test size(hbr_single.M2) == (N_x,)
        @test size(hbr_single.mem) == ()
    end

    for _ in 1:length(X_sample)
        add_data(hbr_single, X_stream())
    end
    @testset "Moments" begin
        # This streaming algorithm should be identical to algorithm C
        @test all(hbr_single.M1 .== M1_ref_C[1,:])
        @test all(hbr_single.M2 .== M2_ref_C[1,:])
    end
end

@testset "OHBR (multiple)" begin
    X_sample = deepcopy(X_small)
    X_stream() = popfirst!(X_sample)


    hbr_multiple = OHBR_multiple(x_edges, N_tau)
    @testset "Structs" begin
        @test hbr_multiple.edges == x_edges
        @test size(hbr_multiple.N) == (N_tau, N_x)
        @test size(hbr_multiple.M1) == (N_tau, N_x)
        @test size(hbr_multiple.M2) == (N_tau, N_x)
        @test size(hbr_multiple.mem) == (N_tau, )
    end

    for _ in 1:length(X_sample)
        add_data(hbr_multiple, X_stream())
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
@testset "Kernels" begin
    boxcar = OnlineMoments.Boxcar()
    @test boxcar(0.0) == 0.5
    @test boxcar(0.5) == 0.5
    @test boxcar(0.99) == 0.5
    @test boxcar(-0.99) == 0.5
    @test boxcar(-1.0) == 0.0
    @test boxcar(1.0) == 0.0

    epaneknikov = OnlineMoments.Epaneknikov()
    @test epaneknikov(0.0) ≈ 5*3*sqrt(5)/100
    @test epaneknikov(sqrt(5)) ≈ 0
    @test epaneknikov(sqrt(5)-0.01) > 0
    @test epaneknikov(sqrt(5)+0.01) == 0
    @test epaneknikov(sqrt(5)+1.0) == 0
    @test epaneknikov(-sqrt(5)+0.01) > 0
    @test epaneknikov(-sqrt(5)-0.01) == 0
    @test epaneknikov(-sqrt(5)-1.0) == 0
end


## Kernel Based Regression, 1D (normal methods)

# Boxcar kernel (can compare)
M1_K_ref_A, M2_K_ref_A = KBR_moments_A(X_small, tau_i_range, x_centers, h, kernel_boxcar)
M1_K_ref_A2, M2_K_ref_A2 = KBR_moments_A2(X_small, tau_i_range, x_centers, h, kernel_boxcar)

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
    end

    @testset "Values" begin
        # Algorithms A and A2 give identical results
        @test all(M1_K_ref_A .== M1_K_ref_A2)
        @test all(M2_K_ref_A .== M2_K_ref_A2)
        @test all(M1_KE_ref_A .== M1_KE_ref_A2)
        @test all(M2_KE_ref_A .== M2_KE_ref_A2)

        # Algorithm A gives almost the same results as (normal) C
        @test all(M1_K_ref_A .≈ M1_ref_C)
        @test all(M2_K_ref_A .≈ M2_ref_C)
    end
end
