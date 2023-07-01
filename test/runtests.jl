using OnlineMoments
using FileIO
using Test
include("./utils_for_tests.jl")

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
N_lag = N_tau + 1
N_x = length(x_centers)

# Shifted test data (for modulo moments)
const modulo_period_large = 100*X_range
X_offset = 1000*modulo_period_large
const X_shift = X_small .- X_offset

# Kernel variables
h = 0.5*step(x_centers)
kernel_boxcar = Boxcar()
kernel_epan = Epaneknikov()

## Initial tests
include("./test_utils.jl")
include("./test_kernels.jl")

## Histogram Based Regression, 1D (normal methods)

# Data, original
M1_ref_A, M2_ref_A = HBR_moments_A(X_small, tau_i_range, x_edges)
M1_ref_B, M2_ref_B = HBR_moments_B(X_small, tau_i_range, x_edges)
M1_ref_C, M2_ref_C = HBR_moments_C(X_small, tau_i_range, x_edges)
M1_ref_C2, M2_ref_C2 = HBR_moments_C2(X_small, tau_i_range, x_edges)

# Data, original, uncorrected second moment
M1_u_ref_A, M2_u_ref_A = HBRu_moments_A(X_small, tau_i_range, x_edges)
M1_u_ref_B, M2_u_ref_B = HBRu_moments_B(X_small, tau_i_range, x_edges)
M1_u_ref_C, M2_u_ref_C = HBRu_moments_C(X_small, tau_i_range, x_edges)
M1_u_ref_C2, M2_u_ref_C2 = HBRu_moments_C2(X_small, tau_i_range, x_edges)

# Data, original, modulo
M1_ref_mod, M2_ref_mod =
    HBR_moments_mod(X_small, tau_i_range, x_edges, modulo_period_large)
M1_ref_mod_shift, M2_ref_mod_shift =
    HBR_moments_mod(X_shift, tau_i_range, x_edges, modulo_period_large)

# Data, original, modulo, uncorrected second moment
M1_u_ref_mod, M2_u_ref_mod =
    HBRu_moments_mod(X_small, tau_i_range, x_edges, modulo_period_large)
M1_u_ref_mod_shift, M2_u_ref_mod_shift =
    HBRu_moments_mod(X_shift, tau_i_range, x_edges, modulo_period_large)

# Run tests
include("./test_HBR.jl")
include("./test_HBR_uncorrected.jl")

## Online Histogram Based Regression, 1D (online methods)

# Original tests
include("./test_OHBR.jl")
include("./test_OHBR_uncorrected.jl")

# Welford tests
include("./test_welford.jl")

## Kernel Based Regression, 1D (normal methods)

# Data, boxcar kernel (can compare)
M1_K_ref_A, M2_K_ref_A = KBR_moments_A(X_small, tau_i_range, x_centers, h, kernel_boxcar)
M1_K_ref_A2, M2_K_ref_A2 = KBR_moments_A2(X_small, tau_i_range, x_centers, h, kernel_boxcar)

# Data, boxcar kernel (can compare), uncorrected second moment
M1_Ku_ref_A, M2_Ku_ref_A = KBRu_moments_A(X_small, tau_i_range, x_centers, h, kernel_boxcar)
M1_Ku_ref_A2, M2_Ku_ref_A2 = KBRu_moments_A2(X_small, tau_i_range, x_centers, h, kernel_boxcar)

# Data, boxcar kernel, modulo
M1_K_ref_mod, M2_K_ref_mod =
    KBR_moments_mod(X_small, tau_i_range, x_centers, h, kernel_boxcar, modulo_period_large)
M1_K_ref_mod_shift, M2_K_ref_mod_shift =
    KBR_moments_mod(X_shift, tau_i_range, x_centers, h, kernel_boxcar, modulo_period_large)

# Data, boxcar kernel, modulo, uncorrected second moment
M1_Ku_ref_mod, M2_Ku_ref_mod =
    KBRu_moments_mod(X_small, tau_i_range, x_centers, h, kernel_boxcar, modulo_period_large)
M1_Ku_ref_mod_shift, M2_Ku_ref_mod_shift =
    KBRu_moments_mod(X_shift, tau_i_range, x_centers, h, kernel_boxcar, modulo_period_large)

# Data epaneknikov kernel (cannot validate)
M1_KE_ref_A, M2_KE_ref_A = KBR_moments_A(X_small, tau_i_range, x_centers, h, kernel_epan)
M1_KE_ref_A2, M2_KE_ref_A2 = KBR_moments_A2(X_small, tau_i_range, x_centers, h, kernel_epan)

# Data epaneknikov kernel (cannot validate), uncorrected second moment
M1_KEu_ref_A, M2_KEu_ref_A = KBRu_moments_A(X_small, tau_i_range, x_centers, h, kernel_epan)
M1_KEu_ref_A2, M2_KEu_ref_A2 = KBRu_moments_A2(X_small, tau_i_range, x_centers, h, kernel_epan)

# Run tests
include("./test_KBR.jl")
include("./test_KBR_uncorrected.jl")

## Online Kernel Based Regression, 1D (online methods)

include("./test_OKBR.jl")
include("./test_OKBR_uncorrected.jl")

## Comparing online algorithms

include("./compare_online.jl")

## Testing turbulence moments
tau1_samples = [4,3,2]
N_tau1 = length(tau1_samples)
tau2 = 6

## Histogram methods
M1_tu_est_A, M2_tu_est_A = HBR_moments_turb_A(X_small, tau1_samples, tau2, x_edges)
M1_tu_est_C, M2_tu_est_C = HBR_moments_turb_C(X_small, tau1_samples, tau2, x_edges)

## Online Histogram methods

include("./test_OHBR_turbulence.jl")

## Online Kernel methods

include("./test_OKBR_turbulence.jl")

## Testing autocorrelation

# Offline method
acf_offline = offline_autocorr(X_small, maximum(tau_i_range))

include("./test_autocorr.jl")
