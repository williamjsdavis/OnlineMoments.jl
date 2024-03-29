module OnlineMoments

# using Statistics: mean
using OnlineStats: AutoCov, fit!, autocov

# Kernels
export Epaneknikov, Boxcar

# Traditional moments
export HBR_moments_A, HBR_moments_B, HBR_moments_C, HBR_moments_C2
export KBR_moments_A, KBR_moments_A2

# Uncorrected moments
export HBRu_moments_A, HBRu_moments_B, HBRu_moments_C, HBRu_moments_C2
export KBRu_moments_A, KBRu_moments_A2

# Online moments
export OHBR_single, OHBR_multiple
export OKBR_single, OKBR_multiple

# Uncorrected moments
export OHBRu_single, OHBRu_multiple
export OKBRu_single, OKBRu_multiple

# Modulo moments
export HBR_moments_mod
export HBRu_moments_mod
export KBR_moments_mod
export KBRu_moments_mod

export OHBR_mod_single, OHBR_mod_multiple
export OHBRu_mod_single, OHBRu_mod_multiple
export OKBR_mod_single, OKBR_mod_multiple
export OKBRu_mod_single, OKBRu_mod_multiple

# Online Welford methods
export OHBR_welford_single
export M2

# Turbulence moments
export HBR_moments_turb_A, HBR_moments_turb_C
export OHBRu_turb_multiple, OKBRu_turb_multiple

# Autocorrelation
export offline_autocorr, online_autocorr
export AutoCov

# Functions
export add_data!

export M1τ, M2τ

include("utils.jl")
include("statistics.jl")
include("kernels.jl")
include("HBR.jl")
include("HBR_uncorrected.jl")
include("KBR.jl")
include("KBR_uncorrected.jl")
include("OHBR.jl")
include("OHBR_welford.jl")
include("OHBR_uncorrected.jl")
include("OKBR.jl")
include("OKBR_uncorrected.jl")
include("XKBR_mod.jl")
include("XKBR_mod_uncorrected.jl")
include("OXBR_turbulence.jl")
include("autocorr.jl")

end # module
