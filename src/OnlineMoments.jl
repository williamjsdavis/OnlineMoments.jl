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

# Online moments
export OHBR_single, OHBR_multiple
export OKBR_single, OKBR_multiple

# Modulo moments
export HBR_moments_mod
export HBRu_moments_mod
export KBR_moments_mod

export OHBR_mod_single, OHBR_mod_multiple
export OKBR_mod_single, OKBR_mod_multiple

# Online Welford methods
export OHBR_welford_single
export M2

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
include("OHBR.jl")
include("OHBR_welford.jl")
include("OKBR.jl")
include("XKBR_mod.jl")
include("autocorr.jl")

end # module
