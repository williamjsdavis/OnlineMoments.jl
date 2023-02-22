module OnlineMoments

# using Statistics: mean

# Kernels
export Epaneknikov, Boxcar

# Traditional moments
export HBR_moments_A, HBR_moments_B, HBR_moments_C, HBR_moments_C2
export KBR_moments_A, KBR_moments_A2
export HBR_moments_A2 # Remove this

# Online moments
export OHBR_single, OHBR_multiple
export OKBR_single, OKBR_multiple

# Modulo moments
export HBR_moments_mod
export KBR_moments_mod

export OHBR_mod_single, OHBR_mod_multiple
export OKBR_mod_single, OKBR_mod_multiple

# Online Welford methods
export OHBR_welford_single
export M2

# Functions
export add_data!

export M1τ, M2τ

include("utils.jl")
include("statistics.jl")
include("kernels.jl")
include("HBR.jl")
include("KBR.jl")
include("OHBR.jl")
include("OHBR_welford.jl")
include("OKBR.jl")
include("XKBR_mod.jl")

end # module
