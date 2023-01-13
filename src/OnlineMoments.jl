module OnlineMoments

# using Statistics: mean

# Kernels
export Epaneknikov, Boxcar

# Traditional moments
export HBR_moments_A, HBR_moments_B, HBR_moments_C, HBR_moments_C2
export KBR_moments_A, KBR_moments_A2
export HBR_moments_A2 # Remove this

# Modulo moments
export HBR_moments_mod
export KBR_moments_mod

# Online moments
export OHBR_single, OHBR_multiple
export OKBR_single, OKBR_multiple
export OHBR_mod_single, OHBR_mod_multiple

export add_data!

export M1τ, M2τ

greet() = print("Hello World!")

include("utils.jl")
include("statistics.jl")
include("kernels.jl")
include("HBR.jl")
include("KBR.jl")
include("OHBR.jl")
include("OKBR.jl")

end # module
