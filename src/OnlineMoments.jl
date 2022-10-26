module OnlineMoments

# using Statistics: mean

export Epaneknikov, Boxcar

export HBR_moments_A, HBR_moments_B, HBR_moments_C, HBR_moments_C2
export KBR_moments_A, KBR_moments_A2
export HBR_moments_A2 # Remove this

export OHBR_single, OHBR_multiple
export OKBR_single

export add_data

greet() = print("Hello World!")

include("utils.jl")
include("kernels.jl")
include("HBR.jl")
include("KBR.jl")
include("OHBR.jl")
include("OKBR.jl")

end # module
