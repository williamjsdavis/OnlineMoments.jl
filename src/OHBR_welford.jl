# Online Histogram Based Regression (OHBR)
# Using Welford's stable algorithm

## Single step algorithm, 1D (for testing)

mutable struct OHBR_welford_single{T<:AbstractRange}
    edges::T
    N::Array{Int64,1}
    M1::Array{Float64,1}
    S::Array{Float64,1}
    mem::Float64
end
function OHBR_welford_single(x_range)
    Nx = length(x_range) - 1
    OHBR_welford_single(
        x_range,
        zeros(Int, Nx),
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        NaN
    )
end

# Calculate second moment
M2(ohbr::OHBR_welford_single) = ohbr.S ./ ohbr.N

# Add data to moments
function add_data!(ohbr::OHBR_welford_single, X_right)
    X_left = ohbr.mem
    if in_range(ohbr.edges, X_left)
        ΔX = X_right - X_left
        i = find_bin(ohbr.edges, X_left)
        M1_old = ohbr.M1[i]

        setindex!(ohbr.N, ohbr.N[i] + 1, i)
        setindex!(
            ohbr.M1,
            update_mean(ohbr.M1[i], ΔX, ohbr.N[i]),
            i
        )
        setindex!(
            ohbr.S,
            update_var_welford(ohbr.S[i], ohbr.M1[i], M1_old, ΔX),
            i
        )
    end
    ohbr.mem = X_right
    return nothing
end