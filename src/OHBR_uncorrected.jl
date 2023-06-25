# Online Histogram Based Regression (OHBR), uncorrected second moment

## Single step algorithm, 1D (for testing)

mutable struct OHBRu_single{T<:AbstractRange}
    edges::T
    N::Array{Int64,1}
    M1::Array{Float64,1}
    M2::Array{Float64,1}
    mem::Float64
end
function OHBRu_single(x_range)
    Nx = length(x_range) - 1
    OHBRu_single(
        x_range,
        zeros(Int, Nx),
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        NaN
    )
end

# Scaled moments
M1τ(ohbr::OHBRu_single, dt) = ohbr.M1 / dt
M2τ(ohbr::OHBRu_single, dt) = ohbr.M2 / dt

# Add data to moments
function add_data!(ohbr::OHBRu_single, X_right)
    X_left = ohbr.mem
    if in_range(ohbr.edges, X_left)
        ΔX = X_right - X_left
        i = find_bin(ohbr.edges, X_left)

        setindex!(ohbr.N, ohbr.N[i] + 1, i)
        setindex!(
            ohbr.M1,
            update_mean(ohbr.M1[i], ΔX, ohbr.N[i]),
            i
        )
        setindex!(
            ohbr.M2,
            update_ss(ohbr.M2[i], ΔX, ohbr.N[i]),
            i
        )
    end
    ohbr.mem = X_right
    return nothing
end
#TODO: Make this return nothing and be consistant

## Multi step algorithm, 1D

#TODO: generalize tau_i
mutable struct OHBRu_multiple{T<:AbstractRange}
    edges::T
    tau_i::UnitRange{Int}
    N::Array{Int64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
    bin_mem::Array{Int64,1}
end
function OHBRu_multiple(x_range, tau_i)
    Nx = length(x_range) - 1
    τ_len = length(tau_i)
    mem = zeros(Float64, τ_len)
    mem .= NaN
    OHBRu_multiple(
        x_range,
        tau_i,
        zeros(Int, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        mem,
        zeros(Int, τ_len)
    )
end

# Scaled moments
M1τ(ohbr::OHBRu_multiple, dt) = ohbr.M1 ./ (dt*ohbr.tau_i)
M2τ(ohbr::OHBRu_multiple, dt) = ohbr.M2 ./ (dt*ohbr.tau_i)

#TODO: Remove mem_tmp?
# Add data to moments
function add_data!(ohbr::OHBRu_multiple, X_right)
    for (i_tau, j_bin) in enumerate(ohbr.bin_mem) if j_bin != 0
        X_left = ohbr.mem[i_tau]
        ΔX = X_right - X_left
        mem_tmp = X_left
        ohbr.mem[i_tau] = ohbr.M1[i_tau,j_bin] # Old mean

        setindex!(ohbr.N, ohbr.N[i_tau,j_bin] + 1, i_tau, j_bin)
        setindex!(
            ohbr.M1,
            update_mean(ohbr.M1[i_tau,j_bin], ΔX, ohbr.N[i_tau,j_bin]),
            i_tau, j_bin
        )
        setindex!(
            ohbr.M2,
            update_ss(ohbr.M2[i_tau,j_bin], ΔX, ohbr.N[i_tau,j_bin]),
            i_tau, j_bin
        )
        ohbr.mem[i_tau] = mem_tmp
    end end
    update_mem!(ohbr.mem, X_right)
    update_mem!(ohbr.bin_mem, get_bin(ohbr.edges, X_right))
    return nothing
end
