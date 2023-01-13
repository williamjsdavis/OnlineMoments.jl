# Online Histogram Based Regression (OHBR)

## Single step algorithm, 1D (for testing)

mutable struct OHBR_single{T<:AbstractRange}
    edges::T
    N::Array{Int64,1}
    M1::Array{Float64,1}
    M2::Array{Float64,1}
    mem::Float64
end
function OHBR_single(x_range)
    Nx = length(x_range) - 1
    OHBR_single(
        x_range,
        zeros(Int, Nx),
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        NaN
    )
end

# Scaled moments
M1τ(ohbr::OHBR_single, dt) = ohbr.M1 / dt
M2τ(ohbr::OHBR_single, dt) = ohbr.M2 / dt

# Add data to moments
function add_data!(ohbr::OHBR_single, X_right)
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
            ohbr.M2,
            update_var(ohbr.M2[i], ohbr.M1[i], M1_old, ΔX, ohbr.N[i]),
            i
        )
    end
    ohbr.mem = X_right
    return nothing
end
#TODO: Make this return nothing and be consistant

## Multi step algorithm, 1D

#TODO: generalize tau_i
mutable struct OHBR_multiple{T<:AbstractRange}
    edges::T
    tau_i::UnitRange{Int}
    N::Array{Int64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
    bin_mem::Array{Int64,1}
end
function OHBR_multiple(x_range, tau_i)
    Nx = length(x_range) - 1
    τ_len = length(tau_i)
    mem = zeros(Float64, τ_len)
    mem .= NaN
    OHBR_multiple(
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
M1τ(ohbr::OHBR_multiple, dt) = ohbr.M1 ./ (dt*ohbr.tau_i)
M2τ(ohbr::OHBR_multiple, dt) = ohbr.M2 ./ (dt*ohbr.tau_i)

# Add data to moments
function add_data!(ohbr::OHBR_multiple, X_right)
    for (i_tau, j_bin) in enumerate(ohbr.bin_mem) if j_bin != 0
        X_left = ohbr.mem[i_tau]
        ΔX = X_right - X_left
        mem_tmp = X_left
        M1_old = ohbr.M1[i_tau,j_bin]
        ohbr.mem[i_tau] = ohbr.M1[i_tau,j_bin] # Old mean

        setindex!(ohbr.N, ohbr.N[i_tau,j_bin] + 1, i_tau, j_bin)
        setindex!(
            ohbr.M1,
            update_mean(M1_old, ΔX, ohbr.N[i_tau,j_bin]),
            i_tau, j_bin
        )
        setindex!(
            ohbr.M2,
            update_var(
                ohbr.M2[i_tau,j_bin],
                ohbr.M1[i_tau,j_bin],
                M1_old, ΔX,
                ohbr.N[i_tau,j_bin]
            ),
            i_tau, j_bin
        )
        ohbr.mem[i_tau] = mem_tmp
    end end
    update_mem!(ohbr.mem, X_right)
    update_mem!(ohbr.bin_mem, get_bin(ohbr.edges, X_right))
    return nothing
end

## Modulo moments, single step algorithm, 1D (for testing)

mutable struct OHBR_mod_single{T<:AbstractRange}
    edges::T
    period::Float64
    N::Array{Int64,1}
    M1::Array{Float64,1}
    M2::Array{Float64,1}
    mem::Float64
end
function OHBR_mod_single(x_range, period)
    Nx = length(x_range) - 1
    OHBR_mod_single(
        x_range,
        period,
        zeros(Int, Nx),
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        NaN
    )
end

# Scaled moments
M1τ(ohbr::OHBR_mod_single, dt) = ohbr.M1 / dt
M2τ(ohbr::OHBR_mod_single, dt) = ohbr.M2 / dt

# Add data to moments
function add_data!(ohbr::OHBR_mod_single, X_right)
    X_left = ohbr.mem
    i = find_mod_bin(ohbr.edges, ohbr.period, X_left)
    if i != 0
        ΔX = X_right - X_left
        M1_old = ohbr.M1[i]

        setindex!(ohbr.N, ohbr.N[i] + 1, i)
        setindex!(
            ohbr.M1,
            update_mean(ohbr.M1[i], ΔX, ohbr.N[i]),
            i
        )
        setindex!(
            ohbr.M2,
            update_var(ohbr.M2[i], ohbr.M1[i], M1_old, ΔX, ohbr.N[i]),
            i
        )
    end
    ohbr.mem = X_right
    return nothing
end

## Modulo moments, multi step algorithm, 1D

#TODO: generalize tau_i
mutable struct OHBR_mod_multiple{T<:AbstractRange}
    edges::T
    tau_i::UnitRange{Int}
    period::Float64
    N::Array{Int64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
    bin_mem::Array{Int64,1}
end
function OHBR_mod_multiple(x_range, tau_i, period)
    Nx = length(x_range) - 1
    τ_len = length(tau_i)
    mem = zeros(Float64, τ_len)
    mem .= NaN
    OHBR_mod_multiple(
        x_range,
        tau_i,
        period,
        zeros(Int, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        mem,
        zeros(Int, τ_len)
    )
end

# Scaled moments
M1τ(ohbr::OHBR_mod_multiple, dt) = ohbr.M1 ./ (dt*ohbr.tau_i)
M2τ(ohbr::OHBR_mod_multiple, dt) = ohbr.M2 ./ (dt*ohbr.tau_i)

# Add data to moments
function add_data!(ohbr::OHBR_mod_multiple, X_right)
    for (i_tau, j_bin) in enumerate(ohbr.bin_mem) if j_bin != 0
        X_left = ohbr.mem[i_tau]
        ΔX = X_right - X_left
        mem_tmp = X_left
        M1_old = ohbr.M1[i_tau,j_bin]
        ohbr.mem[i_tau] = ohbr.M1[i_tau,j_bin] # Old mean

        setindex!(ohbr.N, ohbr.N[i_tau,j_bin] + 1, i_tau, j_bin)
        setindex!(
            ohbr.M1,
            update_mean(M1_old, ΔX, ohbr.N[i_tau,j_bin]),
            i_tau, j_bin
        )
        setindex!(
            ohbr.M2,
            update_var(
                ohbr.M2[i_tau,j_bin],
                ohbr.M1[i_tau,j_bin],
                M1_old, ΔX,
                ohbr.N[i_tau,j_bin]
            ),
            i_tau, j_bin
        )
        ohbr.mem[i_tau] = mem_tmp
    end end
    update_mem!(ohbr.mem, X_right)
    update_mem!(ohbr.bin_mem, find_mod_bin(ohbr.edges, ohbr.period, X_right))
    return nothing
end
