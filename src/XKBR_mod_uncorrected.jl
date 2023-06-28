# Online Histogram/Kernel Based Regression (OHBR/OKBR)

## OHBR modulo moments, single step algorithm, 1D (for testing), uncorrected

mutable struct OHBRu_mod_single{T<:AbstractRange}
    edges::T
    period::Float64
    N::Array{Int64,1}
    M1::Array{Float64,1}
    M2::Array{Float64,1}
    mem::Float64
end
function OHBRu_mod_single(x_range, period)
    Nx = length(x_range) - 1
    OHBRu_mod_single(
        x_range,
        period,
        zeros(Int, Nx),
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        NaN
    )
end

# Scaled moments
M1τ(ohbr::OHBRu_mod_single, dt) = ohbr.M1 / dt
M2τ(ohbr::OHBRu_mod_single, dt) = ohbr.M2 / dt

# Add data to moments
function add_data!(ohbr::OHBRu_mod_single, X_right)
    X_left = ohbr.mem
    i = find_mod_bin(ohbr.edges, ohbr.period, X_left)
    if i != 0
        ΔX = X_right - X_left

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

## OHBR modulo moments, multi step algorithm, 1D

#TODO: generalize tau_i
mutable struct OHBRu_mod_multiple{T<:AbstractRange}
    edges::T
    tau_i::UnitRange{Int}
    period::Float64
    N::Array{Int64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
    bin_mem::Array{Int64,1}
end
function OHBRu_mod_multiple(x_range, tau_i, period)
    Nx = length(x_range) - 1
    τ_len = length(tau_i)
    mem = zeros(Float64, τ_len)
    mem .= NaN
    OHBRu_mod_multiple(
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
M1τ(ohbr::OHBRu_mod_multiple, dt) = ohbr.M1 ./ (dt*ohbr.tau_i)
M2τ(ohbr::OHBRu_mod_multiple, dt) = ohbr.M2 ./ (dt*ohbr.tau_i)

# Add data to moments
function add_data!(ohbr::OHBRu_mod_multiple, X_right)
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
    update_mem!(ohbr.bin_mem, find_mod_bin(ohbr.edges, ohbr.period, X_right))
    return nothing
end

## OKBR modulo moments, single step algorithm, 1D (for testing)

mutable struct OKBRu_mod_single{T<:AbstractRange,K<:Kernel}
    x_eval_points::T
    period::Float64
    w::Array{Float64,1}
    M1::Array{Float64,1}
    M2::Array{Float64,1}
    mem::Float64
    kernel::K
    hinv::Float64
end
function OKBRu_mod_single(x_eval_points, period, kernel, h::Float64)
    Nx = length(x_eval_points)
    hinv = inv(h)
    OKBRu_mod_single(
        x_eval_points,
        period,
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        NaN,
        kernel,
        hinv
    )
end

# Add data
function add_data!(okbr::OKBRu_mod_single, X_right)
    X_left = okbr.mem
    if !isnan(X_left)
        ΔX = X_right - X_left
        for (j_ind, x_eval) in enumerate(okbr.x_eval_points)
            Δx_mod = d_mod(x_eval - X_left, okbr.period)
            K_weight = apply_kernel(Δx_mod, okbr.kernel, okbr.hinv)
            if K_weight > 0.0
                w_old = okbr.w[j_ind]
                M1_old = okbr.M1[j_ind]

                setindex!(okbr.w, w_old + K_weight, j_ind)
                setindex!(
                    okbr.M1,
                    OnlineMoments.update_wmean(okbr.M1[j_ind], w_old, ΔX, K_weight),
                    j_ind
                )
                setindex!(
                    okbr.M2,
                    update_wss(okbr.M2[j_ind], w_old, ΔX, K_weight),
                    j_ind
                )
            end
        end
    end
    okbr.mem = X_right
    return nothing
end

## OHBR modulo moments, multi step algorithm, 1D

#TODO: generalize tau_i
mutable struct OKBRu_mod_multiple{T<:AbstractRange,K<:Kernel}
    x_eval_points::T
    tau_i::UnitRange{Int}
    period::Float64
    w::Array{Int64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
    w_mem::Array{Int64,1}
    kernel::K
    hinv::Float64
end
function OKBRu_mod_multiple(x_eval_points, tau_i, period, kernel, h::Float64)
    Nx = length(x_eval_points)
    τ_len = length(tau_i)
    mem = zeros(Float64, τ_len)
    mem .= NaN
    hinv = inv(h)
    OKBRu_mod_multiple(
        x_eval_points,
        tau_i,
        period,
        zeros(Int, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        mem,
        zeros(Int, τ_len),
        kernel,
        hinv
    )
end

# Scaled moments
M1τ(okbr::OKBRu_mod_multiple, dt) = okbr.M1 ./ (dt*okbr.tau_i)
M2τ(okbr::OKBRu_mod_multiple, dt) = okbr.M2 ./ (dt*okbr.tau_i)

# Add data to moments
function add_data!(okbr::OKBRu_mod_multiple, X_right)
    for (i_tau, X_left) in enumerate(okbr.mem)
        if !isnan(X_left)
            ΔX = X_right - X_left
            for (j_ind, x_eval) in enumerate(okbr.x_eval_points)
                Δx_mod = d_mod(x_eval - X_left, okbr.period)
                #K_weight = okbr.kernel(x_eval - X_left)
                K_weight = apply_kernel(Δx_mod, okbr.kernel, okbr.hinv)
                if K_weight > 0.0
                    w_old = okbr.w[i_tau,j_ind]

                    setindex!(okbr.w, w_old + K_weight, i_tau, j_ind)
                    setindex!(
                        okbr.M1,
                        update_wmean(okbr.M1[i_tau, j_ind], w_old, ΔX, K_weight),
                        i_tau, j_ind
                    )
                    setindex!(
                        okbr.M2,
                        update_wss(okbr.M2[i_tau, j_ind], w_old, ΔX, K_weight),
                        i_tau, j_ind
                    )
                end
            end
        end
    end
    update_mem!(okbr.mem, X_right)
    return nothing
end
