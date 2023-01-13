# Online Kernel Based Regression (OKBR)

## Single step algorithm, 1D (for testing)

mutable struct OKBR_single{T<:AbstractRange}
    x_eval_points::T
    w::Array{Float64,1}
    M1::Array{Float64,1}
    M2::Array{Float64,1}
    mem::Float64
    kernel
end
function OKBR_single(x_eval_points, kernel)
    Nx = length(x_eval_points)
    OKBR_single(
        x_eval_points,
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        NaN,
        kernel
    )
end
function add_data!(okbr::OKBR_single, X_right)
    X_left = okbr.mem
    if !isnan(X_left)
        ΔX = X_right - X_left
        for (j_ind, x_eval) in enumerate(okbr.x_eval_points)
            K_weight = okbr.kernel(x_eval - X_left)
            if K_weight > 0.0
                w_old = okbr.w[j_ind]
                M1_old = okbr.M1[j_ind]

                setindex!(okbr.w, w_old + K_weight, j_ind)
                setindex!(
                    okbr.M1,
                    update_wmean(okbr.M1[j_ind], w_old, ΔX, K_weight),
                    j_ind
                )
                setindex!(
                    okbr.M2,
                    update_wvar(okbr.M2[j_ind], M1_old, w_old, ΔX, okbr.M1[j_ind], K_weight),
                    j_ind
                )
            end
        end
    end
    okbr.mem = X_right
    return nothing
end

## Multi step algorithm, 1D

mutable struct OKBR_multiple{T<:AbstractRange,F<:Function}
    x_eval_points::T
    tau_i::UnitRange{Int}
    w::Array{Float64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
    w_mem::Array{Float64,1}
    kernel::F
end
function OKBR_multiple(x_eval_points, tau_i, kernel)
    #TODO: generalize to array not starting at 1
    Nx = length(x_eval_points)
    τ_len = length(tau_i)
    mem = zeros(Float64, τ_len)
    mem .= NaN
    OKBR_multiple(
        x_eval_points,
        tau_i,
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        mem,
        zeros(Float64, τ_len),
        kernel
    )
end

# Scaled moments
M1τ(ohbr::OKBR_multiple, dt) = ohbr.M1 ./ (dt*ohbr.tau_i)
M2τ(ohbr::OKBR_multiple, dt) = ohbr.M2 ./ (dt*ohbr.tau_i)

# Add data to moments
function add_data!(okbr::OKBR_multiple, X_right)
    for (i_tau, X_left) in enumerate(okbr.mem)
        if !isnan(X_left)
            ΔX = X_right - X_left
            for (j_ind, x_eval) in enumerate(okbr.x_eval_points)
                K_weight = okbr.kernel(x_eval - X_left)
                if K_weight > 0.0
                    w_old = okbr.w[i_tau,j_ind]
                    M1_old = okbr.M1[i_tau, j_ind]

                    setindex!(okbr.w, w_old + K_weight, i_tau, j_ind)
                    setindex!(
                        okbr.M1,
                        update_wmean(okbr.M1[i_tau, j_ind], w_old, ΔX, K_weight),
                        i_tau, j_ind
                    )
                    setindex!(
                        okbr.M2,
                        update_wvar(okbr.M2[i_tau, j_ind], M1_old, w_old, ΔX, okbr.M1[i_tau, j_ind], K_weight),
                        i_tau, j_ind
                    )
                end
            end
        end
    end
    update_mem!(okbr.mem, X_right)
    return nothing
end
