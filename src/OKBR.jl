# Online Kernel Based Regression (OKBR)

## Updating statistics
function update_wmean(x_bar, w, x_new, w_new)
    return x_bar + (x_new - x_bar)*(w_new/(w + w_new))
end
function update_wvar(s2, x_bar, w, x_new, x_bar_new, w_new)
    return (s2*w + w_new*(x_new - x_bar)*(x_new - x_bar_new))/(w + w_new)
end

## Single step algorithm, 1D (for testing)

mutable struct OKBR_single
    x_eval_points::LinRange{Float64}
    N::Array{Float64,1} #NOTE: Channge N to w
    M1::Array{Float64,1}
    M2::Array{Float64,1}
    mem::Float64
    kernel
end
function OKBR_single(x_eval_points::LinRange, kernel)
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
function add_data(OKBR::OKBR_single, X_new)
    if !isnan(OKBR.mem)
        ΔX = X_new - OKBR.mem
        for (j_ind,j_xeval) in enumerate(OKBR.x_eval_points)
            K_weight = OKBR.kernel(j_xeval - OKBR.mem)
            if K_weight > 0.0
                OKBR.mem = OKBR.N[j_ind] # Old weight
                setindex!(OKBR.N, OKBR.N[j_ind] + K_weight, j_ind)
                tmp = OKBR.M1[j_ind]
                setindex!(
                    OKBR.M1,
                    update_wmean(OKBR.M1[j_ind], OKBR.mem, ΔX, K_weight),
                    j_ind
                )
                setindex!(
                    OKBR.M2,
                    update_wvar(OKBR.M2[j_ind], tmp, OKBR.mem, ΔX, OKBR.M1[j_ind], K_weight),
                    j_ind
                )
            end
        end
    end
    OKBR.mem = X_new
end

## Multi step algorithm, 1D

mutable struct OKBR_multiple
    x_eval_points::LinRange{Float64}
    w::Array{Float64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
    w_mem::Array{Float64,1}
    kernel
end
function OKBR_multiple(x_eval_points::LinRange, τ_len::Integer, kernel)
    Nx = length(x_eval_points)
    mem = zeros(Float64, τ_len)
    mem .= NaN
    OKBR_multiple(
        x_eval_points,
        zeros(Int, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        mem,
        zeros(Float64, τ_len),
        kernel
    )
end
function add_data(OKBR::OKBR_multiple, x_data)
    for (i_tau, x_left) in enumerate(OKBR.mem)
        if !isnan(x_left)
            ΔX = x_data - x_left
            for (j_ind,j_xeval) in enumerate(OKBR.x_eval_points)
                K_weight = OKBR.kernel(j_xeval - x_left)
                if K_weight > 0.0
                    mem_tmp = OKBR.mem[i_tau]
                    OKBR.mem[i_tau] = OKBR.w[i_tau,j_ind] # Old weight
                    setindex!(OKBR.w, OKBR.w[i_tau,j_ind] + K_weight, i_tau, j_ind)
                    tmp = OKBR.M1[i_tau, j_ind]
                    setindex!(
                        OKBR.M1,
                        update_wmean(OKBR.M1[i_tau, j_ind], OKBR.mem[i_tau], ΔX, K_weight),
                        i_tau, j_ind
                    )
                    setindex!(
                        OKBR.M2,
                        update_wvar(OKBR.M2[i_tau, j_ind], tmp, OKBR.mem[i_tau], ΔX, OKBR.M1[i_tau, j_ind], K_weight),
                        i_tau, j_ind
                    )
                    OKBR.mem[i_tau] = mem_tmp
                end
            end
        end
    end
    update_mem!(OKBR.mem, x_data)
end