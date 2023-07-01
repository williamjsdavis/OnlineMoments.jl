# (Online) Histogram/Kernel Based Regression (OKBR) for turbulence moments
#NOTES:
# - No single variants for online methods
# - No offline kernel methods


## Offline algorithms

## Multi step algorithm, 1D

""" Histogram Based Regression moments, Algorithm A """
function HBR_moments_turb_A(X, tau1_samples, tau2, edge_vector)
    nτ = length(tau1_samples)
    nx = length(edge_vector) - 1
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)
    
    incr2 = X[(tau2+1):end] - X[1:end-tau2]
    bin_index = map(incr2) do y
        in_range(edge_vector,y) ? find_bin(edge_vector,y) : 0
    end 
    for (i,tau1) in enumerate(tau1_samples)
        incr1 = X[(tau1+1):(end-(tau2-tau1))] - X[1:(end-tau2)]
        
        for j in 1:(length(edge_vector)-1)
            iX = findall(bin_index[1:end] .== j)
            ΔX = incr1[iX] .- incr2[iX]
            M1[i,j] = mean(ΔX)
            residuals = ΔX
            M2[i,j] = mean(residuals.^2)
        end
    end
    return M1, M2
end

""" Histogram Based Regression moments, Algorithm C """
function HBR_moments_turb_C(X, tau1_samples, tau2, edge_vector)
    nτ = length(tau1_samples)
    nx = length(edge_vector) - 1
    nX = length(X)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)

    for (i_left, X_left) in enumerate(X[1:end-tau2])
        incr2 = X[i_left+tau2] - X_left
        if in_range(edge_vector, incr2)
            k = find_bin(edge_vector, incr2)
            for (j,tau1) in enumerate(tau1_samples)
                ii = i_left + j
                if ii <= nX
                    incr1 = X[i_left+tau1] - X_left
                    ΔX = incr1 - incr2
                    setindex!(N, N[j,k] + 1, j, k)
                    setindex!(
                        M1,
                        update_mean(M1[j, k], ΔX, N[j, k]),
                        j, k
                    )
                    setindex!(
                        M2,
                        update_ss(M2[j, k], ΔX, N[j, k]),
                        j, k
                    )
                end
            end
        end
    end

    return M1, M2
end

## Online algorithms

## Online Hisogram algrithms

## Multi step algorithm, 1D

mutable struct OHBRu_turb_multiple{T<:AbstractRange}
    edges::T
    tau_i::UnitRange{Int}
    N::Array{Int64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
end
function OHBRu_turb_multiple(x_range, tau_i)
    Nx = length(x_range) - 1
    τ_len = length(tau_i)
    mem = zeros(Float64, τ_len)
    mem .= NaN
    OHBRu_turb_multiple(
        x_range,
        tau_i,
        zeros(Int, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        mem
    )
end

function add_data!(ohbr::OHBRu_turb_multiple, X_new)
    X_left = ohbr.mem[end]
    incr2 = X_new - X_left
    if in_range(ohbr.edges, incr2)
        j_bin = find_bin(ohbr.edges, incr2) 
        for (i_tau, X_right) in enumerate(view(ohbr.mem,1:length(ohbr.mem)-1))
            incr1 = X_right - X_left
            ΔX = incr1 - incr2
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
        end
    end
    update_mem!(ohbr.mem, X_new)
    return nothing
end

## Online Kernel algrithms

## Multi step algorithm, 1D

mutable struct OKBRu_turb_multiple{T<:AbstractRange}
    x_eval_points::T
    tau_i::UnitRange{Int}
    w::Array{Float64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
    w_mem::Array{Float64,1}
    kernel::Kernel
    hinv::Float64
end
function OKBRu_turb_multiple(x_eval_points, tau_i, kernel, h::Float64)
    Nx = length(x_eval_points)
    τ_len = length(tau_i)
    mem = zeros(Float64, τ_len)
    mem .= NaN
    hinv = inv(h)
    OKBRu_turb_multiple(
        x_eval_points,
        tau_i,
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        mem,
        zeros(Float64, τ_len),
        kernel,
        hinv
    )
end

function add_data!(okbr::OKBRu_turb_multiple, X_new)
    X_left = okbr.mem[end]
    if isnan(X_left)
        update_mem!(okbr.mem, X_new)
        return nothing
    end
    
    incr2 = X_new - X_left
    for (j_ind, x_eval) in enumerate(okbr.x_eval_points)
        K_weight = apply_kernel(x_eval - incr2, okbr.kernel, okbr.hinv)
        if K_weight > 0.0
            for (i_tau, X_right) in enumerate(view(okbr.mem,1:length(okbr.mem)-1))
                incr1 = X_right - X_left
                ΔX = incr1 - incr2
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
    update_mem!(okbr.mem, X_new)
    return nothing
end