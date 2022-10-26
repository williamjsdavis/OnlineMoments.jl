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
    N::Array{Float64,1}
    M1::Array{Float64,1}
    M2::Array{Float64,1}
    mem::Float64
    kernel
end
function OKBR_single(x_range::LinRange, kernel)
    Nx = length(x_range)
    OKBR_single(
        x_range,
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        NaN,
        kernel
    )
end
function add_data_old(KBR::OKBR_single, X_new)
    if !isnan(KBR.mem)
        ΔX = X_new - KBR.mem

        for (j_ind,j_xeval) in enumerate(KBR.x_eval_points)
            K_weight = KBR.kernel(j_xeval - KBR.mem)
            KBR.mem = KBR.M1[j_ind] # Old mean?
            setindex!(KBR.N, KBR.N[j_ind] + 1, j_ind)
            setindex!(
                KBR.M1,
                update_mean!(KBR.M1[j_ind], ΔX, KBR.N[j_ind]),
                j_ind
            )
            setindex!(
                KBR.M2,
                update_var!(KBR.M2[j_ind], KBR.M1[j_ind], KBR.mem, ΔX, KBR.N[j_ind]),
                j_ind
            )
            #N[i_ind,j_ind] += K_weight
            #M1[i_ind,j_ind] += K_weight * ΔX
            #M2[i_ind,j_ind] += K_weight * ΔX*ΔX
        end
    end
    KBR.mem = X_new
end
function add_data(KBR, X_new)
    if !isnan(KBR.mem)
        ΔX = X_new - KBR.mem
        for (j_ind,j_xeval) in enumerate(KBR.x_eval_points)
            K_weight = KBR.kernel(j_xeval - KBR.mem)
            if K_weight > 0.0
                KBR.mem = KBR.N[j_ind] # Old weight
                setindex!(KBR.N, KBR.N[j_ind] + K_weight, j_ind)
                tmp = KBR.M1[j_ind]
                setindex!(
                    KBR.M1,
                    update_wmean(KBR.M1[j_ind], KBR.mem, ΔX, K_weight),
                    j_ind
                )
                setindex!(
                    KBR.M2,
                    update_wvar(KBR.M2[j_ind], tmp, KBR.mem, ΔX, KBR.M1[j_ind], K_weight),
                    j_ind
                )
            end
        end
    end
    KBR.mem = X_new
end
