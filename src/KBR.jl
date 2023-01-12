# Kernel Based Regression (KBR)

function KBR_moments_A(X, tau_i_range, x_range, h, kernel::Kernel)
    n = length(X)
    nτ = length(tau_i_range)
    nx = length(x_range)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)

    hinv = inv(h)
    kernel_scaled(x) = hinv*kernel(hinv*x)

    for (i_left, X_left) in enumerate(X[1:end-1])
        for (i_ind,i_tau) in enumerate(tau_i_range)
            i_right = i_left + i_tau
            if i_right <= n
                ΔX = X[i_right] - X_left
                for (j_ind,x_eval) in enumerate(x_range)
                    K_weight = kernel_scaled(x_eval - X_left)
                    N[i_ind,j_ind] += K_weight
                    M1[i_ind,j_ind] += K_weight * ΔX
                    M2[i_ind,j_ind] += K_weight * ΔX*ΔX
                end
            end
        end
    end

    #NOTE: No Bessel correction
    for i in 1:length(N)
        M2[i] = (M2[i] - (M1[i]*M1[i])/N[i]) / N[i]
        M1[i] = M1[i] / N[i]
    end

    return M1, M2
end
function KBR_moments_A2(X, tau_i_range, x_range, h, kernel)
    n = length(X)
    nτ = length(tau_i_range)
    nx = length(x_range)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)

    hinv = inv(h)
    kernel_scaled(x) = hinv*kernel(hinv*x)

    for (i_left, X_left) in enumerate(X[1:end-1])
        for (j_ind, x_eval) in enumerate(x_range)
            K_weight = kernel_scaled(x_eval - X_left)
            for (i_ind,i_tau) in enumerate(tau_i_range)
                i_right = i_left + i_tau
                if i_right <= n
                    ΔX = X[i_right] - X_left
                    N[i_ind,j_ind] += K_weight
                    M1[i_ind,j_ind] += K_weight * ΔX
                    M2[i_ind,j_ind] += K_weight * ΔX*ΔX
                end
            end
        end
    end

    #NOTE: No Bessel correction
    for i in 1:length(N)
        M2[i] = (M2[i] - (M1[i]*M1[i])/N[i]) / N[i]
        M1[i] = M1[i] / N[i]
    end

    return M1, M2
end

## Modulo moments
d_mod(x,n) = min(mod(x,n),mod(-x,n))
function KBR_moments_mod(X, tau_i_range, x_range, h, kernel, period)
    n = length(X)
    nτ = length(tau_i_range)
    nx = length(x_range)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)

    hinv = inv(h)
    kernel_scaled(x) = hinv*kernel(hinv*d_mod(x,period))

    for (i_left, X_left) in enumerate(X[1:end-1])
        for (i_ind,i_tau) in enumerate(tau_i_range)
            i_right = i_left + i_tau
            if i_right <= n
                ΔX = X[i_right] - X_left
                for (j_ind,x_eval) in enumerate(x_range)
                    K_weight = kernel_scaled(x_eval - X_left)
                    N[i_ind,j_ind] += K_weight
                    M1[i_ind,j_ind] += K_weight * ΔX
                    M2[i_ind,j_ind] += K_weight * ΔX*ΔX
                end
            end
        end
    end

    #NOTE: No Bessel correction
    for i in 1:length(N)
        M2[i] = (M2[i] - (M1[i]*M1[i])/N[i]) / N[i]
        M1[i] = M1[i] / N[i]
    end

    return M1, M2
end
