# Histogram Based Regression (HBR), uncorrected

# Algorithm A: loop over moment elements (no variance correction)
function HBRu_moments_A(X, tau_vector, edge_vector)
    nτ = length(tau_vector)
    nx = length(edge_vector) - 1
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)
    bin_index = map(X) do y
        in_range(edge_vector, y) ? find_bin(edge_vector, y) : 0
    end
    for (i,tau) in enumerate(tau_vector)
        for j in 1:(length(edge_vector)-1)
            iX = findall(bin_index[1:end-tau] .== j)
            ΔX = X[iX .+ tau] .- X[iX]
            M1[i,j] = mean(ΔX)
            M2[i,j] = mean(ΔX.^2)
        end
    end
    return M1, M2
end

# Algorithm B: map over X elements (no variance correction)
function HBRu_moments_B(X::AbstractVector, tau_vector, edge_vector::AbstractRange)
    ti_grid = make_grid(tau_vector, edge_vector)

    bin_index = map(X) do y
        in_range(edge_vector, y) ? find_bin(edge_vector, y) : 0
    end

    bin_groups = map(1:(length(edge_vector)-1)) do i
        findall(bin_index[1:end-maximum(tau_vector)].==i)
    end

    M1M2 = map(ti_grid) do ti
        get_moments_uncorrected(X, bin_groups[ti[2]], ti[1])
    end

    M1 = broadcast(first,M1M2)
    M2 = broadcast(last,M1M2)

    return M1, M2
end
function make_grid(tau_vector::AbstractVector, edge_vector::LinRange)
    return broadcast((x,y)->(x,y), tau_vector, (1:edge_vector.lendiv)')
end
function make_grid(tau_vector::AbstractVector, edge_vector::AbstractRange)
    return broadcast((x,y)->(x,y), tau_vector, (1:(length(edge_vector)-1))')
end
function get_moments_uncorrected(X, iX, tau)
    ΔX = X[iX .+ tau] .- X[iX]
    M1 = mean(ΔX)
    M2 = mean(ΔX.^2)
    return M1, M2
end

#TODO: Fix badly labelled variables here
# Algorithm C: Loop over X elements (no variance correction)
function HBRu_moments_C(X, tau_vector, edge_vector)
    nτ = length(tau_vector)
    nx = length(edge_vector) - 1
    nX = length(X)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)
    mem = 0.0

    for (i_left, X_left) in enumerate(X[1:end-1])
        if in_range(edge_vector, X_left)
            k = find_bin(edge_vector, X_left)
            for (j,tau) in enumerate(tau_vector)
                ii = i_left + j
                if ii <= nX
                    ΔX = X[ii] - X_left
                    mem = M1[j, k]
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

# Algorithm C: Loop over X elements (different element order) (no variance correction)
function HBRu_moments_C2(X, tau_vector, edge_vector)
    nτ = length(tau_vector)
    nx = length(edge_vector) - 1
    nX = length(X)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)
    mem = 0.0

    # Bins
    b = broadcast(x->get_bin(edge_vector, x), X)

    for (i_ind, i_tau) in enumerate(tau_vector)
        for (i_left, X_left) in enumerate(X[1:end-i_tau])
            j_bin = b[i_left]
            if j_bin != 0
                ΔX = X[i_left+i_tau] - X_left
                setindex!(N, N[i_ind, j_bin] + 1, i_ind, j_bin)
                setindex!(
                    M1,
                    update_mean(M1[i_ind, j_bin], ΔX, N[i_ind, j_bin]),
                    i_ind, j_bin
                )
                setindex!(
                    M2,
                    update_ss(M2[i_ind, j_bin], ΔX, N[i_ind, j_bin]),
                    i_ind, j_bin
                )
            end
        end
    end

    return M1, M2
end

#TODO: Fix badly labelled variables here
## Modulo moments (no variance correction)
function HBRu_moments_mod(X, tau_vector, edge_vector, period)
    nτ = length(tau_vector)
    nx = length(edge_vector) - 1
    nX = length(X)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)

    for (i_left, X_left) in enumerate(X[1:end-1])
        k = find_mod_bin(edge_vector, period, X_left)
        if k != 0
            for (j,tau) in enumerate(tau_vector)
                ii = i_left + j
                if ii <= nX
                    ΔX = X[ii] - X_left
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
