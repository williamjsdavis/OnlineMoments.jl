# Histogram Based Regression (HBR)

# Algorithm A: loop over moment elements
function HBR_moments_A(X,tau_vector,edge_vector)
    nτ = length(tau_vector)
    nx = length(edge_vector) - 1
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)
    bin_index = map(X) do y
        in_range(edge_vector,y) ? find_bin(edge_vector,y) : 0
    end
    for (i,tau) in enumerate(tau_vector)
        for j in 1:edge_vector.lendiv
            iX = findall(bin_index[1:end-tau] .== j)
            inc = X[iX .+ tau] .- X[iX]
            M1[i,j] = mean(inc)
            residuals = inc .- M1[i,j]
            M2[i,j] = mean(residuals.^2)
        end
    end
    return M1, M2
end

# Algorithm A: loop over moment elements (no variance correction)
function HBR_moments_A2(X,tau_vector,edge_vector)
    nτ = length(tau_vector)
    nx = length(edge_vector) - 1
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)
    bin_index = map(X) do y
        in_range(edge_vector,y) ? find_bin(edge_vector,y) : 0
    end
    for (i,tau) in enumerate(tau_vector)
        for j in 1:edge_vector.lendiv
            iX = findall(bin_index[1:end-tau] .== j)
            inc = X[iX .+ tau] .- X[iX]
            M1[i,j] = mean(inc)
            M2[i,j] = mean(inc.^2)
        end
    end
    return M1, M2
end

# Algorithm B: map over X elements
function HBR_moments_B(X::AbstractVector,tau_vector,edge_vector::LinRange)
    ti_grid = make_grid(tau_vector,edge_vector)

    bin_index = map(X) do y
        in_range(edge_vector,y) ? find_bin(edge_vector,y) : 0
    end

    bin_groups = map(1:edge_vector.lendiv) do i
        findall(bin_index[1:end-maximum(tau_vector)].==i)
    end

    M1M2 = map(ti_grid) do ti
        get_moments(X,bin_groups[ti[2]],ti[1])
    end

    M1 = broadcast(first,M1M2)
    M2 = broadcast(last,M1M2)

    return M1, M2
end
function make_grid(tau_vector::AbstractVector,edge_vector::LinRange)
    return broadcast((x,y)->(x,y),tau_vector,(1:edge_vector.lendiv)')
end
function get_moments(X, iX, tau)
    inc = X[iX .+ tau] .- X[iX]
    M1 = mean(inc)
    residuals = inc .- M1
    M2 = mean(residuals.^2)
    return M1, M2
end

# Algorithm C: Loop over X elements
update_mean!(x_bar, x_new, n) = x_bar + (x_new - x_bar)/n
update_var!(s2, x_bar, x_bar_old, x_new, n) = s2 + ((x_new - x_bar)*(x_new - x_bar_old) - s2)/n
function HBR_moments_C(X, tau_vector, edge_vector)
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
                    inc = X[ii] - X_left
                    mem = M1[j, k]
                    setindex!(N, N[j,k] + 1, j, k)
                    setindex!(
                        M1,
                        update_mean!(M1[j, k], inc, N[j, k]),
                        j, k
                    )
                    setindex!(
                        M2,
                        update_var!(M2[j, k], M1[j, k], mem, inc, N[j, k]),
                        j, k
                    )
                end
            end
        end
    end

    return M1, M2
end

# Algorithm C: Loop over X elements (different element order)
function HBR_moments_C2(X, tau_vector, edge_vector)
    nτ = length(tau_vector)
    nx = length(edge_vector) - 1
    nX = length(X)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)
    mem = 0.0

    # Bins
    b = broadcast(x->get_bin(edge_vector, x), X)

    for (j,tau) in enumerate(tau_vector)
        for (i_left, X_left) in enumerate(X[1:end-j])
            bin = b[i_left]
            if bin != 0
                inc = X[i_left+j] - X_left
                mem = M1[j, bin]
                setindex!(N, N[j, bin] + 1, j, bin)
                setindex!(
                    M1,
                    update_mean!(M1[j, bin], inc, N[j, bin]),
                    j, bin
                )
                setindex!(
                    M2,
                    update_var!(M2[j, bin], M1[j, bin], mem, inc, N[j, bin]),
                    j, bin
                )
            end
        end
    end

    return M1, M2
end