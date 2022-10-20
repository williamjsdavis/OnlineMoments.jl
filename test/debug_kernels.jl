using OnlineMoments
using JLD2, FileIO
include("./load_test_data.jl")

# Test settings
const tau_i_range = 1:1:6
const x_edges = LinRange(-0.0,0.4,5)

# Load test data
const X_small = load_small_data()
#X_small = X_small[1:10]

# Derived test variables
x_centers = 0.5*(x_edges[1:end-1] + x_edges[2:end])
N_tau = length(tau_i_range)
N_x = length(x_centers)

h = 0.5*step(x_centers)
kernel = Boxcar()

##

function myHBR_moments_C(X, tau_vector, edge_vector)
    nτ = length(tau_vector)
    nx = length(edge_vector) - 1
    nX = length(X)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)
    data = [[] for i in 1:nτ, j in 1:nx]
    mem = 0.0

    for (i_left, X_left) in enumerate(X[1:end-1])
        println("Data: X = $X_left")
        if OnlineMoments.in_range(edge_vector, X_left)
            println("In range")
            k = OnlineMoments.find_bin(edge_vector, X_left)
            println("Bin: k = $k")
            for (j,tau) in enumerate(tau_vector)
                ii = i_left + j
                if ii <= nX
                    inc = X[ii] - X_left
                    push!(data[j, k], inc)
                    println("Diff: ΔX = $inc")
                    mem = M1[j, k]
                    setindex!(N, N[j,k] + 1, j, k)
                    setindex!(
                        M1,
                        OnlineMoments.update_mean!(M1[j, k], inc, N[j, k]),
                        j, k
                    )
                    setindex!(
                        M2,
                        OnlineMoments.update_var!(M2[j, k], M1[j, k], mem, inc, N[j, k]),
                        j, k
                    )
                end
            end
        end
    end

    return M1, M2, data
end
function myKBR_moments_A(X,tau_i_range,x_range,h,kernel)
    n = length(X)
    nτ = length(tau_i_range)
    nx = length(x_range)
    N = zeros(nτ,nx)
    M1 = zeros(nτ,nx)
    M2 = zeros(nτ,nx)
    data = [[] for i in 1:nτ, j in 1:nx]

    hinv = inv(h)
    kernel_scaled(x) = hinv*kernel(hinv*x)

    for (i_left, X_left) in enumerate(X[1:end-1])
        println("Data: X = $X_left")
        for (i_ind,i_tau) in enumerate(tau_i_range)
            i_right = i_left + i_tau
            if i_right <= n
                ΔX = X[i_right] - X_left
                println("Diff: ΔX = $ΔX")
                for (j_ind,j_xeval) in enumerate(x_range)
                    println("x eval: x = $j_xeval")
                    K_weight = kernel_scaled(j_xeval - X_left)
                    if K_weight != 0
                        push!(data[i_ind,j_ind], ΔX)
                    end
                    println("K_weight: K = $K_weight")
                    N[i_ind,j_ind] += K_weight
                    M1[i_ind,j_ind] += K_weight * ΔX
                    M2[i_ind,j_ind] += K_weight * ΔX*ΔX
                end
            end
        end
    end

    for i in 1:length(N)
        # This line seems to be the problem
        #M2[i] = (M2[i] - (M1[i]*M1[i])/N[i]) / (N[i] - 1)
        M2[i] = (M2[i] - (M1[i]*M1[i])/N[i]) / N[i]
        M1[i] = M1[i] / N[i]
    end

    return M1, M2, data
end

i = 100
X_new = X_small[1:i]

M1_ref_C, M2_ref_C, data_ref_C = myHBR_moments_C(X_new, tau_i_range, x_edges)
println("New")
M1_K_ref_A, M2_K_ref_A, data_K_ref_A = myKBR_moments_A(X_new, tau_i_range, x_centers, h, kernel)

M1_diff = M1_K_ref_A .- M1_ref_C
M2_diff = M2_K_ref_A .- M2_ref_C

display(M1_diff)
display(M2_diff)

display(M2_ref_C)
display(M2_K_ref_A)

display(data_ref_C)
display(data_K_ref_A)

## Variance algorithms

function var_naive(X)
    n, s, ss = 0, 0, 0
    for x in X
        n += 1
        s += x
        ss += x*x
    end
    return (ss - (s*s)/n) / n
end
