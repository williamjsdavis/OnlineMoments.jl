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

X_sample = deepcopy(X_small)
X_stream() = popfirst!(X_sample)

hinv = inv(h)
kernel_scaled(x) = hinv*kernel(hinv*x)

## KBR

N = 50

# Boxcar kernel (can compare)
M1_K_ref_A, M2_K_ref_A = KBR_moments_A(X_small[1:N], tau_i_range, x_centers, h, kernel)

display(M1_K_ref_A[1,:])

## OKBR
update_wmean(x_bar, w, x_new, w_new) = x_bar + (x_new - x_bar)*(w_new/(w + w_new))
update_wvar(s2, x_bar, w, x_new, x_bar_new, w_new) = (s2*w + w_new*(x_new - x_bar)*(x_new - x_bar_new))/(w + w_new)
function myadd_data!(KBR, X_new)
    if !isnan(KBR.mem)
        ΔX = X_new - KBR.mem
        for (j_ind,j_xeval) in enumerate(KBR.x_eval_points)
            K_weight = KBR.kernel(j_xeval - KBR.mem)
            if K_weight > 0.0
                KBR.mem = KBR.M1[j_ind] # Old mean?
                setindex!(KBR.N, KBR.N[j_ind] + K_weight, j_ind)
                setindex!(
                    KBR.M1,
                    OnlineMoments.update_mean(KBR.M1[j_ind], K_weight * ΔX, KBR.N[j_ind]),
                    j_ind
                )
                setindex!(
                    KBR.M2,
                    OnlineMoments.update_var(KBR.M2[j_ind], KBR.M1[j_ind], KBR.mem, K_weight * ΔX*ΔX, KBR.N[j_ind]),
                    j_ind
                )
                #N[i_ind,j_ind] += K_weight
                #M1[i_ind,j_ind] += K_weight * ΔX
                #M2[i_ind,j_ind] += K_weight * ΔX*ΔX
            end
        end
    end
    KBR.mem = X_new
end
function myadd_dataw!(KBR, X_new)
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
                #update_wvar(s2, x_bar, w, x_new, x_bar_new, w_new)
            end
        end
    end
    KBR.mem = X_new
end

## Test single

hbr_single = OHBR_single(x_edges)
kbr_single = OKBR_single(x_centers, kernel_scaled)

X_data = X_stream()

add_data!(hbr_single, X_data)
myadd_data!(kbr_single, X_data)

## Test N

N = 7

X_sample = deepcopy(X_small)
X_stream() = popfirst!(X_sample)

hbr_single = OHBR_single(x_edges)
kbr_single = OKBR_single(x_centers, kernel_scaled)

M1_K_ref_A, M2_K_ref_A = KBR_moments_A(X_small[1:N], tau_i_range, x_centers, h, kernel)
for _ in 1:N
    local X_data = X_stream()
    add_data!(hbr_single, X_data)
    myadd_dataw!(kbr_single, X_data)
end

display(M1_K_ref_A[1,:])
display(hbr_single.M1)
display(kbr_single.M1)

display(M2_K_ref_A[1,:])
display(hbr_single.M2)
display(kbr_single.M2)

## Testing again

N = 100

X_sample = deepcopy(X_small)
X_stream() = popfirst!(X_sample)

hbr_single = OHBR_single(x_edges)
kbr_single = OKBR_single(x_centers, kernel_scaled)

M1_K_ref_A, M2_K_ref_A = KBR_moments_A(X_small[1:N], tau_i_range, x_centers, h, kernel)
for _ in 1:N
    local X_data = X_stream()
    add_data!(hbr_single, X_data)
    myadd_dataw!(kbr_single, X_data)
end

display(M1_K_ref_A[1,:] .- hbr_single.M1)
display(M1_K_ref_A[1,:] .- kbr_single.M1)

display(M2_K_ref_A[1,:] .- hbr_single.M2)
display(M2_K_ref_A[1,:] .- kbr_single.M2)

## Testing multi OKBR

function myadd_datam!(OKBR, x_data)
    #println(" ")
    #println(x_data)
    for (i_tau, x_left) in enumerate(OKBR.mem)
        if !isnan(x_left)
            ΔX = x_data - x_left
            #println((i_tau,HBR.mem[i_tau]))
            for (j_ind,j_xeval) in enumerate(OKBR.x_eval_points)
                K_weight = OKBR.kernel(j_xeval - x_left)
                if K_weight > 0.0
                    mem_tmp = OKBR.mem[i_tau]
                    OKBR.mem[i_tau] = OKBR.w[i_tau,j_ind] # Old weight
                    setindex!(OKBR.w, OKBR.w[i_tau,j_ind] + K_weight, i_tau, j_ind)
                    tmp = OKBR.M1[i_tau, j_ind]
                    setindex!(
                        OKBR.M1,
                        OnlineMoments.update_wmean(OKBR.M1[i_tau, j_ind], OKBR.mem[i_tau], ΔX, K_weight),
                        i_tau, j_ind
                    )
                    setindex!(
                        OKBR.M2,
                        OnlineMoments.update_wvar(OKBR.M2[i_tau, j_ind], tmp, OKBR.mem[i_tau], ΔX, OKBR.M1[i_tau, j_ind], K_weight),
                        i_tau, j_ind
                    )
                    OKBR.mem[i_tau] = mem_tmp
                end
            end
        end
    end
    OnlineMoments.update_mem!(OKBR.mem, x_data)
end

kbr_multiple = OKBR_multiple(x_centers, N_tau, kernel_scaled)

X_sample = deepcopy(X_small)
X_stream() = popfirst!(X_sample)

X_data_sample = X_stream()
myadd_datam!(kbr_multiple, X_data_sample)
X_data_sample = X_stream()
myadd_datam!(kbr_multiple, X_data_sample)

## Testing again

hbr_multiple = OHBR_multiple(x_edges, N_tau)
kbr_multiple = OKBR_multiple(x_centers, N_tau, kernel_scaled)

X_sample = deepcopy(X_small)
X_stream() = popfirst!(X_sample)

N = 100
M1_K_ref_A, M2_K_ref_A = KBR_moments_A(X_small[1:N], tau_i_range, x_centers, h, kernel)
for _ in 1:N
    local X_data = X_stream()
    add_data!(hbr_multiple, X_data)
    myadd_datam!(kbr_multiple, X_data)
end

println("M1 full")
display(M1_K_ref_A)
display(hbr_multiple.M1)
display(kbr_multiple.M1)
display(M1_K_ref_A .- kbr_multiple.M1)

println("M2 full")
display(M2_K_ref_A)
display(hbr_multiple.M2)
display(kbr_multiple.M2)
display(M2_K_ref_A .- kbr_multiple.M2)

## Testing separately

X_sample = deepcopy(X_small)
X_stream() = popfirst!(X_sample)

X_data_sample = X_stream()
add_data!(hbr_multiple, X_data_sample)
myadd_datam!(kbr_multiple, X_data_sample)
