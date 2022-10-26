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
function myadd_data(KBR, X_new)
    if !isnan(KBR.mem)
        ΔX = X_new - KBR.mem
        for (j_ind,j_xeval) in enumerate(KBR.x_eval_points)
            K_weight = KBR.kernel(j_xeval - KBR.mem)
            if K_weight > 0.0
                KBR.mem = KBR.M1[j_ind] # Old mean?
                setindex!(KBR.N, KBR.N[j_ind] + K_weight, j_ind)
                setindex!(
                    KBR.M1,
                    OnlineMoments.update_mean!(KBR.M1[j_ind], K_weight * ΔX, KBR.N[j_ind]),
                    j_ind
                )
                setindex!(
                    KBR.M2,
                    OnlineMoments.update_var!(KBR.M2[j_ind], KBR.M1[j_ind], KBR.mem, K_weight * ΔX*ΔX, KBR.N[j_ind]),
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
function myadd_dataw(KBR, X_new)
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

add_data(hbr_single, X_data)
myadd_data(kbr_single, X_data)

### Test N

N = 7

X_sample = deepcopy(X_small)
X_stream() = popfirst!(X_sample)

hbr_single = OHBR_single(x_edges)
kbr_single = OKBR_single(x_centers, kernel_scaled)

M1_K_ref_A, M2_K_ref_A = KBR_moments_A(X_small[1:N], tau_i_range, x_centers, h, kernel)
for _ in 1:N
    X_data = X_stream()
    add_data(hbr_single, X_data)
    myadd_dataw(kbr_single, X_data)
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
    X_data = X_stream()
    add_data(hbr_single, X_data)
    myadd_dataw(kbr_single, X_data)
end

display(M1_K_ref_A[1,:] .- hbr_single.M1)
display(M1_K_ref_A[1,:] .- kbr_single.M1)

display(M2_K_ref_A[1,:] .- hbr_single.M2)
display(M2_K_ref_A[1,:] .- kbr_single.M2)