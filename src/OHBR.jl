# Online Histogram Based Regression (OHBR)

## Single step algorithm, 1D (for testing)

mutable struct OHBR_single
    edges::LinRange{Float64}
    N::Array{Int64,1}
    M1::Array{Float64,1}
    M2::Array{Float64,1}
    mem::Float64
end
function OHBR_single(x_range::LinRange)
    Nx = length(x_range) - 1
    OHBR_single(
        x_range,
        zeros(Int, Nx),
        zeros(Float64, Nx),
        zeros(Float64, Nx),
        NaN
    )
end
function add_data!(HBR::OHBR_single, X_new)
    if in_range(HBR.edges, HBR.mem)
        ΔX = X_new - HBR.mem
        i = find_bin(HBR.edges, HBR.mem)
        HBR.mem = HBR.M1[i]
        setindex!(HBR.N, HBR.N[i] + 1, i)
        setindex!(
            HBR.M1,
            update_mean(HBR.M1[i], ΔX, HBR.N[i]),
            i
        )
        setindex!(
            HBR.M2,
            update_var(HBR.M2[i], HBR.M1[i], HBR.mem, ΔX, HBR.N[i]),
            i
        )
    end
    HBR.mem = X_new
    return nothing
end
#TODO: Make this return nothing and be consistant

## Multi step algorithm, 1D

mutable struct OHBR_multiple
    edges::LinRange{Float64}
    N::Array{Int64,2}
    M1::Array{Float64,2}
    M2::Array{Float64,2}
    mem::Array{Float64,1}
    bin_mem::Array{Int64,1}
end
function OHBR_multiple(x_range::LinRange, τ_len::Integer)
    #TODO: make τ_len consistent with length or array
    Nx = length(x_range) - 1
    mem = zeros(Float64, τ_len)
    mem .= NaN
    OHBR_multiple(
        x_range,
        zeros(Int, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        zeros(Float64, τ_len, Nx),
        mem,
        zeros(Int, τ_len)
    )
end
function add_data!(OHBR::OHBR_multiple, x_data)
    for (i_tau, i_bin) in enumerate(OHBR.bin_mem) if i_bin != 0
        Δx = x_data - OHBR.mem[i_tau]
        mem_tmp = OHBR.mem[i_tau]
        OHBR.mem[i_tau] = OHBR.M1[i_tau,i_bin] # Old mean
        setindex!(OHBR.N, OHBR.N[i_tau,i_bin] + 1, i_tau, i_bin)
        setindex!(
            OHBR.M1,
            update_mean(OHBR.M1[i_tau,i_bin], Δx, OHBR.N[i_tau,i_bin]),
            i_tau, i_bin
        )
        setindex!(
            OHBR.M2,
            update_var(OHBR.M2[i_tau,i_bin], OHBR.M1[i_tau,i_bin], OHBR.mem[i_tau], Δx, OHBR.N[i_tau,i_bin]),
            i_tau, i_bin
        )
        OHBR.mem[i_tau] = mem_tmp
    end end
    update_mem!(OHBR.mem, x_data)
    update_mem!(OHBR.bin_mem, get_bin(OHBR.edges, x_data))
    return nothing
end
