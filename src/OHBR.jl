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
function add_data(HBR::OHBR_single, X_new)
    if in_range(HBR.edges, HBR.mem)
        ΔX = X_new - HBR.mem
        i = find_bin(HBR.edges, HBR.mem)
        HBR.mem = HBR.M1[i]
        setindex!(HBR.N, HBR.N[i] + 1, i)
        setindex!(
            HBR.M1,
            update_mean!(HBR.M1[i], ΔX, HBR.N[i]),
            i
        )
        setindex!(
            HBR.M2,
            update_var!(HBR.M2[i], HBR.M1[i], HBR.mem, ΔX, HBR.N[i]),
            i
        )
    end
    HBR.mem = X_new
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
function add_data(HBR::OHBR_multiple, x_data)
    #println(' ')
    #println(x_data)
    for (i_tau, i_bin) in enumerate(HBR.bin_mem) if i_bin != 0
        Δx = x_data - HBR.mem[i_tau]
        #println((i_tau,HBR.mem[i_tau]))
        mem_tmp = HBR.mem[i_tau]
        HBR.mem[i_tau] = HBR.M1[i_tau,i_bin] # Old mean
        setindex!(HBR.N, HBR.N[i_tau,i_bin] + 1, i_tau, i_bin)
        #println((HBR.N, (i_tau,i_bin), HBR.N[i_tau,i_bin]))
        setindex!(
            HBR.M1,
            update_mean!(HBR.M1[i_tau,i_bin], Δx, HBR.N[i_tau,i_bin]),
            i_tau, i_bin
        )
        setindex!(
            HBR.M2,
            update_var!(HBR.M2[i_tau,i_bin], HBR.M1[i_tau,i_bin], HBR.mem[i_tau], Δx, HBR.N[i_tau,i_bin]),
            i_tau, i_bin
        )
        HBR.mem[i_tau] = mem_tmp
    end end
    update_mem!(HBR.mem, x_data)
    update_mem!(HBR.bin_mem, get_bin(HBR.edges, x_data))
end
