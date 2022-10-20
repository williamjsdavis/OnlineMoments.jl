# Unused functions and etc...

#TODO: Adding data to multiple OHBR, with N vector instead of matrix
function add_data_N(HBR::MyHBRmultiple3, x_data)
    println(' ')
    println(x_data)
    if in_range(HBR.edges, x_data)
        k = find_bin(HBR.edges, x_data)
        for (i_tau, i_bin) in enumerate(HBR.bin_mem)
            setindex!(HBR.N, HBR.N[i_tau,k] + 1, i_tau, k)
        end
    end
    for (i_tau, i_bin) in enumerate(HBR.bin_mem) if i_bin != 0
        Δx = x_data - HBR.mem[i_tau]
        println((i_tau,HBR.mem[i_tau]))
        mem_tmp = HBR.mem[i_tau]
        HBR.mem[i_tau] = HBR.M1[i_tau,i_bin] # Old mean
        println((i_tau,i_bin,HBR.N))
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
