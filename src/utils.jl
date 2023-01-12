# Utility functions

# Binning (for linear ranges)
in_range(r::LinRange, x::Real) = (x >= r.start) && (x <= r.stop)
find_bin(r::LinRange, x::Real) = floor(Int, ((x-r.start)/(r.stop-r.start)*r.lendiv + 1))

in_range(r::StepRangeLen, x::Real) = (x >= first(r)) && (x <= last(r))
find_bin(r::StepRangeLen, x::Real) = floor(Int, ((x-first(r))/(last(r)-first(r))*(length(r) - 1) + 1))

# Modulo binning
d_plus(a,b,n) = mod(b-a,n)
"Half-open interval"
is_in_interval_mod(a,b,n,x) = d_plus(a,x,n) < d_plus(a,b,n)
function find_mod_bin(edge_vector, n, X)
    for (i,a) in enumerate(edge_vector[1:end-1])
        b = edge_vector[i+1]
        if is_in_interval_mod(a,b,n,X)
            return i
        end
    end
    return 0
end

# Test
#map(x->(find_bin(range(0.0, 0.4, 5), x), find_bin(LinRange(0.0, 0.4, 5), x)), range(-0.1,0.5,10))

get_bin(r,x) = in_range(r,x) ? find_bin(r,x) : 0

# Stack-like memory
function update_mem!(mem::Vector, data)
    for j in Iterators.reverse(Base.OneTo(length(mem)-1))
        mem[j+1] = mem[j]
    end
    mem[1] = data
end
