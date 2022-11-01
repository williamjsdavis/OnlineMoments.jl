# Utility functions

# Binning (for linear ranges)
in_range(r::LinRange, x::Real) = (x >= r.start) && (x <= r.stop)
find_bin(r::LinRange, x::Real) = floor(Int, ((x-r.start)/(r.stop-r.start)*r.lendiv + 1))

in_range(r::StepRangeLen, x::Real) = (x >= first(r)) && (x <= last(r))
find_bin(r::StepRangeLen, x::Real) = floor(Int, ((x-first(r))/(last(r)-first(r))*(length(r) - 1) + 1))

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
