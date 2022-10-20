# Utility functions

# Statistics
mean(x) = sum(x)/length(x)

# Binning (for linear ranges)
in_range(r::LinRange,x::Real) = (x >= r.start) && (x <= r.stop)
find_bin(r::LinRange,x::Real) = floor(Int, ((x-r.start)/(r.stop-r.start)*r.lendiv + 1))
get_bin(r::LinRange,x::Real) = in_range(r,x) ? find_bin(r,x) : 0

# Stack-like memory
function update_mem!(mem::Vector, data)
    for j in Iterators.reverse(Base.OneTo(length(mem)-1))
        mem[j+1] = mem[j]
    end
    mem[1] = data
end
