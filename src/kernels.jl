# Kernels
#TODO: Add kernel scale into struct

abstract type Kernel end
#TODO: Add "Support" pair to Kernel? For easily checking if a value is in the support?

# Epaneknikov
struct Epaneknikov <: Kernel
end

#NOTE: Maybe change this to support [-1,1]?
const epan_scale = 3*sqrt(5)/100

function (k::Epaneknikov)(x)
    x2 = x*x
    if x2 < 5.
        return epan_scale*(5. - x2)
    else
        return 0.
    end
end

# Boxcar
struct Boxcar <: Kernel
end

(k::Boxcar)(x) = (abs(x) < 1.) ? 0.5 : 0.
