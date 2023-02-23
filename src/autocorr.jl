# Online autocorrelation

## Offline method

"Offline autocorrelation"
function offline_autocorr(X,lags::Int)
    N = length(X)
    X_mean = sum(X)/N
    X_demean = X .- X_mean
    
    acf = zeros(lags+1)
    for i in 1:N
        for j in 1:(lags+1)
            k = i + j - 1
            k > N && break
            acf[j] += X_demean[k]*X_demean[i]
        end
    end
    for j in 1:(lags+1)
        acf[j] /= N - j + 1
    end
    return acf
end

## Online method
#TODO: Add this properly without OnlineStats

"Add data point to sample autocorrelation"
add_data!(acf::AutoCov, X_right) = fit!(acf,X_right)

online_autocorr(acf::AutoCov) = autocov(acf)
