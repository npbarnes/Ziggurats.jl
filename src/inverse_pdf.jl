_meanoffset(y, σ2) = √(-2σ2 * log(√(2π*σ2)*y))
ipdf_right(dist::Normal, y) = mean(dist) + _meanoffset(y, var(dist))
ipdf_left(dist::Normal, y) =  mean(dist) - _meanoffset(y, var(dist))

function ipdf_right(dist::Truncated{<:Normal}, y)
    ut = dist.untruncated
    @assert dist.upper > mean(ut)

    L = max(dist.lower, mean(ut))
    R = dist.upper

    @assert y <= pdf(dist, L)

    if y <= pdf(dist, R)
        return R
    else
        return ipdf_right(ut, y*dist.tp)
    end
end

ipdf_right(dist::Exponential, y) = -scale(dist)*log(scale(dist)*y)
