_meanoffset(y, σ2) = √(-2σ2 * log(√(2π * σ2) * y))
ipdf_right(dist::Normal, y) = mean(dist) + _meanoffset(y, var(dist))
ipdf_left(dist::Normal, y) = mean(dist) - _meanoffset(y, var(dist))

function ipdf_right(dist::Truncated{<:Normal}, y)
    ut = dist.untruncated
    mode_pd = pdf(dist, mode(dist))
    if y > mode_pd
        throw(DomainError(
            y,
            "y must be less than or equal to the density at the mode = $mode_pd"
        ))
    end

    y < pdf(dist, maximum(dist)) ? maximum(dist) : ipdf_right(ut, y * dist.tp)
end

function ipdf_left(dist::Truncated{<:Normal}, y)
    ut = dist.untruncated
    mode_pd = pdf(dist, mode(dist))
    if y > mode_pd
        throw(DomainError(
            y,
            "y must be less than or equal to the density at the mode = $mode_pd"
        ))
    end

    y < pdf(dist, minimum(dist)) ? minimum(dist) : ipdf_left(ut, y * dist.tp)
end

ipdf_right(dist::Exponential, y) = -scale(dist) * log(scale(dist) * y)
