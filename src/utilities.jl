function between(a, b, x)
    l, r = minmax(a, b)
    l <= x <= r
end

function regularize_domain(domain)
    unique(sort!(collect(promote(float.(domain)...))))
end
