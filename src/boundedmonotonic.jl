function buildziggurat_bounded!(x, y, L, R, pdf, ipdf, A)
    Lpd = pdf(L)
    Rpd = pdf(R)
    isincreasing = Lpd < Rpd

    mode = isincreasing ? R : L
    mode_pd = isincreasing ? Rpd : Lpd

    x[1] = isincreasing ? L : R
    y[1] = A/abs(x[1] - mode)
    
    for i in eachindex(x)[begin:end-1]
        if y[i] >= mode_pd
            return nothing
        end
        x[i+1] = ipdf(y[i])
        y[i+1] = y[i] + A/abs(x[i] - mode)
    end

    x, y
end
