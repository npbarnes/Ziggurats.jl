module Normal

using SpecialFunctions
export f, finv, F, tailsample

f(x) = 1/√(2π) * exp(-x^2/2)
finv(y) = √(-2log(√(2π) * y))
F(x) = 1/2 * (1 + erf(x/√2))
function tailsample(x1)
    while true
        x = -log(rand()) / x1
        y = -log(rand())
        if 2y > x^2
            return x + x1
        end
    end
end

end # module