module LogBeta

import SpecialFunctions

"""
    logbeta(α::Float64, β::Float64)

Logarithm of the Beta function, log B(α, β).
"""
logbeta(α::Float64, β::Float64) = SpecialFunctions.logbeta(α, β)

"""
    logbeta(α::Float64, β::Float64, x::Float64)

Logarithm of the incomplete Beta function, log B(α, β, x).
"""
function logbeta(α::Float64, β::Float64, x::Float64)
    @assert α > 0 && β > 0 && 0 ≤ x ≤ 1
    if x == 0.0
        return 0.0
    elseif x == 1.0
        return SpecialFunctions.logbeta(α, β)
    elseif x ≤ (α + 1) / (α + β + 2)
        return α * log(x) + β * log(1 - x) - log(α) + lbeta_cf(α, β, x)
    else
        return log_sub(logbeta(α, β), logbeta(β, α, 1 - x))
    end
end





"""
    logbeta(α::Float64, β::Float64, x₁::Float64, x₂::Float64)
Logarithm of the incomplete Beta function,
    log B(α, β, x₁, x₂) = log(B(α, β, x₂) - B(α, β, x₁))
"""

function logbeta(α::Float64, β::Float64, x₁::Float64, x₂::Float64)
    @assert α > 0 && β > 0 && 0 ≤ x₁ ≤ x₂ ≤ 1
    if x₁ == 0
        return logbeta(α, β, x₂)
    elseif x₂ == 1
        return logbeta(β, α, 1 - x₁)
    elseif α ≤ 1 || β ≤ 1
        return log(betar(α, β, x₁, x₂)) + logbeta(α, β)
    end

    mode = (α - 1) / (α + β - 2)

    if x₂ ≤ mode
        return log_sub(logbeta(α, β, x₂), logbeta(α, β, x₁))
    elseif mode ≤ x₁
        return log_sub(logbeta(β, α, 1 - x₁), logbeta(β, α, 1 - x₂))
    else
        return log_sub(logbeta(α, β), log_add(logbeta(α, β, x₁), logbeta(β, α, 1 - x₂)))
    end
end



"""
    betar(α::Float64, β::Float64, x::Float64)
    betar(α::Float64, β::Float64, x₁::Float64, x₂::Float64)

Computes the regularized incomplete Beta function,

    Br(α, β, x) = B(α, β, x) / B(α, β)
    Br(α, β, x₂, x₁) = (B(α, β, x₂) - B(α, β, x₁)) / B(α, β) = Br(α, β, x₂) - Br(α, β, x₁)
"""
function betar end

betar(α::Float64, β::Float64, x::Float64) = exp(logbeta(α, β, x) - logbeta(α, β))

function betar(α::Float64, β::Float64, x₁::Float64, x₂::Float64)
    @assert α > 0 && β > 0 && 0 ≤ x₁ ≤ x₂ ≤ 1
    if x₁ == 0
        return betar(α, β, x₂)
    elseif x₂ == 1
        return betar(β, α, 1 - x₁)
    elseif α ≤ 1 || β ≤ 1
        return betar(α, β, x₂) - betar(α, β, x₁)
    end

    mode = (α - 1) / (α + β - 2)

    if x₂ ≤ mode
        return betar(α, β, x₂) - betar(α, β, x₁)
    elseif x₁ ≥ mode
        return betar(β, α, 1 - x₁) - betar(β, α, 1 - x₂)
    else
        return 1 - betar(α, β, x₁) - betar(β, α, 1 - x₂)
    end
end

"""
    log1exp(x::Float64)

Stable computation of log(1 - exp(-x)), for x ≥ 0.
"""
function log1exp(x::Float64)
    @assert x ≥ 0
    if x > log(2); log(-expm1(-x)) else log1p(-exp(-x)) end
end

"""
    log_sub(logx::Float64, logy::Float64)

Stable computation of log(|x-y|) from log(x) and log(y)
"""
function log_sub(logx::Float64, logy::Float64)
    logmin, logmax = minmax(logx, logy)
    if logmin > -Inf
        return logmax + log1exp(logmax - logmin)
    else
        return logmax
    end
end

"""
    log_add(logx::Float64, logy::Float64)

Stable computation of log(x + y) from log(x) and log(y)
"""
function log_add(logx::Float64, logy::Float64)
    logmin, logmax = minmax(logx, logy)
    if -Inf < logmin && logmax < Inf
        return logmax + log1p(exp(logmin - logmax))
    else
        return logmax
    end
end



"""
    lbeta_cf(α::Float64, β::Float64, x::Float64)

Computes the logarithm of the incomplete Beta function, by a continued fraction expansion.
See Eq. 6.4.5 of Numerical Recipes 3rd Edition.
Assumes that x ≤ (α + 1) / (α + β + 2).
This is an internal function, you are not supposed to call it.
"""
function lbeta_cf(α::Float64, β::Float64, x::Float64)
    @assert α > 0 && β > 0 && 0 ≤ x ≤ 1
    @assert x ≤ (α + 1) / (α + β + 2)

    # Evaluates the continued fraction for the incomplete Beta function by the modified Lentz's method

    MAXIT = 10^3
    ϵ = 1e-10
    fpmin = 1e-30

    # some constant factors
    qab = α + β
    qap = α + 1.0
    qam = α - 1.0

    # first iteration of Lentz's method
    C = 1.0
    D = 1.0 - qab * x / qap
    if abs(D) < fpmin; D = fpmin end
    D = 1.0 / D
    @assert D > 0

    # function estimate
    log_cf = log(D)

    for m = 1:MAXIT
        m2::Int = 2m

        # one step (the even one) of the recurrence. The index here is 2m
        aa = m * (β - m) * x / ((qam + m2) * (α + m2))

        D = 1.0 + aa * D
        if abs(D) < fpmin; D = fpmin end
        D = 1.0 / D

        C = 1.0 + aa / C
        if abs(C) < fpmin; C = fpmin end

        log_cf += log(D) + log(C)

        # next step (the odd one) of the recurrence. The index here is 2m + 1
        aa = -(α + m) * (qab + m) * x / ((α + m2) * (qap + m2))

        D = 1.0 + aa * D
        if abs(D) < fpmin; D = fpmin end
        D = 1.0 / D

        C = 1.0 + aa / C
        if abs(C) < fpmin; C = fpmin end

        log_del = log(D) + log(C)
        log_cf += log_del

        if abs(log_del) < ϵ; return log_cf end
    end

    error("beta continued fraction did not converge after $MAXIT iterations")
end


end # module
