

struct MCConfidenceInterval
    # magic constants
    C₀::Float64
    C₁::Float64
    param::Parameter
    # mutables
    w::Array{Float64, 1}
    M::Array{Float64, 1}
    # Stores the bounds
    upper::Array{Float64, 1}
    lower::Array{Float64, 1}
    uppers::Array{Float64, 1}
    lowers::Array{Float64, 1}
    # Ctor
    function MCConfidenceInterval(parameter::Parameter, nz)
        C₀, C₁ = 0.7164, 31.395
        w = zeros(nz)
        M = zeros(nz)
        w = parameter.ead ./ sum(parameter.ead)
        new(C₀, C₁, parameter, w, M, zeros(nz), zeros(nz), zeros(nz), zeros(nz))
    end
end


function track_ci(ci::MCConfidenceInterval, pnc, μ, σ, p, i)
    ci.M = ci.w .* sum(ci.param.lgc .* pnc; dims=2)
    ρ = sum(@. abs((ci.param.weights - ci.M) / σ)^3 * ci.pnc)
    bound = min(ci.C₀ * ρ, ci.C₁ * ρ / (1 + abs((ci.param.l-μ)/σ)^3))

    ci.lower[i] = p - bound
    ci.upper[i] = p + bound
    ci.lowers[i] = mean(ci.lower[1:i]) + invtcdf(0.95; ν=i) * sqrt(var(ci.lower[1:i])/i)
    ci.uppers[i] = mean(ci.upper[1:i]) + invtcdf(0.95; ν=i) * sqrt(var(ci.upper[1:i])/i)
end
