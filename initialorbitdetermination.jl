include("time.jl")
include("utils.jl")
include("nutation.jl")
include("coordinate_systems.jl")
include("eop.jl")
include("state.jl")
include("lagrange.jl")

using PolynomialRoots
using StaticArrays
using ForwardDiff

function site_track(utc_dt, gd::GeodeticCoordinate, ρ, β, e, ρ̇, β̇, ė; eop = EOP(0, 0, 0, 0, 0, 0, 0))
    x_SEZ, v_SEZ = azel_to_SEZ(ρ, β, e, ρ̇, β̇, ė)
    x_ECEF, v_ECEF = SEZ_to_ECEF(gd, x_SEZ, v_SEZ)

    tt_jd = dateToJD(TAItoTT(UTCtoTAI(utc_dt, 29)))
    ut1_jd = dateToJD(UTCtoUT1(utc_dt, eop.ΔUT1))

    return ITRF_to_GCRF(tt_jd, ut1_jd; x_p = eop.x, y_p = eop.y, ΔX = eop.ΔX, ΔY = eop.ΔY, LOD = eop.ΔLOD)(x_ECEF, v_ECEF)
end

abstract type Observations end

struct ObservationsRhoAzEl{T <: Real} <: Observations
    t :: DateTime
    r :: AbstractArray
    ρ :: Nullable{T}
    β :: T
    e :: T
    ρ̇ :: Nullable{T}
    β̇ :: Nullable{T}
    ė :: Nullable{T}
end

struct ObservationsRhoRaDec{T <: Real} <: Observations
    t :: DateTime
    r :: AbstractArray
    ρ :: Nullable{T}
    α :: T
    δ :: T
    ρ̇ :: Nullable{T}
    α̇ :: Nullable{T}
    δ̇ :: Nullable{T}
end

ObservationsRhoRaDec(t, r, ρ::T, α::T, δ::T, ρ̇::T, α̇::T, δ̇::T) where T <: Real = ObservationsRhoRaDec(t, r, Nullable{T}(ρ), α, δ, Nullable{T}(ρ̇), Nullable{T}(α̇), Nullable{T}(δ̇))
ObservationsRhoRaDec(t, r, ρ::T, α::T, δ::T) where T <: Real = ObservationsRhoRaDec(t, r, Nullable{T}(ρ), α, δ, Nullable{T}(), Nullable{T}(), Nullable{T}())

ObservationsRhoRaDec(t, r, α::T, δ::T, α̇::T, δ̇::T) where T <: Real = ObservationsRhoRaDec(t, r, Nullable{T}(), α, δ, Nullable{T}(), Nullable{T}(α̇), Nullable{T}(δ̇))
ObservationsRhoRaDec(t, r, α::T, δ::T) where T <: Real = ObservationsRhoRaDec(t, r, Nullable{T}(), α, δ, Nullable{T}(), Nullable{T}(), Nullable{T}())

function is_real_positive(z)
    isapprox(imag(z), 0.0, atol=1e-10) && real(z) > 0
end

function is_almost_real(z)
    isapprox(imag(z), 0.0, atol=1e-10)
end

function laplace_method(obs; μ = 398600.4418)

    middle_i = convert(Int, ceil(length(obs) / 2))

    Ls = [radec_to_state(1, o.α, o.δ) for o in obs]
    rs = [o.r for o in obs]
    ts = [dateToJD(o.t) * 86400 for o in obs]

    t_m = ts[middle_i]
    τs = ts - t_m

    L = lagrange(Ls, τs, 0.)
    L̇ = lagrange_d1(Ls, τs, 0.)
    L̈ = lagrange_d2(Ls, τs, 0.)

    r = lagrange(rs, τs, 0.)
    ṙ = lagrange_d1(rs, τs, 0.)
    r̈ = lagrange_d2(rs, τs, 0.)

    D = 2 * det(hcat(L, L̇, L̈))
    D1 = det(hcat(L, L̇, r̈))
    D2 = det(hcat(L, L̇, r))

    C = dot(L, r)

    p1_coeffs = [
    -4 * μ^2 * D2^2 / D^2,
    0., 0.,
    μ * (4 * C * D2 / D - 8 * D1 * D2 / D^2),
    0., 0.,
    (4 * C * D1 / D - 4 * D1^2 / D^2 - norm(r)^2),
    0.,
    1.]

    croots = roots(p1_coeffs, polish = true)
    rroots = map(real, filter(is_real_positive, croots))

    if length(rroots) > 1
        error("more then 1 real positive root found")
    end

    r_mag = rroots[1]
    ρ_mag = -2D1 / D - 2μ / r_mag^3 * D2 / D

    r_m = ρ_mag * L + rs[middle_i]

    D3 = det(hcat(L, r̈, L̈))
    D4 = det(hcat(L, r, L̈))

    ρ̇_mag = - D3 / D - μ / r_mag^3 * D4 / D

    v_m = ρ̇_mag * L + ρ_mag * L̇ + ṙ

    return StateVector(r_m, v_m, JDtoDate(t_m / 86400))
end

function gauss_method(obs; μ = 398600.4418)

    Ls = [radec_to_state(1, o.α, o.δ) for o in obs]
    rs = [o.r for o in obs]
    ts = [dateToJD(o.t) * 86400 for o in obs]

    t_m = ts[2]
    τs = ts - t_m

    a_1 = τs[3] / (τs[3] - τs[1])
    a_3 = - τs[1] / (τs[3] - τs[1])

    a_1u = τs[3] * ((τs[3] - τs[1])^2 - τs[3]^2) / (6 * (τs[3] - τs[1]))
    a_3u = - τs[1] * ((τs[3] - τs[1])^2 - τs[1]^2) / (6 * (τs[3] - τs[1]))

    L_det = det(hcat(Ls...))

    L_inv = @SMatrix [
    (  Ls[2][2] * Ls[3][3] - Ls[3][2] * Ls[2][3]) (- Ls[1][2] * Ls[3][3] + Ls[3][2] * Ls[1][3]) (  Ls[1][2] * Ls[2][3] - Ls[2][2] * Ls[1][3]) ;
    (- Ls[2][1] * Ls[3][3] + Ls[3][1] * Ls[2][3]) (  Ls[1][1] * Ls[3][3] - Ls[3][1] * Ls[1][3]) (- Ls[1][1] * Ls[2][3] + Ls[2][1] * Ls[1][3]) ;
    (  Ls[2][1] * Ls[3][2] - Ls[3][1] * Ls[2][2]) (- Ls[1][1] * Ls[3][2] + Ls[3][1] * Ls[1][2]) (  Ls[1][1] * Ls[2][2] - Ls[2][1] * Ls[1][2]) ]

    L_inv = L_inv.' / L_det

    M = L_inv * hcat(rs...)

    M[2,1]

    d_1 = M[2, 1] * a_1 - M[2, 2] + M[2, 3] * a_3
    d_2 = M[2, 1] * a_1u + M[2, 3] * a_3u

    C = dot(Ls[2], rs[2])

    p1_coeffs = [
    - μ^2 * d_2^2,
    0., 0.,
    -2μ * (C * d_2 + d_1 * d_2),
    0., 0.,
    -(d_1^2 + 2C * d_1 + norm(rs[2])^2),
    0.,
    1.]

    croots = roots(p1_coeffs, polish = true)
    rroots = map(real, filter(is_real_positive, croots))

    if length(rroots) > 1
        error("more then 1 real positive root found")
    end

    r_mag = rroots[1]

    u = μ / r_mag^3

    c = SVector(a_1 + a_1u * u, -1., a_3 + a_3u * u)

    c_ρ = -M * c
    ρ = c_ρ ./ c

    r = ρ .* Ls + rs

    # refine position vectors ?!
    # Gibbs & Herrick-Gibbs

    return r
end

function double_r(obs; t_m = +1., ϵ = 1e-6, max_iter = 100, μ = 398600.4418)
    @assert length(obs) === 3

    Ls = [radec_to_state(1., o.α, o.δ) for o in obs]
    rs = [o.r for o in obs]
    ts = [dateToJD(o.t) * 86400 for o in obs]

    r1 = isnull(obs[1].ρ) ? 12756.274 : sqrt(obs[1].ρ^2 + 2 * obs[1].ρ * (Ls[1] ⋅ rs[1]) + norm(rs[1])^2)
    r2 = isnull(obs[2].ρ) ? 12820.055 : sqrt(obs[2].ρ^2 + 2 * obs[2].ρ * (Ls[2] ⋅ rs[2]) + norm(rs[2])^2)
    r3 = isnull(obs[3].ρ) ? 12000.000 : sqrt(obs[3].ρ^2 + 2 * obs[3].ρ * (Ls[3] ⋅ rs[3]) + norm(rs[3])^2)

    τs = ts - ts[2]

    c = 2Ls .⋅ rs

    calc_F = function (r_mag)
        ρ = (-c + sqrt.(c.^2 - 4(norm.(rs).^2 - (r_mag).^2))) / 2
        r = ρ .* Ls + rs

        W = cross(r[1], r[2]) / (norm(r[1]) * norm(r[2]))

        ρ[3] = -(rs[3] ⋅ W) / (Ls[3] ⋅ W)

        r[3] = ρ[3] .* Ls[3] + rs[3]

        cos_Δν = zeros(MMatrix{3, 3, Float64})
        sin_Δν = zeros(MMatrix{3, 3, Float64})
        Δν = zeros(MMatrix{3, 3, Float64})

        for j in 2:3
            for k in 1:2
                cos_Δν[j, k] = (r[j] ⋅ r[k]) / (norm(r[j]) * norm(r[k]))
                if cos_Δν[j, k] > 1.
                    # warn(j, " ", k, " ", cos_Δν[j, k])
                    cos_Δν[j, k] = 1.0
                end
                sin_Δν[j, k] = t_m * sqrt(1 - cos_Δν[j, k]^2)
                Δν[j, k] = atan2(sin_Δν[j, k], cos_Δν[j, k])
            end
        end

        if Δν[3, 1] > π
            c_1 = norm(r[2]) / norm(r[1]) * sin_Δν[3, 2] / sin_Δν[3, 1]
            c_3 = norm(r[2]) / norm(r[3]) * sin_Δν[2, 1] / sin_Δν[3, 1]
            p = (c_1 * norm(r[1]) + c_3 * norm(r[3]) - norm(r[2])) / (c_1 + c_3 - 1)
        else
            c_1 = norm(r[1]) / norm(r[2]) * sin_Δν[3, 1] / sin_Δν[3, 2]
            c_3 = norm(r[1]) / norm(r[3]) * sin_Δν[2, 1] / sin_Δν[3, 2]
            p = (c_3 * norm(r[3]) - c_1 * norm(r[2]) + norm(r[1])) / (-c_1 + c_3 + 1)
        end

        e_cos_ν = p ./ norm.(r) - 1

        if Δν[2, 1] ≉ π
            e_sin_ν2 = (-cos_Δν[2, 1] * e_cos_ν[2] + e_cos_ν[1]) / sin_Δν[2, 1]
        else
            e_sin_ν2 = ( cos_Δν[3, 2] * e_cos_ν[2] - e_cos_ν[3]) / sin_Δν[3, 2]
        end

        e = sqrt(e_cos_ν[2]^2 + e_sin_ν2^2)

        a = p / (1 - e^2)
        n = sqrt(μ / a^3)

        S = norm(r[2]) / p * sqrt(1 - e^2) * e_sin_ν2
        C = norm(r[2]) / p * (e^2 + e_cos_ν[2])

        ΔE32 = atan2( norm(r[3]) / sqrt(a * p) * sin_Δν[3, 2] - norm(r[3]) / p * (1 - cos_Δν[3, 2]) * S,
        1 - norm(r[2]) * norm(r[3]) / (a * p) * (1 - cos_Δν[3, 2]))

        ΔE21 = atan2( norm(r[1]) / sqrt(a * p) * sin_Δν[2, 1] + norm(r[1]) / p * (1 - cos_Δν[2, 1]) * S,
        1 - norm(r[2]) * norm(r[1]) / (a * p) * (1 - cos_Δν[2, 1]))

        ΔM32 =  ΔE32 + 2S * sin(ΔE32 / 2)^2 - C * sin(ΔE32)
        ΔM12 = -ΔE21 + 2S * sin(ΔE21 / 2)^2 + C * sin(ΔE21)

        local F1 = τs[1] - ΔM12 / n
        local F2 = τs[3] - ΔM32 / n

        return F1, F2, r, a, ΔE32
    end

    for iter in 1:max_iter
        # println("Iteration $iter")

        r_mag = [r1, r2, r3]
        # println("r_mag: ", r_mag)
        F1, F2, r, a, ΔE32 = calc_F(r_mag)
        # println("F1: ", F1, " F2: ", F2)

        dr1 = r1 * 0.005
        r_mag_dr1 = [r1 + dr1, r2, r3]
        F1_dr1, F2_dr1, r, a, ΔE32 = calc_F(r_mag_dr1)
        # println("F1: ", F1_dr1, " F2: ", F2_dr1)

        dF1_dr1 = (F1_dr1 - F1) / dr1
        dF2_dr1 = (F2_dr1 - F2) / dr1
        # println("dr1: ", dr1, " ", dF1_dr1, " ", dF2_dr1)

        dr2 = r2 * 0.005
        r_mag_dr2 = [r1, r2 + dr1, r3]
        F1_dr2, F2_dr2, r, a, ΔE32 = calc_F(r_mag_dr2)
        # println("F1: ", F1_dr2, " F2: ", F2_dr2)

        dF1_dr2 = (F1_dr2 - F1) / dr2
        dF2_dr2 = (F2_dr2 - F2) / dr2
        # println("dr2: ", dr2, " ", dF1_dr2, " ", dF2_dr2)

        Δ = dF1_dr1 * dF2_dr2 - dF2_dr1 * dF1_dr2
        Δ1 = dF2_dr2 * F1 - dF1_dr2 * F2
        Δ2 = dF1_dr1 * F2 - dF2_dr1 * F1

        Δr1 = - Δ1 / Δ
        Δr2 = - Δ2 / Δ

        r1 = r1 + Δr1
        r2 = r2 + Δr2

        if abs(Δr1) < ϵ && abs(Δr2) < ϵ
            break
        end

        if iter === max_iter
            error("Did not converge after $iter iterations")
        end
    end

    F1, F2, r, a, ΔE32 = calc_F([r1, r2, r3])

    Q = sqrt(F1^2 + F2^2)

    f = 1 - a / norm(r[2]) * (1 - cos(ΔE32))
    g = τs[3] - sqrt(a^3 / μ) * (ΔE32 - sin(ΔE32))

    v_2 = (r[3] - f * r[2]) / g

    return StateVector(r[2], v_2, obs[2].t)
end

function gibbs_method(r, t_2; μ = 398600.4418)
    @assert length(r) === 3

    Z_12 = cross(r[1], r[2])
    Z_23 = cross(r[2], r[3])
    Z_31 = cross(r[3], r[1])

    α_cop = asin((Z_23 ⋅ r[1]) / (norm(Z_23) * norm(r[1])))
    if α_cop > deg2rad(5.)
        warn("Vectors not coplanar (α_cop = ", rad2deg(α_cop), "°)")
    end

    α_12 = acos((r[1] ⋅ r[2]) / (norm(r[1]) * norm(r[2])))
    α_23 = acos((r[2] ⋅ r[3]) / (norm(r[2]) * norm(r[3])))
    if α_12 < deg2rad(1.) || α_23 < deg2rad(1.)
        warn("Vectors very close together (α_12 = ", rad2deg(α_12), "°, α_23 = ", rad2deg(α_23), "°)")
    end

    N = norm(r[1]) * Z_23 + norm(r[2]) * Z_31 + norm(r[3]) * Z_12
    D = Z_12 + Z_23 + Z_31
    S = (norm(r[2]) - norm(r[3])) * r[1] + (norm(r[3]) - norm(r[1])) * r[2] + (norm(r[1]) - norm(r[2])) * r[3]

    B = cross(D, r[2])

    L_g = sqrt(μ / norm(N) / norm(D))

    v_2 = L_g / norm(r[2]) * B + L_g * S

    return StateVector(r[2], v_2, t_2)
end

function herrick_gibbs_method(r, t; μ = 398600.4418)
    @assert length(r) === 3

    Δt_31 = (t[3] - t[1]) * 86400
    Δt_32 = (t[3] - t[2]) * 86400
    Δt_21 = (t[2] - t[1]) * 86400

    Z_23 = cross(r[2], r[3])

    α_cop = asin((Z_23 ⋅ r[1]) / (norm(Z_23) * norm(r[1])))
    if α_cop > deg2rad(5.)
        warn("Vectors not coplanar (α_cop = ", rad2deg(α_cop), "°)")
    end

    α_12 = acos((r[1] ⋅ r[2]) / (norm(r[1]) * norm(r[2])))
    α_23 = acos((r[2] ⋅ r[3]) / (norm(r[2]) * norm(r[3])))
    if α_12 > deg2rad(5.) || α_23 > deg2rad(5.)
        warn("Vectors not close together (α_12 = ", rad2deg(α_12), "°, α_23 = ", rad2deg(α_23), "°)")
    end

    v_2 = - Δt_32 * (1 / (Δt_21 * Δt_31) + μ / (12 * norm(r[1])^3) ) * r[1] + (Δt_32 - Δt_21) * (1 / (Δt_21 * Δt_32) + μ / (12 * norm(r[2])^3) ) * r[2] + Δt_21 * (1 / (Δt_32 * Δt_31) + μ / (12 * norm(r[3])^3) ) * r[3]

    return StateVector(r[2], v_2, t[2])
end

function lambert_minumum_energy(r0, r1; μ = 398600.4418)

    cos_Δν = (r0 ⋅ r1) / (norm(r0) * norm(r1))
    Δν = acos(cos_Δν)

    c = sqrt(norm(r0)^2 + norm(r1)^2 - 2 * norm(r0) * norm(r1) * cos_Δν)
    s = (norm(r0) + norm(r1) + c) / 2

    a_min = s/2
    p_min = norm(r0) * norm(r1) / c * (1 - cos_Δν)
    e_min = sqrt(1 - 2 * p_min / s)

    α_min = π
    β_min = 2 * asin(sqrt((s - c) / s))

    t_min_amin_1 = sqrt(a_min^3 / μ) * (α_min - (β_min - sin(β_min)))
    t_min_amin_2 = sqrt(a_min^3 / μ) * (α_min + (β_min - sin(β_min)))

    t_min_abs = sqrt(2 / μ) * (s^1.5 - (s - c)^1.5) / 3

    v_0 = sqrt(μ * p_min) / (norm(r0) * norm(r1) * sin(Δν)) * (r1 - (1 - norm(r1) / p_min * (1 - cos_Δν)) * r0)

    return v_0, t_min_amin_1, t_min_amin_2, t_min_abs, a_min, p_min, e_min
end

function newton_raphson_method(fn, y0 = 1.; max_iter = 100, ϵ = 1e-15)
    yn = y0
    fi = 1.
    for i in 1:max_iter
        fio = fi
        fi = fn(yn)
        Δfi = fi - fio

        dfi = ForwardDiff.derivative(fn, yn)
        Δyn = - fi ./ dfi
        yn = yn + Δyn

        if abs(Δyn) < ϵ
            # println("Iterations: $i y: $yn Δy: $Δyn Δf: $Δfi f: $fi df: $dfi")
            # break
            return yn
        end

        if i === max_iter
            error("Failed to converge after $i iterations, y: $yn Δy: $Δyn Δf: $Δfi f: $fi df: $dfi")
        end
    end
end

function lambert_gauss(r0, r1, Δt; t_m = +1., elliptic = true, ϵ = 1e-9, max_iter = 100, μ = 398600.4418)

    cos_Δν = (r0 ⋅ r1) / (norm(r0) * norm(r1))
    sin_Δν = t_m * sqrt(1 - cos_Δν^2)
    Δν = atan2(sin_Δν, cos_Δν)

    l = (norm(r0) + norm(r1)) / (4 * sqrt(norm(r0) * norm(r1)) * cos(Δν / 2)) - 0.5
    m = μ * Δt^2 / (2 * sqrt(norm(r0) * norm(r1)) * cos(Δν / 2))^3

    println("l: $l, m: $m")

    function y_fn(y)
        x_1 = m / y^2 - l
        g = 2 * asin(sqrt(x_1))
        x_2 = (2g - sin(2g)) / sin(g)^3

        return y - 1. - x_2 * (l + x_1)
    end

    y = newton_raphson_method(y_fn, 1.)

    if elliptic
        cos_ΔE_2 = 1 - 2x_1
        p = (norm(r0) * norm(r1) * (1 - cos_Δν)) / (norm(r0) + norm(r1) - 2 * sqrt(norm(r0) * norm(r1)) * cos(Δν / 2) * cos_ΔE_2)
    else
        cosh_ΔH_2 = 1 - 2x_1
        p = (norm(r0) * norm(r1) * (1 - cos_Δν)) / (norm(r0) + norm(r1) - 2 * sqrt(norm(r0) * norm(r1)) * cos(Δν / 2) * cosh_ΔH_2)
    end

    f = 1 - norm(r1) / p * (1 - cos_Δν)
    g = norm(r0) * norm(r1) * sin_Δν / sqrt(μ * p)

    ḟ = sqrt(1 / p) * tan(Δν / 2) * ((1 - cos_Δν) / p - 1/norm(r1) - 1/norm(r0))
    ġ = 1 - norm(r0) / p * (1 - cos_Δν)

    v0 = (r1 - f * r0) / g
    v1 = (ġ * r1 - r0) / g

    return v0, v1
end

function lambert_univeral(r0, r1, Δt; t_m = +1., ψ_n = 0., ψ_up = 4π^2, ψ_low = -4π^2, ϵ = 1e-6, max_iter = 500, μ = 398600.4418)

    cos_Δν = (r0 ⋅ r1) / (norm(r0) * norm(r1))
    sin_Δν = t_m * sqrt(1 - cos_Δν^2)
    Δν = atan2(sin_Δν, cos_Δν)

    A = t_m * sqrt(norm(r0) * norm(r1) * (1 + cos_Δν))

    if isapprox(A, 0., atol = 1e-12)
        error("A = 0: Can’t calculate the orbit.")
    end

    c_2 = 1/2
    c_3 = 1/6

    # println(A)

    local y_n
    for iter in 1:max_iter
        y_n = norm(r0) + norm(r1) + A * (ψ_n * c_3 - 1) / sqrt(c_2)

        if A > 0. && y_n < 0.
            ψ_low /= 2
            if ψ_low > ψ_up
                ψ_low, ψ_up = ψ_up, ψ_low
            end
            ψ_n = (ψ_up + ψ_low) / 2
            # println(y_n, " ", ψ_n, " ", ψ_up, " ", ψ_low, " readj")
            continue
        end

        # println(y_n, " ", ψ_n, " ", ψ_up, " ", ψ_low, "")

        χ_n = sqrt(y_n / c_2)

        Δt_n = (χ_n^3 * c_3 + A * sqrt(y_n)) / sqrt(μ)

        if Δt_n ≤ Δt
            ψ_low = ψ_n
        else
            ψ_up = ψ_n
        end

        # println(Δt_n, " ", ψ_n, " ", ψ_up, " ", ψ_low, " Δt adj")

        ψ_n = (ψ_up + ψ_low) / 2

        c_2, c_3 = ψ_to_c2c3(ψ_n)

        dt = Δt - Δt_n
        if abs(dt) < ϵ || ψ_up ≈ ψ_low
            break
        end

        if iter === max_iter
            error("Did not converge after $iter iterations (dt = $dt)")
        end
    end

    f = 1 - y_n / norm(r0)
    g = A * sqrt(y_n / μ)
    ġ = 1 - y_n / norm(r1)

    v0 = (r1 - f * r0) / g
    v1 = (ġ * r1 - r0) / g

    return v0, v1
end

function lambert_batin_ξ(x, n_max = 5)
    η = x / (sqrt(1 + x) + x)^2
    s = 1.0
    for n in n_max:-1:4
        c_η = n^2 / ((2n)^2 - 1)
        s = 1 + c_η * η / s
    end
    return 8 * (sqrt(1 + x) + 1) / (3 + 1 / (5 + η + s))
end

function lambert_batin(r0, r1, Δt; t_m = +1., elliptic = true, ϵ = 1e-12, max_iter = 100, μ = 398600.4418)

    cos_Δν = (r0 ⋅ r1) / (norm(r0) * norm(r1))
    sin_Δν = t_m * sqrt(1 - cos_Δν^2)
    Δν = atan2(sin_Δν, cos_Δν)

    c = sqrt(norm(r0)^2 + norm(r1)^2 - 2 * norm(r0) * norm(r1) * cos_Δν)
    s = (norm(r0) + norm(r1) + c) / 2
    ϵ_ = (norm(r1)  - norm(r0)) / norm(r0)

    tan2_2w = ϵ_^2 / 4 / (sqrt(norm(r1) / norm(r0)) + norm(r1) / norm(r0) * (2 + sqrt(norm(r1) / norm(r0))))

    r_op = sqrt(norm(r0) * norm(r1)) * (cos(Δν/4)^2 + tan2_2w)

    if 0 <= Δν <= π
        l = (sin(Δν/4)^2 + tan2_2w) / (sin(Δν/4)^2 + tan2_2w + cos(Δν/2))
    else
        l = (cos(Δν/4)^2 + tan2_2w - cos(Δν/2)) / (cos(Δν/4)^2 + tan2_2w)
    end

    m = μ * Δt^2 / 8r_op^3

    x = elliptic ? l : 0.

    local y
    for iter in 1:max_iter
        ξ_ = lambert_batin_ξ(x, 4)
        h1 = (l + x)^2 / (1 + 2x + l) * (1 + 3x + ξ_) / (4x + ξ_ * (3 + x))
        h2 = m / (1 + 2x + l) * (x - l + ξ_) / (4x + ξ_ * (3 + x))

        croots = roots([-h2, 0., -1. - h1, 1], polish = true)
        rroots = map(real, filter(is_almost_real, croots))

        if length(rroots) > 1
            error("more then 1 real positive root found")
        end

        y = rroots[1]
        x_ = x
        x = sqrt(((1 - l) / 2)^2 + m / y^2) - (1 + l) / 2

        if abs(x - x_) < ϵ
            break
        end

        if iter === max_iter
            error("Did not converge after $iter iterations")
        end
    end

    a = μ * Δt^2 / (16 * r_op^2 * x * y^2)

    if a > 0
        α_min = π
        β_min = 2 * asin(sqrt((s - c) / 2a))

        if Δν > π
            β_min = -β_min
        end

        a_min = s / 2
        t_min = sqrt(a_min^3 / μ) * (α_min - (β_min - sin(β_min)))

        α_e = 2 * asin(sqrt(s / 2a))
        β_e = 2 * asin(sqrt((s - c) / 2a))

        if Δt > t_min
            α_e = 2π - α_e
        end

        ΔE = α_e - β_e

        f = 1 - a / norm(r0) * (1 - cos(ΔE))
        g = Δt - sqrt(a^3 / μ) * (ΔE - sin(ΔE))
        ġ = 1 - a / norm(r1) * (1 - cos(ΔE))

    else
        α_h = 2 * asinh(sqrt(s / -2a))
        β_h = 2 * asinh(sqrt((s - c) / -2a))

        ΔH = α_h - β_h

        f = 1 - a / norm(r0) * (1 - cosh(ΔH))
        g = Δt - sqrt(-a^3 / μ) * (sinh(ΔH) - ΔH)
        ġ = 1 - a / norm(r1) * (1 - cosh(ΔH))
    end

    v0 = (r1 - f * r0) / g
    v1 = (ġ * r1 - r0) / g
    return v0, v1
end

function hit_central_body(r_i, r_t, v_1, v_2; R = 6378.136)
    if r_i ⋅ v_1 < 0. && r_t ⋅ v_2 > 0

        a = 1 / (2 / norm(r_i) - norm(v_1)^2 / μ)

        h_t = norm(cross(r_i, v_1))
        p = h_t^2 / μ

        e = sqrt((a - p) / a)
        r_p = a * (1 - e)

        if r_p ≤ R
            return true
        end
    end
    return false
end

function target(r_i, r_t, v_i, v_t, Δt_transfer, Δt_wait = 0.; way = :both, prop = true)
    s_i2 = propagate_kepler(StateVector(r_i, v_i, 0.), Δt_wait)
    s_t2 = propagate_kepler(StateVector(r_t, v_t, 0.), Δt_wait + Δt_transfer * (prop ? 1 : 0))

    if way === :short
        v_1, v_2 = lambert_univeral(s_i2.r, s_t2.r, Δt_transfer; t_m = +1.)
        Δv1 = v_1 - s_i2.v
        Δv2 = s_t2.v - v_2
        Δv_total = norm(Δv1) + norm(Δv2)

    elseif way === :long
        v_1, v_2 = lambert_univeral(s_i2.r, s_t2.r, Δt_transfer; t_m = -1.)
        Δv1 = v_1 - s_i2.v
        Δv2 = s_t2.v - v_2
        Δv_total = norm(Δv1) + norm(Δv2)

    elseif way === :both
        v_1_s, v_2_s = lambert_univeral(s_i2.r, s_t2.r, Δt_transfer; t_m = +1.)
        Δv1_s = v_1_s - s_i2.v
        Δv2_s = s_t2.v - v_2_s
        Δv_s = norm(Δv1_s) + norm(Δv2_s)

        v_1_l, v_2_l = lambert_univeral(s_i2.r, s_t2.r, Δt_transfer; t_m = -1.)
        Δv1_l = v_1_l - s_i2.v
        Δv2_l = s_t2.v - v_2_l
        Δv_l = norm(Δv1_l) + norm(Δv2_l)

        if Δv_s < Δv_l
            Δv1, Δv2 = Δv1_s, Δv2_s
            Δv_total = Δv_s
            v_1, v_2 = v_1_s, v_2_s
            way = :short
        else
            Δv1, Δv2 = Δv1_l, Δv2_l
            Δv_total = Δv_l
            v_1, v_2 = v_1_l, v_2_l
            way = :long
        end
    end

    hit = hit_central_body(s_i2.r, s_t2.r, v_1, v_2)

    return Δv1, Δv2, Δv_total, hit, way
end
