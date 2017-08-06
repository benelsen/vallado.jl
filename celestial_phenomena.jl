using StaticArrays

include("nutation.jl")
include("time.jl")
include("utils.jl")
include("coordinate_systems.jl")
include("state.jl")

function sun_ecliptic(ut1_jd)
  ut1_jc = julianYears(ut1_jd)

  # Explanatory Supplement to the Astronomical Almanac - 1991
  # Paul Kenneth Seidelmann
  # p 484, eq 9.311-2
  # λ_M = rem2pi(deg2rad(280.460 + 36_000.770 * ut1_jc))
  # M = rem2pi(deg2rad(357.528 + 35_999.050 * ut1_jc))
  # λ_ecl = λ_M + 0.03342306 * sin(M) + 0.00034907 * sin(2M)
  # ϵ = deg2rad(23.4393 - 0.01300 * ut1_jc)

  # Fundamentals of Astrodynamics and Applications - 2013
  # David A. Vallado
  # p 280, alg 29
  # M = rem2pi(deg2rad(357.529_109_2 + 35_999.050_34 * ut1_jc))
  M = as_to_rad_mod( @evalpoly(ut1_jc, 1_287_104.793_048,   129_596_581.048_1,  -0.553_2,  0.000_136, -0.000_011_49) )
  λ_M = rem2pi(deg2rad(280.460_618_4 + 36_000.770_053_61 * ut1_jc))
  λ_ecl = λ_M + deg2rad(1.914_666_471 * sin(M) + 0.019_994_643 * sin(2M))
  r = 1.000_140_612 - 0.016_708_617 * cos(M) - 0.000_139_589 * cos(2M)

  return r, λ_ecl, 0.0
end

function sun_radec(ut1_jd)
  r, λ_ecl, ϕ_ecl = sun_ecliptic(ut1_jd)
  ϵ = obliquity_IAU_2000(julianYears(ut1_jd))
  α, δ = ecliptic_to_radec(λ_ecl, ϕ_ecl, ϵ)
  return r, α, δ
end

function sun_vector(ut1_jd)
  r, α, δ = sun_radec(ut1_jd)
  return radec_to_state(r, α, δ)
end

function JDtoGMST0(jd)
  T = julianYears(jd)
  θGMST = mod(@evalpoly(T, 24_110.548_41, 8_640_184.812_866, 0.093_104, -6.2e-6), 86400) # [s]
  return θGMST / 86400 * 2π
  # θGMST = @evalpoly(T, 4.894961212823059, 230121.67531542323, 6.770713944903336e-6, -4.5087672343186846e-10)
  # return mod2pi(θGMST)
end

function sunrise_sunset(ut1_jd, λ, ϕ, σ_d = dms_to_d(90, 50))
  σ = deg2rad(σ_d)

  ut1_jd_sunrise = floor(ut1_jd + 0.5) - 0.5 +  6/24 - λ/360
  r_sr, α_sr, δ_sr = sun_radec(ut1_jd_sunrise)
  LHA_sunrise = 2π - acos((cos(σ) - sin(δ_sr) * sin(ϕ)) / (cos(δ_sr) * cos(ϕ)))
  GMST_sr = JDtoGMST0(floor(ut1_jd + 0.5) - 0.5 +  6/24)

  println(rad2deg(LHA_sunrise), " ", rad2deg(α_sr), " ", rad2deg(GMST_sr))
  UT_sunrise = mod2pi(LHA_sunrise + α_sr - GMST_sr)

  ut1_jd_sunset  = floor(ut1_jd + 0.5) - 0.5 + 18/24 - λ/360
  r_ss, α_ss, δ_ss = sun_radec(ut1_jd_sunset)
  LHA_sunset = acos((cos(σ) - sin(δ_ss) * sin(ϕ)) / (cos(δ_ss) * cos(ϕ)))
  GMST_ss = JDtoGMST0(floor(ut1_jd + 0.5) - 0.5 + 18/24)

  println(rad2deg(LHA_sunset), " ", rad2deg(α_ss), " ", rad2deg(GMST_ss))
  UT_sunset = mod2pi(LHA_sunset + α_ss - GMST_ss)

  return UT_sunrise, UT_sunset
end

function moon_ecliptic(tdb_jd)
  tdb_jc = julianYears(tdb_jd)

  M_m, M_s, u_sm, D_s = lunisolar_fundamental_arguments(tdb_jc)

  # Meeus 1991:132
  λ_M = deg2rad(218.32 + 481_267.883 * tdb_jc)

  # Green 1988:174
  λ_ecl = λ_M + deg2rad(6.29 * sin(M_m) - 1.27 * sin(M_m - 2D_s) + 0.66 * sin(2D_s) + 0.21 * sin(2M_m) - 0.19 * sin(M_s) - 0.11 * sin(2u_sm))
  ϕ_ecl = deg2rad(5.13 * sin(u_sm) + 0.28 * sin(M_m + u_sm) - 0.28 * sin(u_sm - M_m) - 0.17 * sin(u_sm - 2D_s))
  p = deg2rad(0.9508 + 0.0518 * cos(M_m) + 0.0095 * cos(M_m - 2D_s) + 0.0078 * cos(2D_s) + 0.0028 * cos(2M_m))

  r_m = 6378.137 / sin(p)

  return r_m, mod2pi(λ_ecl), mod2pi(ϕ_ecl)
end

function moon_radec(tdb_jd)
  r, λ_ecl, ϕ_ecl = moon_ecliptic(tdb_jd)
  ϵ = obliquity_IAU_2000(julianYears(tdb_jd))
  α, δ = ecliptic_to_radec(λ_ecl, ϕ_ecl, ϵ)
  return r, mod2pi(α), δ
end

function moon_vector(tdb_jd)
  r, α, δ = moon_radec(tdb_jd)
  return radec_to_state(r, α, δ)
end

function moon_riseset(ut1_jd, λ, ϕ, set = false)
  const h_m = 0.00233

  ut1_jd_0 = floor(ut1_jd + 0.5) - 0.5
  ut1_jd_t = ut1_jd_0 + (set ? 18/24 : 0.0)
  t = (set ? 18/24 : 0.0)

  Δut = 0.0
  i = 0

  while abs(Δut) > (15 / 86400) || i === 0
    i += 1

    r, α, δ = moon_radec(ut1_jd_t)
    θ_GMST = JDtoGMST(ut1_jd_t)
    println("θ_GMST:", rad2deg(θ_GMST))
    println("α:", rad2deg(α))

    GHA_i = θ_GMST - α
    LHA = GHA_i + λ

    if i === 1
      ΔGHA = deg2rad(347.81)
    else
      ΔGHA = (GHA_i - GHA) / Δut
      println("GHA: ", rad2deg(GHA))
    end

    if ΔGHA < 0
      ΔGHA += 2π / Δut
    end
    println("GHA_i: ", rad2deg(GHA_i))
    println("ΔGHA: ", rad2deg(ΔGHA))


    x_i = (h_m - sin(δ) * sin(ϕ)) / (cos(δ) * cos(ϕ))

    println("x_i: ", x_i)

    if set
      LHA_i = acos(x_i)
    else
      LHA_i = 2π - acos(x_i)
    end

    println("LHA: ", rad2deg(LHA))
    println("LHA_i: ", rad2deg(LHA_i))

    Δut = (LHA_i - LHA) / ΔGHA
    println("Δut: ", Δut)

    if Δut < -0.5
      Δut += 2π / ΔGHA
    elseif Δut > 0.5
      Δut -= 2π / ΔGHA
    end
    println("Δut: ", Δut, " ", Δut * 24)

    GHA = GHA_i

    t += Δut
    println("t: ", t, " ", t * 24)

    ut1_jd_t = ut1_jd_0 + t
  end

  return t, ut1_jd_t

end

function moonrise_moonset(ut1_jd)
  ut_mr = moon_riseset(ut1_jd, false)
  ut_ms = moon_riseset(ut1_jd, false)

  return ut_mr, ut_ms
end

function moon_phase(tdb_jd)
  r_m, λ_m, ϕ_m = moon_ecliptic(tdb_jd)
  r_s, λ_s, ϕ_s = sun_ecliptic(tdb_jd)

  phase = λ_s - λ_m

  return phase, (1 - cos(phase)) / 2
end

# Meeus 1991:202-204
# J2000
function jupter_elements(tdb_jc)
  a = @evalpoly(tdb_jc, 5.202_603_191,  0.191_3e-6)
  e = @evalpoly(tdb_jc, 0.048_494_85, 163.244e-6, -471.9e-9, -1.97e-9)
  i = @evalpoly(tdb_jc, 1.303_270, -1.987_2e-3, 33.18e-6, 92e-9)
  Ω = @evalpoly(tdb_jc, 100.464_441, 176.682_8e-3, 903.87e-6, -7.032e-6)
  ω̃ = @evalpoly(tdb_jc, 14.331_309, 215.552_5e-3, 722.52e-6, -4.590e-6)
  λ = @evalpoly(tdb_jc, 34.351_484, 3_034.905_674_6, -85.01e-6, 4e-9)

  return a, e, (deg2rad(i)), (deg2rad(Ω)), (deg2rad(ω̃)), (deg2rad(λ))
end

function planet_ecliptic(tdb_jd)
  tdb_jc = julianYears(tdb_jd)

  a, e, i, Ω, ω̃, λ = jupter_elements(tdb_jc)
  M = λ - ω̃
  println(rad2deg(M))
  ω = ω̃ - Ω
  println(rad2deg(ω))
  p = a * (1 - e^2)
  println(p)
  E = kepler_equation_ell(M, e, 1e-16)
  ν = EBH_to_ν(E, e)
  println(rad2deg(ν))

  elem = StateKepler(p, e, i, Ω, ω, ν, tdb_jc)

  state = elementsToState(elem; μ = 1)

  return state
end

function shadow_geometric(r_s, r)
  if dot(r_s, r) < 0.0
    σ = acos(dot(-r_s, r) / (norm(r_s) * norm(r)))

    h_satellite = norm(r) * cos(σ)
    v_satellite = norm(r) * sin(σ)

    x = R_primary / sin(α_penumbra)
    v_penumbra = tan(α_penumbra) * (x + h_satellite)

    y = R_primary / sin(α_umbra)
    v_umbra = tan(α_umbra) * (y - h_satellite)

    if v_satellite ≤ v_umbra
      return :umbra
    elseif v_satellite ≤ v_penumbra
      return :penumbra
    end
  end
  return :no_shadow
end

function sight(r1, r2; ellipsoidal = true)
  const e = 0.081_819_221_456
  const R = 6378.137

  if ellipsoidal
    scale = SVector(1, 1, 1/sqrt(1 - e^2))
    r1 = scale .* r1
    r2 = scale .* r2
  end

  τ = (norm(r1)^2 - dot(r1, r2)) / (norm(r1)^2 + norm(r2)^2 - 2 * dot(r1, r2))

  if τ < 0 || τ > 1
    return true
  elseif ((1 - τ) * norm(r1)^2 + τ * dot(r1, r2)) ≥ R^2
    return true
  end

  return false
end

function light(r1, ut1_jd)
  const AU = 149597870.700

  r_s = sun_vector(ut1_jd)

  return sight(r1, r_s)
end

function ground_illumination(el, ut1_jd, body = :sun)
  x = 2 * el / π

  if body === :sun
    if el ≥ deg2rad(20)
      L1 = @evalpoly(x,  3.74,   3.97,   -4.07,     1.47)
    elseif el ≥ deg2rad(5)
      L1 = @evalpoly(x,  3.05,  13.28,  -45.98,    64.33)
    elseif el ≥ deg2rad(-0.8)
      L1 = @evalpoly(x,  2.88,  22.26, -207.64,  1034.30)
    elseif el ≥ deg2rad(-5)
      L1 = @evalpoly(x,  2.88,  21.81, -258.11,  -858.36)
    elseif el ≥ deg2rad(-12)
      L1 = @evalpoly(x,  2.70,  12.17, -431.69, -1899.83)
    elseif el ≥ deg2rad(-18)
      L1 = @evalpoly(x, 13.84, 262.72, 1447.42,  2797.93)
    else
      L1 = -3.19064
    end

    return 10^L1

  elseif body === :moon
    if el ≥ deg2rad(20)
      L1 = @evalpoly(x, -1.95,  4.06,   -4.24,    1.56)
    elseif el ≥ deg2rad(5)
      L1 = @evalpoly(x, -2.58, 12.58,  -42.58,   59.06)
    elseif el ≥ deg2rad(-0.8)
      L1 = @evalpoly(x, -2.79, 24.27, -252.95, 1321.29)
    else
      L1 = -10.0
    end

    r_s, α_s, δ_s = sun_radec(ut1_jd)
    r_m, α_m, δ_m = moon_radec(ut1_jd)

    E = acos(sin(δ_s) * sin(δ_m) + cos(δ_s) * cos(δ_m) * cos(α_s - α_m))
    f = π - E
    L2 = -8.68e-3 * f - 2.2e-9 * f^4

    p = asin(6378.137 / r_m)
    L3 = 2 * log10(p / 0.951)

    return 10^(L1 + L2 + L3)
  end
end
