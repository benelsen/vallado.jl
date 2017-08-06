using StaticArrays

include("utils.jl")
include("time.jl")
include("matrix.jl")

include("nutation_data.jl")

function GCRF_to_ITRF(tt_jd::Real, ut1_jd::Real; x_p::Real = 0.0, y_p::Real = 0.0, ΔX::Real = 0.0, ΔY::Real = 0.0, LOD::Real = 0.0)
  tt_jc = julianYears(tt_jd)

  W = polar_motion(tt_jc, x_p, y_p)
  R = CIRS_to_TIRS_CIO(ut1_jd)
  Q = precession_nutation_IAU_2000_2006_CIO(tt_jc; ΔX = ΔX, ΔY = ΔY)

  T_GCRF_to_TIRS = R * transpose(Q)

  ωE = SVector{3}(0, 0, (7.292_115_146_706_979e-5 * (1 - LOD/86400)))

  function _transform(x::AbstractArray{Float64,1})
    return W * T_GCRF_to_TIRS * x
  end

  function _transform(x::AbstractArray{Float64,1}, v::AbstractArray{Float64,1})
    x_TIRS = T_GCRF_to_TIRS * x
    v_TIRS = (T_GCRF_to_TIRS * v - cross(ωE, x_TIRS))
    return (W * x_TIRS, W * v_TIRS)
  end

  function _transform(x::AbstractArray{Float64,1}, v::AbstractArray{Float64,1}, a::AbstractArray{Float64,1})
    x_TIRS = T_GCRF_to_TIRS * x
    v_TIRS = (T_GCRF_to_TIRS * v - cross(ωE, x_TIRS))
    a_TIRS = (T_GCRF_to_TIRS * a - cross(cross(ωE, ωE), x_TIRS) - 2 * cross(ωE, v_TIRS))

    return (W * x_TIRS, W * v_TIRS, W * a_TIRS)
  end

  return _transform
end

function ITRF_to_GCRF(tt_jd::Real, ut1_jd::Real; x_p::Real = 0.0, y_p::Real = 0.0, ΔX::Real = 0.0, ΔY::Real = 0.0, LOD::Real = 0.0)
  tt_jc = julianYears(tt_jd)

  W = polar_motion(tt_jc, x_p, y_p)
  R = CIRS_to_TIRS_CIO(ut1_jd)
  Q = precession_nutation_IAU_2000_2006_CIO(tt_jc; ΔX = ΔX, ΔY = ΔY)

  T_TIRS_to_GCRF = Q * transpose(R)

  ωE = SVector{3}(0, 0, (7.292_115_146_706_979e-5 * (1 - LOD/86400)))

  function _transform(x)
    return T_TIRS_to_GCRF * transpose(W) * x
  end

  function _transform(x, v)
    x_TIRS = transpose(W) * x
    v_TIRS = transpose(W) * v
    return (
      T_TIRS_to_GCRF * x_TIRS,
      T_TIRS_to_GCRF * (v_TIRS + cross(ωE, x_TIRS))
    )
  end

  function _transform(x, v, a)
    x_TIRS = transpose(W) * x
    v_TIRS = transpose(W) * v
    a_TIRS = transpose(W) * a

    return (
      T_TIRS_to_GCRF * x_TIRS,
      T_TIRS_to_GCRF * (v_TIRS + cross(ωE, x_TIRS)),
      T_TIRS_to_GCRF * (a_TIRS + cross(cross(ωE, ωE), x_TIRS) + 2 * cross(ωE, v_TIRS))
    )
  end

  return _transform
end

function J2000_to_ITRF(tt_jd::Real, ut1_jd::Real; x_p::Real = 0.0, y_p::Real = 0.0, Δψ::Real = 0.0, Δϵ::Real = 0.0, LOD::Real = 0.0)
  tt_jc = julianYears(tt_jd)
  ut1_jc = julianYears(ut1_jd)

  W = polar_motion_IAU_1980(x_p, y_p)
  N, P, R, Ṙ = precession_nutation_IAU_1980(tt_jc, ut1_jc; gcrf_correction_Δψ = Δψ, gcrf_correction_Δϵ = Δϵ)

  function _transform(x)
    return W * R * N * P * x
  end

  function _transform(x, v)
    return (
      W * R * N * P * x,
      W * (Ṙ * N * P * x + R * N * P * v)
    )
  end

  return _transform
end

function ITRF_to_J2000(tt_jd::Real, ut1_jd::Real; x_p::Real = 0.0, y_p::Real = 0.0, Δψ::Real = 0.0, Δϵ::Real = 0.0, LOD::Real = 0.0)
  tt_jc = julianYears(tt_jd)
  ut1_jc = julianYears(ut1_jd)

  W = polar_motion_IAU_1980(x_p, y_p)
  N, P, R, Ṙ = precession_nutation_IAU_1980(tt_jc, ut1_jc; gcrf_correction_Δψ = 0.0, gcrf_correction_Δϵ = 0.0, LOD = LOD)

  function _transform(x)
    return P.' * N.' * R.' * W.' * x
  end

  function _transform(x, v)
    return (
      P.' * N.' * R.' * W.' * x,
      P.' * N.' * (Ṙ.' * W.' * x + R.' * W.' * v)
    )
  end

  # function _transform(x, v, a)
  #   x_TIRS = transpose(W) * x
  #   v_TIRS = transpose(W) * v
  #   a_TIRS = transpose(W) * a
  #
  #   return (
  #     T_TIRS_to_GCRF * x_TIRS,
  #     T_TIRS_to_GCRF * (v_TIRS + cross(ωE, x_TIRS)),
  #     T_TIRS_to_GCRF * (a_TIRS + cross(cross(ωE, ωE), x_TIRS) + 2 * cross(ωE, v_TIRS))
  #   )
  # end

  return _transform
end

as_to_rad(x) = deg2rad(x / 3600.0)
as_to_rad_mod(x) = (deg2rad(rem(x, 1296000) / 3600.0))
rad_to_as(x) = rad2deg(x) * 3600.0

μas_to_rad(x) = deg2rad(x / 3600e6)
rad_to_μas(x) = rad2deg(x) * 3600e6

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.5.7 - Equation 5.32
function GMST_IAU_2006(ut1_jd, tt_jc)
  era = earth_rotation_angle(ut1_jd)
  eo = as_to_rad( @evalpoly(tt_jc, 0.014_506, 4_612.156_534, 1.391_581_7, -0.000_000_44, -0.000_029_956, -0.000_000_036_8) )
  return rem2pi(era + eo)
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.7.2 - Equation 5.43
# Simon et al. (1994)
function lunisolar_fundamental_arguments(tdb_jc)
  # Mean Anomaly of the Moon
  M_m  = as_to_rad_mod( @evalpoly(tdb_jc,   485_868.249_036, 1_717_915_923.217_8,  31.879_2,  0.051_635, -0.000_244_70) )

  # Mean Anomaly of the Sun
  M_s  = as_to_rad_mod( @evalpoly(tdb_jc, 1_287_104.793_048,   129_596_581.048_1,  -0.553_2,  0.000_136, -0.000_011_49) )

  # Mean Argument of Latitude of the Moon
  u_sm = as_to_rad_mod( @evalpoly(tdb_jc,   335_779.526_232, 1_739_527_262.847_8, -12.751_2, -0.001_037,  0.000_004_17) )

  # Mean Elongation of the Moon from the Sun
  D_s  = as_to_rad_mod( @evalpoly(tdb_jc, 1_072_260.703_692, 1_602_961_601.209_0,  -6.370_6,  0.006_593, -0.000_031_69) )

  # Mean Longitude of the ascending Node of the Moon
  Ω_m  = as_to_rad_mod( @evalpoly(tdb_jc,   450_160.398_036,    -6_962_890.543_1,   7.472_2,  0.007_702, -0.000_059_39) )

  return @SVector [M_m, M_s, u_sm, D_s, Ω_m]
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.7.3 - Equation 5.44
function planetary_fundamental_arguments(tdb_jc)
  # Mean Heliocentric Longitudes of the Planets
  # Souchay et al. (1999)
  λ_M1 = @evalpoly(tdb_jc, 4.402_608_842, 2608.790_314_157_4) # Mercury
  λ_M2 = @evalpoly(tdb_jc, 3.176_146_697, 1021.328_554_621_1) # Venus
  λ_M3 = @evalpoly(tdb_jc, 1.753_470_314,  628.307_584_999_1) # Earth
  λ_M4 = @evalpoly(tdb_jc, 6.203_480_913,  334.061_242_670_0) # Mars
  λ_M5 = @evalpoly(tdb_jc, 0.599_546_497,   52.969_096_264_1) # Jupiter
  λ_M6 = @evalpoly(tdb_jc, 0.874_016_757,   21.329_910_496_0) # Saturn
  λ_M7 = @evalpoly(tdb_jc, 5.481_293_872,    7.478_159_856_7) # Uranus
  λ_M8 = @evalpoly(tdb_jc, 5.311_886_287,    3.813_303_563_8) # Neptune

  # General Precession
  # Kinoshita and Souchay (1990)
  p_λ  = @evalpoly(tdb_jc, 0.0          ,    0.024_381_75   , 0.000_005_386_91)

  return @SVector [λ_M1, λ_M2, λ_M3, λ_M4, λ_M5, λ_M6, λ_M7, λ_M8, p_λ]
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.5.3 - Equation 5.44
# Capitaine et al. (2000)
function earth_rotation_angle(ut1_jd)
  Tu = ut1_jd - 2_451_545.0
  frac_of_day = rem(ut1_jd, 1.0)
  return 2π * rem(0.779_057_273_264_0 + 0.002_737_811_911_354_48 * Tu + frac_of_day, 1.0)
end

# Slightly less accurate:
# function earth_rotation_angle(ut1_jd)
#   Tu = ut1_jd - 2_451_545.0
#   return 2π * rem(0.779_057_273_264_0 + 1.002_737_811_911_354_48 * Tu, 1.0)
# end

function sum_harmonic_components(tt_jc, coeffs, coeffs0, idx, fund_args)

  for j in 0:(size(idx, 1) - 2)
    X_j = 0.0

    for i in ( idx[j+2] - 1 ):-1:idx[j+1]
      a_i = 0.0
      for k in 18:-1:5
        a_i += coeffs[i, k] * fund_args[k-4]
      end
      X_j += coeffs[i, 3] * sin(a_i) + coeffs[i, 4] * cos(a_i)
    end

    coeffs0[j+1] += X_j
  end

  if size(coeffs0, 1) == 6
    X_ = @evalpoly(tt_jc, coeffs0[1], coeffs0[2], coeffs0[3], coeffs0[4], coeffs0[5], coeffs0[6])
  elseif size(coeffs0, 1) == 2
    X_ = @evalpoly(tt_jc, coeffs0[1], coeffs0[2])
  else
    throw(error("Size of coeffs0"))
  end

  return μas_to_rad(X_)
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.5.4 - Equation 5.16
# Capitaine et al. (2003) "Expressions for IAU 2000 precession quantities" - Eq 49-51
# What about changes from Capitaine et al. (2006) "Improvement of the IAU?2000 precession model"?
function precession_nutation_IAU_2000_2006_CIO_XYs(tt_jc; ΔX = 0.0, ΔY = 0.0)
  # P03
  X_poly_coeff = @MVector [-16_617.0, +2_004_191_898.00,    -429_782.90, -198_618.34,     +7.578,  +5.928_5]
  Y_poly_coeff = @MVector [ -6_951.0,        -25_896.00, -22_407_274.70,   +1_900.59, +1_112.526,  +0.135_8]
  s_poly_coeff = @MVector [    +94.0,         +3_808.65,        -122.68,  -72_574.11,    +27.980, +15.620_0]

  # P03_rev1
  # X_poly_coeff = @MVector [-16_617.0, +2_004_191_804.00,    -429_755.80, -198_618.29,     +7.575,  +5.928_5]
  # Y_poly_coeff = @MVector [ -6_951.0,        -24_867.00, -22_407_272.70,   +1_900.26, +1_112.525,  +0.135_8]

  # P03_rev2
  # X_poly_coeff = @MVector [-16_617.0, +2_004_192_130.00,    -429_775.20, -198_618.39,     +7.576,  +5.928_5]
  # Y_poly_coeff = @MVector [ -6_951.0,        -25_817.00, -22_407_280.10,   +1_900.46, +1_112.526,  +0.135_8]

  X_idx = @SVector [1, 1307, 1560, 1596, 1600, 1601]
  Y_idx = @SVector [1, 963, 1240, 1270, 1275, 1276]
  s_idx = @SVector [1, 34, 37, 62, 66, 67]


  fund_args = vcat(lunisolar_fundamental_arguments(tt_jc), planetary_fundamental_arguments(tt_jc))

  X = sum_harmonic_components(tt_jc, nutation_2000A_X_coeffs, X_poly_coeff, X_idx, fund_args) + ΔX
  Y = sum_harmonic_components(tt_jc, nutation_2000A_Y_coeffs, Y_poly_coeff, Y_idx, fund_args) + ΔY
  s = sum_harmonic_components(tt_jc, nutation_2000A_s_coeffs, s_poly_coeff, s_idx, fund_args) - X * Y / 2

  return X, Y, s
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.4.4 - Equation 5.10
# Kaplan (2005) - Eq 6.18 - page 65
# CIRS => GCRF
function precession_nutation_IAU_2000_2006_CIO(tt_jc; ΔX = 0.0, ΔY = 0.0)
  X, Y, s = precession_nutation_IAU_2000_2006_CIO_XYs(tt_jc; ΔX = ΔX, ΔY = ΔY)

  # d = atan(sqrt((X^2 + Y^2) / (1 - X^2 - Y^2)))
  # a = 1 / (1 + cos(d))
  a = 1/2 + 1/8 * (X^2 + Y^2)

  PN = @SMatrix [ (1 - a * X^2) (-a * X * Y)  (X)                   ;
                  (-a * X * Y)  (1 - a * Y^2) (Y)                   ;
                  (-X)          (-Y)          (1 - a * (X^2 + Y^2)) ]

  return PN * rot3(s)
end


# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.6.4 - Equation 5.40
# Capitaine, Wallace, Chapront (2003)
function obliquity_IAU_2000(tt_jc)
  return as_to_rad( @evalpoly(tt_jc, 84_381.406, -46.836_769, -0.000_183_1, 0.002_003_40, -0.576e-6, -43.4e-9) )
end


# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.6.1 - Equation 5.35
function nutation_IAU_2000A_R06_EQU_ang(tt_jc)
  Δψ_poly_coeff = [0.0, 0.0]
  Δψ_idx = [1, 1321, 1358]

  Δϵ_poly_coeff = [0.0, 0.0]
  Δϵ_idx = [1, 1037, 1056]

  fund_args = vcat(lunisolar_fundamental_arguments(tt_jc), planetary_fundamental_arguments(tt_jc))

  Δψ = sum_harmonic_components(tt_jc, nutation_2000_R06_Dpsi_coeffs, Δψ_poly_coeff, Δψ_idx, fund_args)
  Δϵ = sum_harmonic_components(tt_jc, nutation_2000_R06_Deps_coeffs, Δϵ_poly_coeff, Δϵ_idx, fund_args)

  return Δψ, Δϵ
end

# N GCRF => MOD
function nutation_IAU_2000A_R06_EQU(Δψ, Δϵ, ϵ_A)
  return rot1(- ϵ_A - Δϵ) * rot3(-Δψ) * rot1(ϵ_A)
end

# N GCRF => MOD
function nutation_IAU_2000A_R06_EQU(tt_jc)
  Δψ, Δϵ = nutation_IAU_2000A_R06_EQU_ang(tt_jc)
  ϵ_A = obliquity_IAU_2000(tt_jc)
  return rot1(-ϵ_A - Δϵ) * rot3(-Δψ) * rot1(ϵ_A)
end


# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.6.4 - Equation 5.39-40
# Capitaine et al. (2003) "Expressions for IAU 2000 precession quantities" - Eq 37/39
function precession_2006_EQU_ang(tt_jc)
  ϵ_0 = 84_381.406

  ψ_a = @evalpoly(tt_jc, 0,  5_038.481_507, -1.079_006_9, -0.001_140_45,  0.132_851e-3, -95.1e-9)
  ω_a = @evalpoly(tt_jc, ϵ_0,   -0.025_754,  0.051_262_3, -0.007_725_03, -4.67e-7,        3.337e-7)
  χ_a = @evalpoly(tt_jc, 0,     10.556_403, -2.381_429_2, -0.001_211_97,  1.706_63e-4,   -5.60e-8)

  return as_to_rad(ψ_a), as_to_rad(ω_a), as_to_rad(χ_a)
end

function precession_2006_EQU_rev1_ang(tt_jc)
  ϵ_0 = 84_381.406

  ψ_a = @evalpoly(tt_jc, 0,  5_038.481_270, -1.078_996_9, -0.001_140_38,  0.132_851e-3, -95.1e-9)
  ω_a = @evalpoly(tt_jc, ϵ_0,   -0.024_734,  0.051_262_2, -0.007_725_01, -4.67e-7,        3.337e-7)
  χ_a = @evalpoly(tt_jc, 0,     10.556_403, -2.381_429_2, -0.001_211_97,  1.706_63e-4,   -5.60e-8)

  return as_to_rad(ψ_a), as_to_rad(ω_a), as_to_rad(χ_a)
end

function precession_2006_EQU_rev2_ang(tt_jc)
  ϵ_0 = 84_381.406

  ψ_a = @evalpoly(tt_jc, 0,  5_038.482_090, -1.078_992_1, -0.001_140_40,  0.132_851e-3, -95.1e-9)
  ω_a = @evalpoly(tt_jc, ϵ_0,   -0.023_675,  0.051_262_2, -0.007_725_01, -4.67e-7,        3.337e-7)
  χ_a = @evalpoly(tt_jc, 0,     10.556_403, -2.381_429_2, -0.001_211_97,  1.706_63e-4,   -5.60e-8)

  return as_to_rad(ψ_a), as_to_rad(ω_a), as_to_rad(χ_a)
end

# P MOD => GCRF
function precession_2006_EQU(ψ_a, ω_a, χ_a, ϵ_0)
  return rot3(χ_a) * rot1(-ω_a) * rot3(-ψ_a) * rot1(ϵ_0)
end

# P MOD => GCRF
function precession_2006_EQU(tt_jc)
  ϵ_0 = as_to_rad(84_381.406)

  ψ_a, ω_a, χ_a = precession_2006_EQU_ang(tt_jc)
  return rot3(χ_a) * rot1(-ω_a) * rot3(-ψ_a) * rot1(ϵ_0)
end


# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.6.4 - Equation 5.40
# Includes bias
function precession_2006_EQU_FW_ang(tt_jc)
  γ = @evalpoly(tt_jc,     -0.052_928,   10.556_378, 0.493_204_4, -312.38e-6,  -2.788e-6,  2.60e-8)
  φ = @evalpoly(tt_jc, 84_381.412_819,  -46.811_016, 0.051_126_8,  532.89e-6,  -0.440e-6, -1.76e-8)
  ψ = @evalpoly(tt_jc,     -0.041_775, 5038.481_484, 1.558_417_5, -185.22e-6, -26.452e-6, -1.48e-8)

  return as_to_rad(γ), as_to_rad(φ), as_to_rad(ψ)
end

# P MOD => GCRF
function precession_2006_EQU_FW(γ, φ, ψ, ϵ_A)
  return rot1(-ϵ_A) * rot3(-ψ) * rot1(φ) * rot3(γ)
end

# P MOD => GCRF
function precession_2006_EQU_FW(tt_jc)
  ϵ_A = obliquity_IAU_2000(tt_jc)

  γ, φ, ψ = precession_2006_EQU_FW_ang(tt_jc)
  return rot1(-ϵ_A) * rot3(-ψ) * rot1(φ) * rot3(γ)
end


# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.6.1 - Equation 5.33
# B GCRF => ITRF
function frame_bias_IAU_2006()
  δα_0 = as_to_rad(-0.014_6)
  ξ_0 = as_to_rad(-0.016_617_0)
  η_0 = as_to_rad(-0.006_819_2)
  return rot1(-η_0) * rot2(ξ_0) * rot3(δα_0)
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.5.7 - Equation 5.30
function GAST_IAU_2000A_2006(tt_jc, ut1_jd, Δψ, ϵ_A)
  ERA = earth_rotation_angle(ut1_jd)
  EO = equation_of_the_origins_IAU_2000A_2006(tt_jc, Δψ, ϵ_A)
  return ERA - EO
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.5.7 - Equation 5.31 - Table 5.2e
function equation_of_the_origins_IAU_2000A_2006(tt_jc, Δψ, ϵ_A)
  EO_poly_coeff = [0.0, 0.0]
  EO_idx = [1, 34, 35]

  fund_args = vcat(lunisolar_fundamental_arguments(tt_jc), planetary_fundamental_arguments(tt_jc))

  EO_p = @evalpoly(tt_jc, 0.014_506, 4612.156_534, 1.391_581_7, -0.000_000_44, -0.000_029_956, -0.000_000_036_8)
  EO_eq = Δψ * cos(ϵ_A)
  EO_np = sum_harmonic_components(tt_jc, nutation_2000A_GST_coeffs, EO_poly_coeff, EO_idx, fund_args)

  EO = -as_to_rad(EO_p) - EO_eq - EO_np
end

function precession_nutation_IAU_2000_2006_EQU(tt_jc, ut1_jd)
  ϵ_A = obliquity_IAU_2000(tt_jc)

  Δψ, Δϵ = nutation_IAU_2000A_R06_EQU_ang(tt_jc)
  N = nutation_IAU_2000A_R06_EQU(Δψ, Δϵ, ϵ_A)

  ψ_a, ω_a, χ_a = precession_2006_EQU_ang(tt_jc)
  ϵ_0 = as_to_rad(84_381.406)
  P = precession_2006_EQU(ψ_a, ω_a, χ_a, ϵ_0)

  B = frame_bias_IAU_2006()

  R = CIRS_to_TIRS_EQU(tt_jc, ut1_jd, Δψ, ϵ_A)

  return N, P, B, R
end

function precession_nutation_IAU_2000_2006_EQU_FW(tt_jc, ut1_jd)
  ϵ_A = obliquity_IAU_2000(tt_jc)

  Δψ, Δϵ = nutation_IAU_2000A_R06_EQU_ang(tt_jc)
  N = nutation_IAU_2000A_R06_EQU(Δψ, Δϵ, ϵ_A)

  γ, φ, ψ = precession_2006_EQU_FW_ang(tt_jc)
  P = precession_2006_EQU_FW(γ, φ, ψ, ϵ_A)

  R = CIRS_to_TIRS_EQU(tt_jc, ut1_jd, Δψ, ϵ_A)

  return N, P, eye(P), R
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.5.2 - Equation 5.13
function TIO_locator_IAU_2006(tt_jc)
  return as_to_rad(-0.000_047) * tt_jc
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.4.1 - Equation 5.3
# W TIRS => ITRF
function polar_motion(tt_jc, xp, yp)
  s′ = TIO_locator_IAU_2006(tt_jc)
  return rot2(-xp) * rot1(-yp) * rot3(s′)
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.4.2 - Equation 5.5
# R CIRS => TIRS
function CIRS_to_TIRS_CIO(ut1_jd)
  Θ = earth_rotation_angle(ut1_jd)
  return rot3(Θ)
end

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.4 / 5.5.7
# R CIRS => TIRS
function CIRS_to_TIRS_EQU(tt_jc, ut1_jd, Δψ, ϵ_A)
  GAST = GAST_IAU_2000A_2006(tt_jc, ut1_jd, Δψ, ϵ_A)
  return rot3(GAST)
end

## IAU-1976/FK5

function lunisolar_fundamental_arguments_1980(tdb_jc)
  # Mean Anomaly of the Moon
  M_m  = as_to_rad_mod( @evalpoly(tdb_jc, 134.962_981_51, (1325 * 360 + 198.867_398_1),  0.008_697_2,  1.78e-5) * 3600)

  # Mean Anomaly of the Sun
  M_s  = as_to_rad_mod( @evalpoly(tdb_jc, 357.527_723_33, (  99 * 360 + 359.050_340_0), -0.000_160_3, -3.3e-6) * 3600)

  # Mean Argument of Latitude of the Moon
  u_sm = as_to_rad_mod( @evalpoly(tdb_jc,  93.271_910_28, (1342 * 360 +  82.017_538_1), -0.003_682_5,  3.1e-6) * 3600)

  # Mean Elongation of the Moon from the Sun
  D_s  = as_to_rad_mod( @evalpoly(tdb_jc, 297.850_363_06, (1236 * 360 + 307.111_480_0), -0.001_914_2,  5.3e-6) * 3600)

  # Mean Longitude of the ascending Node of the Moon
  Ω_m  = as_to_rad_mod( @evalpoly(tdb_jc, 125.044_522_22,-(   5 * 360 + 134.136_260_8),  0.002_070_8,  2.2e-6) * 3600)

  return @SVector [M_m, M_s, u_sm, D_s, Ω_m]
end

function precession_IAU_1980_ang(tt_jc)
  ζ = @evalpoly(tt_jc, 0, 2306.2181,  0.301_88,  0.017_998)
  Θ = @evalpoly(tt_jc, 0, 2004.3109, -0.426_65, -0.041_833)
  z = @evalpoly(tt_jc, 0, 2306.2181,  1.094_68,  0.018_203)
  return as_to_rad(ζ), as_to_rad(Θ), as_to_rad(z)
end

# P GCRF => ITRF
function precession_IAU_1980(ζ, Θ, z)
  return rot3(-z) * rot2(Θ) * rot3(-ζ)
end

function nutation_IAU_1980_EQU_ang(tt_jc)
  fund_args = lunisolar_fundamental_arguments_1980(tt_jc)

  Δψ = 0.0
  Δϵ = 0.0

  for i in 1:106
    a_i = dot(nutation_coeffs_arguments[i, :], fund_args)

    Δψ += (nutation_coeffs_longitude[i, 1] + nutation_coeffs_longitude[i, 2] * tt_jc) * sin(a_i)
    Δϵ += (nutation_coeffs_obliquity[i, 1] + nutation_coeffs_obliquity[i, 2] * tt_jc) * cos(a_i)
  end

  return μas_to_rad(Δψ * 100), μas_to_rad(Δϵ * 100)
end

function nutation_IAU_1980_EQU(Δψ, Δϵ, ϵ_A)
  return rot1(- ϵ_A - Δϵ) * rot3(-Δψ) * rot1(ϵ_A)
end

function obliquity_IAU_1980(tt_jc)
  return as_to_rad(@evalpoly(tt_jc, 84_381.448, -46.8150, -0.000_59, 0.001_813))
end

# W GCRF => ITRF
function polar_motion_IAU_1980(xp, yp)
  return rot2(-xp) * rot1(-yp)
end

function equation_of_the_equinoxes_IAU_1982(Δψ, ϵ_A, tt_jc; post1997 = true)
  if (post1997)
    Ω_m = as_to_rad_mod( @evalpoly(tt_jc, 125.044_522_22,-(   5 * 360 + 134.136_260_8),  0.002_070_8,  2.2e-6) * 3600)
    return Δψ * cos(ϵ_A) + as_to_rad(0.002_64 * sin(Ω_m) + 0.000_063 * sin(2 * Ω_m))
  else
    return Δψ * cos(ϵ_A)
  end
end

function GMST_IAU_1980(ut1_jc)
  return mod(@evalpoly(ut1_jc, 67_310.548_41, (876_600 * 3600 + 8_640_184.812_866), 0.093_104, -6.2e-6), 86400) / 86400 * 2π
end

function GAST_IAU_1982(ut1_jc, Δψ, ϵ_A, tt_jc)
  return GMST_IAU_1980(ut1_jc) + equation_of_the_equinoxes_IAU_1982(Δψ, ϵ_A, tt_jc)
end

function TOD_to_PEF_1980(tt_jc, ut1_jc, Δψ, ϵ_A)
  GAST_1982 = GAST_IAU_1982(ut1_jc, Δψ, ϵ_A, tt_jc)
  return rot3(GAST_1982)
end

function precession_nutation_IAU_1980(tt_jc, ut1_jc; gcrf_correction_Δψ = 0.0, gcrf_correction_Δϵ = 0.0, LOD = 0.0)
  ϵ_A = obliquity_IAU_1980(tt_jc)

  Δψ, Δϵ = nutation_IAU_1980_EQU_ang(tt_jc)
  Δψ += gcrf_correction_Δψ
  Δϵ += gcrf_correction_Δϵ
  N = nutation_IAU_1980_EQU(Δψ, Δϵ, ϵ_A)

  ζ, Θ, z = precession_IAU_1980_ang(tt_jc)
  P = precession_IAU_1980(ζ, Θ, z)

  θ_GAST = GAST_IAU_1982(ut1_jc, Δψ, ϵ_A, tt_jc)
  R = rot3(θ_GAST)

  ωE = 7.292_115_146_706_979e-5 * (1 - LOD/86400)

  Ṙ = @SMatrix [
    (-ωE * sin(θ_GAST)) (-ωE * cos(θ_GAST)) 0 ;
    ( ωE * cos(θ_GAST)) (-ωE * sin(θ_GAST)) 0 ;
    0                   0                   0 ]

  return N, P, R, Ṙ.'
end

## TEME

function TOD_to_TEME(Δψ, ϵ_A, tt_jc)
  eq_1982 = equation_of_the_equinoxes_IAU_1982(Δψ, ϵ_A, tt_jc; post1997 = true) # post 1997 corrections ???
  return rot3(eq_1982)
end

function PEF_to_TEME(ut1_jc)
  θ_GMST = GMST_IAU_1980(ut1_jc)
  return rot3(-θ_GMST)
end

## FK4

function FK4_to_J2000()
  return @SMatrix [ 0.999_925_679_495_687_7 -0.011_181_483_220_466_2 -0.004_859_003_815_359_2 ;
                    0.011_181_483_239_171_7  0.999_937_484_893_313_5 -0.000_027_162_594_714_2 ;
                    0.004_859_003_772_314_3 -0.000_027_170_293_744_0  0.999_988_194_602_374_2]
end

## IERS 2010 Utils

function polar_motion_libration(ut1_jd, tt_jd)
  tt_jc = julianYears(tt_jd)

  gmst = GMST_IAU_2006(ut1_jd, tt_jc)

  fundamental_args = SVector{6,Float64}(mod2pi(gmst + π), lunisolar_fundamental_arguments(tt_jc)...)

  ΔX = 0.0
  ΔY = 0.0

  for i in size(UT1_LOD_libration_coeffs, 1):-1:1
    a_i = dot(polar_motion_libration_coeffs[i, 1:6], fundamental_args)
    sin_a_i = sin(a_i)
    cos_a_i = cos(a_i)

    ΔX += polar_motion_libration_coeffs[i,  8] * sin_a_i + polar_motion_libration_coeffs[i,  9] * cos_a_i
    ΔY += polar_motion_libration_coeffs[i, 10] * sin_a_i + polar_motion_libration_coeffs[i, 11] * cos_a_i
  end

  ΔX += polar_motion_libration_secular_coeffs[1] * tt_jc * 100
  ΔY += polar_motion_libration_secular_coeffs[2] * tt_jc * 100

  return ΔX, ΔY
end

function UT1_LOD_libration(ut1_jd, tt_jd)
  tt_jc = julianYears(tt_jd)

  gmst = GMST(ut1_jd, tt_jc)

  fundamental_args = SVector{6,Float64}(mod2pi(gmst + π), lunisolar_fundamental_arguments(tt_jc)...)

  ΔUT1 = 0.0
  ΔLOD = 0.0

  for i in size(UT1_LOD_libration_coeffs, 1):-1:1
    a_i = dot(UT1_LOD_libration_coeffs[i, 1:6], fundamental_args)
    sin_a_i = sin(a_i)
    cos_a_i = cos(a_i)

    ΔUT1 += UT1_LOD_libration_coeffs[i,  8] * sin_a_i + UT1_LOD_libration_coeffs[i,  9] * cos_a_i
    ΔLOD += UT1_LOD_libration_coeffs[i, 10] * sin_a_i + UT1_LOD_libration_coeffs[i, 11] * cos_a_i
  end

  return ΔUT1, ΔLOD
end

function polar_motion_ocean_tides(tt_jc)
  return ΔX, ΔY
end

function UT1_LOD_ocean_tides(tt_jc)
  return ΔUT1, ΔLOD
end
