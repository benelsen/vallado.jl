using StaticArrays

include("utils.jl")
include("matrix.jl")
include("time.jl")

abstract type SphericalCoordinate <: FieldVector{3, Float64} end

struct GeocentricCoordinate <: SphericalCoordinate
  λ::Float64
  φ::Float64
  r::Float64
end

struct ReducedCoordinate <: SphericalCoordinate
  λ::Float64
  φ::Float64
  r::Float64
end

struct GeodeticCoordinate <: SphericalCoordinate
  λ::Float64
  φ::Float64
  h::Float64
end

struct CartesianCoordinate <: FieldVector{3, Float64}
  x::Float64
  y::Float64
  z::Float64
end

function geocentric_to_cartesian(sc::GeocentricCoordinate)
  return CartesianCoordinate(sc.r * [
    cos(sc.φ) * cos(sc.λ)
    cos(sc.φ) * sin(sc.λ)
    sin(sc.φ)
  ])
end

function cartesian_to_geocentric(cc::CartesianCoordinate)
  return GeocentricCoordinate(
    atan2(cc.y, cc.x),
    atan2(cc.z, hypot(cc.x, cc.y)),
    norm(cc)
  )
end

geodetic_to_cartesian(λ, ϕ, h) = geodetic_to_cartesian(GeodeticCoordinate(λ, ϕ, h))

function geodetic_to_cartesian(sc::GeodeticCoordinate)
  R = 6378.137
  e = 0.081_819_191_0435

  dn = sqrt(1 - (e * sin(sc.φ))^2)
  C = R / dn
  S = R * (1 - e^2) / dn
  r_δ = (C + sc.h) * cos(sc.φ)
  r_K = (S + sc.h) * sin(sc.φ)

  return CartesianCoordinate(r_δ * cos(sc.λ), r_δ * sin(sc.λ), r_K)
end

function cartesian_to_geodetic(cc::CartesianCoordinate)
  throw(NotImplementedException())
end

function reduced_to_cartesian(sc::ReducedCoordinate)
  throw(NotImplementedException())
end

function cartesian_to_reduced(cc::CartesianCoordinate)
  throw(NotImplementedException())
end

function SEZ_to_ECEF(sc::GeodeticCoordinate, ρ, ρ̇)
  r_site = geodetic_to_cartesian(sc)

  R_sez_to_ecef = rot3(-sc.λ) * rot2(-(π/2 - sc.φ))

  return R_sez_to_ecef * ρ + r_site, R_sez_to_ecef * ρ̇
end

function ECEF_to_SEZ(sc::GeodeticCoordinate, r_site, r, v)
  R_ecef_to_sez = transpose(rot3(-sc.λ) * rot2(-(π/2 - sc.φ)))

  return R_ecef_to_sez * (r - r_site), R_ecef_to_sez * v
end

function ECEF_to_SEZ(sc::GeodeticCoordinate, r, v)
  r_site = geodetic_to_cartesian(sc)

  R_ecef_to_sez = transpose(rot3(-sc.λ) * rot2(-(π/2 - sc.φ)))

  return R_ecef_to_sez * (r - r_site), R_ecef_to_sez * v
end

function IJKt_to_IJK(r_site, v_site, r, v)
  return r + r_site, v + v_site
end

function IJK_to_IJKt(r_site, v_site, r, v)
  return r - r_site, v - v_site
end

# Needs velocity conversions
function PQW_to_IJK(el, r, v)
  R_PQW_to_IJK = rot3(-el.Ω) * rot1(-el.i) * rot3(el.ω)
  return R_PQW_to_IJK * r, R_PQW_to_IJK * v
end

function IJK_to_PQW(el, r)
  R_IJK_to_PQW = transpose(rot3(-el.Ω) * rot1(-el.i) * rot3(el.ω))
  return R_IJK_to_PQW * r
end


function PQW_to_RSW(el, r)
  R_PQW_to_RSW = rot3(el.ν)
  return R_PQW_to_RSW * r
end

function RSW_to_PQW(el, r)
  R_RSW_to_PQW = rot3(-el.ν)
  return R_RSW_to_PQW * r
end


function EQW_to_ECI(el, r)
  f_r = el.i > 90 ? -1 : +1
  R_EQW_to_ECI = rot3(-el.Ω) * rot1(-el.i) * rot3(f_r * el.Ω)
  return R_EQW_to_ECI * r
end

function ECI_to_EQW(el, r)
  f_r = el.i > 90 ? -1 : +1
  R_ECI_to_EQW = rot3(-f_r * el.Ω) * rot1(el.i) * rot3(el.Ω)
  return R_ECI_to_EQW * r
end


function NTW_to_ECI(el, r)
  R_NTW_to_ECI = rot3(-el.Ω) * rot1(-el.i) * rot3(-el.u) * rot3(el.φ_fpa)
  return R_NTW_to_ECI * r
end

function ECI_to_NTW(el, r)
  R_ECI_to_NTW = transpose(rot3(-el.Ω) * rot1(-el.i) * rot3(-el.u) * rot3(el.φ_fpa))
  return R_ECI_to_NTW * r
end


function RSW_to_ECI(el, r)
  R_RSW_to_ECI = rot3(-el.Ω) * rot1(-el.i) * rot3(-el.u)
  return R_RSW_to_ECI * r
end

function ECI_to_RSW(el, r)
  R_ECI_to_RSW = transpose(rot3(-el.Ω) * rot1(-el.i) * rot3(-el.u))
  return R_ECI_to_RSW * r
end


function SEZ_to_body(att, r)
  R_SEZ_to_body = rot3(att.yaw) * rot1(att.roll) * rot2(att.pitch)
  return R_SEZ_to_body * r
end

function body_to_SEZ(att, r)
  R_body_to_SEZ = transpose(rot3(att.yaw) * rot1(att.roll) * rot2(att.pitch))
  return R_body_to_SEZ * r
end

# Algorithm 12 - Vallado et al. 2013 page 172 (199)
# Astronomical Almanac (1992:K12)
function ECEF_to_geodetic(cc::CartesianCoordinate, max_iter = 20)
  R = 6378.137
  e = 0.081_819_191_0435

  r_δ = hypot(cc.x, cc.y)

  λ = α = atan2(cc.y, cc.x)
  φ = δ = atan(cc.z / r_δ)

  for i in 1:max_iter
    φ_old = φ
    C = R / sqrt(1 - (e * sin(φ))^2)
    φ = atan((cc.z + C * e^2 * sin(φ)) / r_δ)

    if abs(φ_old - φ) < 1e-12
      break
    end
  end

  if abs(φ) < deg2rad(89)
    C = R / sqrt(1 - (e * sin(φ))^2)
    h_ell = r_δ / cos(φ) - C
  else
    S = R * (1 - e^2) / sqrt(1 - (e * sin(φ))^2)
    h_ell = cc.z / sin(φ) - S
  end

  return GeodeticCoordinate(λ, φ, h_ell)
end

# ECEF_to_geodetic(CartesianCoordinate(6524.834, 6862.875, 6448.296))

# Algorithm 13 - Vallado et al. 2013 page 173 (200)
# Borokowski (1989)
function ECEF_to_geodetic_direct(cc::CartesianCoordinate)
  a = 6378.137
  e = 0.081_819_191_0435

  r_δ = hypot(cc.x, cc.y)
  b = a * sqrt(1 - e^2) * sign(cc.z)

  E = (b * cc.z - (a^2 - b^2)) / (a * r_δ)
  F = (b * cc.z + (a^2 - b^2)) / (a * r_δ)

  P = 4/3 * (E * F + 1)
  Q = 2 * (E^2 - F^2)
  D = P^3 + Q^2

  if D > 0
    ν = cbrt(sqrt(D) - Q) - cbrt(sqrt(D) + Q)
  else
    ν = 2 * sqrt(-P) * cos(1/3 * acos(Q / (P * sqrt(-P))))
  end

  G = 0.5 * (sqrt(E^2 + ν) + E)
  t = sqrt(G^2 + (F - ν * G) / (2G - E)) - G

  λ = α = atan2(cc.y, cc.x)
  φ = atan(a * (1 - t^2) / (2b * t))
  h_ell = (r_δ - a * t) * cos(φ) + (cc.z - b) * sin(φ)

  return GeodeticCoordinate(λ, φ, h_ell)
end

# ECEF_to_geodetic_direct(CartesianCoordinate(6524.834, 6862.875, 6448.296))

function radec_to_state(r, α, δ)
  return r .* SVector{3}(
    cos(δ) * cos(α),
    cos(δ) * sin(α),
    sin(δ)
  )
end

function radec_to_state(r, α, δ, ṙ, α̇, δ̇)
  x = radec_to_state(r, α, δ)
  v = SVector{3}(
    (ṙ * cos(δ) * cos(α) - r * sin(δ) * cos(α) * δ̇ - r * cos(δ) * sin(α) * α̇),
    (ṙ * cos(δ) * sin(α) - r * sin(δ) * sin(α) * δ̇ + r * cos(δ) * cos(α) * α̇),
    (ṙ * sin(δ) + r * cos(δ) * δ̇)
  )
  return x, v
end

function state_to_radec(x)
  r = norm(x)

  α = atan2(x[2], x[1])
  δ = asin(x[3] / r)

  return r, α, δ
end

function state_to_radec(x, v)
  r = norm(x)
  r_xy = hypot(x[1], x[2])

  if r_xy ≉ 0
    α = atan2(x[2], x[1])
  else
    α = atan2(v[2], v[1])
  end

  δ = asin(x[3] / r)

  ṙ = dot(x, v) / r
  α̇ = - (v[1] * x[2] - v[2] * x[1]) / (x[1]^2 + x[2]^2)
  δ̇ = (v[3] - ṙ * x[3] / r) / r_xy

  return r, α, δ, ṙ, α̇, δ̇
end

function azel_to_SEZ(ρ, β, e)
  return SVector{3}(
    -ρ * cos(e) * cos(β),
     ρ * cos(e) * sin(β),
     ρ * sin(e)
  )
end

function azel_to_SEZ(ρ, β, e, ρ̇, β̇, ė)
  x = azel_to_SEZ(ρ, β, e)
  v = SVector{3}(
    (-ρ̇ * cos(e) * cos(β) + ρ * sin(e) * cos(β) * ė + ρ * cos(e) * sin(β) * β̇),
    ( ρ̇ * cos(e) * sin(β) - ρ * sin(e) * sin(β) * ė + ρ * cos(e) * cos(β) * β̇),
    ( ρ̇ * sin(e) + ρ * cos(e) * ė)
  )
  return x, v
end

function SEZ_to_azel(x)
  ρ = norm(x)

  e = asin(x[3] / ρ)
  β = atan2(x[2], -x[1])

  return ρ, β, e
end

function SEZ_to_azel(x, v)
  ρ = norm(x)
  ρ_SE = hypot(x[1], x[2])

  e = asin(x[3] / ρ)

  if abs(e) ≉ π
    β = atan2(x[2], -x[1])
  else
    β = atan2(v[2], -v[1])
  end

  ρ̇ = dot(x, v) / ρ
  β̇ = (v[1] * x[2] - v[2] * x[1]) / ρ_SE^2
  ė = (v[3] - ρ̇ * sin(e)) / ρ_SE

  return ρ, β, e, ρ̇, β̇, ė
end

function geocentric_to_topocentric(x_site, v_site, x, v)
  return x - x_site, v - v_site
end

function geocentric_to_topocentric(x_site, x)
  return x - x_site
end

function topocentric_to_geocentric(x_site, v_site, ρ, ρ̇)
  return ρ + x_site, ρ̇ + v_site
end

function topocentric_to_geocentric(x_site, ρ)
  return ρ + x_site
end

function radec_to_azel(α_t, δ_t, ϕ, θ_LST)
  LHA = θ_LST - α_t

  e = asin(sin(ϕ) * sin(δ_t) + cos(ϕ) * cos(δ_t) * cos(LHA))
  β = atan2(-sin(LHA) * cos(δ_t) / cos(e), (sin(δ_t) - sin(e) * sin(ϕ)) / (cos(e) * cos(ϕ)))

  return β, e
end

function azel_to_radec(β, e, ϕ, θ_LST)
  δ_t = asin(sin(e) * sin(ϕ) + cos(e) * cos(ϕ) * cos(β))

  LHA = atan2(-sin(β) * cos(e) / cos(δ_t), (cos(ϕ) * sin(e) - sin(ϕ) * cos(β) * cos(e)) / cos(δ_t))

  α_t = θ_LST - LHA

  return α_t, δ_t
end

ϵ_0 = deg2rad(84_381.406 / 3600)

function ecliptic_to_radec(λ, ϕ, ϵ = ϵ_0)
  δ = asin(sin(ϕ) * cos(ϵ) + cos(ϕ) * sin(ϵ) * sin(λ))
  α = atan2((-sin(ϕ) * sin(ϵ) + cos(ϕ) * cos(ϵ) * sin(λ)) / cos(δ), (cos(ϕ) * cos(λ)) / cos(δ))
  return α, δ
end

function ecliptic_to_radec(λ, ϕ, λ̇, ϕ̇, ϵ = ϵ_0)
  δ = asin(sin(ϕ) * cos(ϵ) + cos(ϕ) * sin(ϵ) * sin(λ))
  α = atan2((-sin(ϕ) * sin(ϵ) + cos(ϕ) * cos(ϵ) * sin(λ)) / cos(δ), (cos(ϕ) * cos(λ)) / cos(δ))

  δ̇ = (cos(ϕ) * cos(ϵ) * ϕ̇ - sin(ϕ) * sin(ϵ) * sin(λ) * ϕ̇ + cos(ϕ) * sin(ϵ) * cos(λ) * λ̇) / cos(δ)

  if sin(λ) ≉ 0
    α̇ = (sin(ϕ) * cos(λ) * ϕ̇ + cos(ϕ) * sin(λ) * λ̇ - sin(δ) * cos(α) * δ̇) / (cos(δ) * sin(α))
  else
    α̇ = (-cos(ϕ) * sin(ϵ) * ϕ̇ - sin(ϕ) * cos(ϵ) * sin(λ) * ϕ̇ + cos(ϕ) * cos(ϵ) * cos(λ) * λ̇ + sin(δ) * sin(α) * δ̇) / (cos(δ) * cos(α))
  end

  return α, δ, α̇, δ̇
end

function radec_to_ecliptic(α, δ, ϵ = ϵ_0)
  ϕ = asin(-cos(δ) * sin(α) * sin(ϵ) + sin(δ) * cos(ϵ))
  λ = atan2((cos(δ) * sin(α) * cos(ϵ) + sin(δ) * sin(ϵ)) / cos(ϕ), (cos(δ) * cos(α)) / cos(ϕ))
  return λ, ϕ
end

function radec_to_ecliptic(α, δ, α̇, δ̇, ϵ = ϵ_0)
  ϕ = asin(-cos(δ) * sin(α) * sin(ϵ) + sin(δ) * cos(ϵ))
  λ = atan2((cos(δ) * sin(α) * cos(ϵ) + sin(δ) * sin(ϵ)) / cos(ϕ), (cos(δ) * cos(α)) / cos(ϕ))

  ϕ̇ = (sin(δ) * sin(α) * sin(ϵ) * δ̇ - cos(δ) * cos(α) * sin(ϵ) * α̇ + cos(δ) * cos(ϵ) * δ̇) / cos(ϕ)

  if sin(λ) ≉ 0
    λ̇ = - (sin(ϕ) * cos(λ) * ϕ̇ - sin(δ) * cos(α) * δ̇ - cos(δ) * sin(α) * α̇) / (cos(ϕ) * sin(λ))
  else
    λ̇ = (sin(ϕ) * sin(λ) * ϕ̇ - sin(δ) * sin(α) * cos(ϵ) * δ̇ + cos(δ) * cos(α) * cos(ϵ) * α̇ + cos(δ) * sin(ϵ) * δ̇) / (cos(ϕ) * cos(λ))
  end

  return λ, ϕ, λ̇, ϕ̇
end

"""
  ecliptic_to_cartesian(r, λ, ϕ)

Convert ecliptic spherical coordinates to ecliptic cartestian coordinates.
"""
function ecliptic_to_cartesian(r, λ, ϕ)
  return r .* SVector{3}(
    cos(ϕ) * cos(λ),
    cos(ϕ) * sin(λ),
    sin(ϕ)
  )
end
