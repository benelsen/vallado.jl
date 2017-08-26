
using StaticArrays
using SpecialFunctions

include("matrix.jl")

abstract type OrbitalState end
abstract type StateElements <: OrbitalState end

const Epoch = Union{DateTime, Real}

const μ = 398600.4418
const e_min = 1e-6
const i_min = 1e-6

struct StateVector <: OrbitalState
  r::SVector{3, Float64}
  v::SVector{3, Float64}
  t::Epoch
end

struct StateKepler <: StateElements
  p::Float64
  e::Float64
  i::Float64
  Ω::Float64
  ω::Float64
  ν::Float64
  t::Epoch
end

struct StateKeplerEquatorial <: StateElements
  p::Float64
  e::Float64
  i::Float64
  ω_true::Float64
  ν::Float64
  t::Epoch
end

struct StateKeplerCircular <: StateElements
  p::Float64
  e::Float64
  i::Float64
  Ω::Float64
  u::Float64
  t::Epoch
end

struct StateKeplerCircularEquatorial <: StateElements
  p::Float64
  e::Float64
  i::Float64
  λ_true::Float64
  t::Epoch
end

struct TwoLineElements <: StateElements
  n::Float64
  e::Float64
  i::Float64
  Ω::Float64
  ω::Float64
  M::Float64
  ṅ::Float64
  n̈::Float64
  B_star::Float64
  t::Epoch
end

struct StateEquinoctial <: StateElements
  a_f::Float64
  a_g::Float64
  L::Float64
  Χ::Float64
  Ψ::Float64
  f_r::Int64
  t::Epoch
end

struct StateFlightElements <: StateElements
  α::Float64
  δ::Float64
  r::Float64
  v::Float64
  φ_fpa::Float64
  β::Float64
  t::Epoch
end

struct StateGeographicElements <: StateElements
  φ_gc::Float64
  λ::Float64
  r::Float64
  v::Float64
  φ_fpa::Float64
  β::Float64
  t::Epoch
end

struct StateDelaunay <: StateElements
  M::Float64
  ω::Float64
  Ω::Float64
  L_d::Float64
  h::Float64
  H_d::Float64
  t::Epoch
end

struct StatePoincare <: StateElements
  λ_M::Float64
  g_p::Float64
  h_p::Float64
  L_d::Float64
  G_p::Float64
  H_p::Float64
  t::Epoch
end

import Base.convert

function convert(::Type{StateDelaunay}, el::StateKepler)
  @assert el.e < 1

  a = el.p / (1 - el.e^2)
  E = ν_to_EBH(ν, e)
  M = E - el.e * sin(E)

  L_d = sqrt(μ * a)
  h = sqrt(μ * el.p)
  H_d = sqrt(μ * a * (1 - el.e^2)) * cos(el.i)

  return StateDelaunay(M, el.ω, el.Ω, L_d, h, H_d, el.t)
end

function convert(::Type{StatePoincare}, el::StateKepler)
  @assert el.e < 1

  a = el.p / (1 - el.e^2)
  E = ν_to_EBH(ν, e)
  M = E - el.e * sin(E)

  λ_M = M + el.ω + el.Ω
  L_d = sqrt(μ * a)
  g_p = sqrt(2 * L_d * (1 - sqrt(1 - el.e^2))) * cos(el.ω + el.Ω)
  h_p = sqrt(2 * L_d * sqrt(1 - el.e^2) * (1 - cos(el.i))) * cos(el.Ω)
  G_p = sqrt(2 * L_d * (1 - sqrt(1 - el.e^2))) * sin(el.ω + el.Ω)
  H_p = sqrt(2 * L_d * sqrt(1 - el.e^2) * (1 - cos(el.i))) * sin(el.Ω)

  return StatePoincare(λ_M, g_p, h_p, L_d, G_p, H_p, el.t)
end

function convert(::Type{StateEquinoctial}, el::StateKepler)
  @assert el.e < 1

  a = el.p / (1 - el.e^2)
  E = ν_to_EBH(el.ν, el.e)
  M = E - el.e * sin(E)

  f_r = el.i > π / 2 ? -1 : 1

  a_f = el.e * cos(el.ω + f_r * el.Ω)
  a_g = el.e * sin(el.ω + f_r * el.Ω)
  L = M + el.ω + f_r * el.Ω
  Χ = tan(el.i / 2)^f_r * sin(el.Ω)
  Ψ = tan(el.i / 2)^f_r * cos(el.Ω)

  return StateEquinoctial(a_f, a_g, L, Χ, Ψ, f_r, el.t)
end

function convert(::Type{StateFlightElements}, el::StateKepler)
  rv = elementsToState(el)

  r_mag = norm(rv.r)
  v_mag = norm(rv.v)

  δ = atan2(rv.r[3] / r_mag, hypot(rv.r[1], rv.r[2]) / r_mag)
  α = atan2(rv.r[2] / hypot(rv.r[1], rv.r[2]), rv.r[1] / hypot(rv.r[1], rv.r[2]))

  v_sez = rot2(π/2 - δ) * rot3(α) * rv.v

  # φ_fpa = atan2(hypot(v_sez[1], v_sez[2]) / v_mag, v_sez[3] / v_mag)
  φ_fpa = atan2(el.e * sin(el.ν), 1 + el.e * cos(el.ν))
  β = atan2(v_sez[2] / hypot(v_sez[1], v_sez[2]), -v_sez[1] / hypot(v_sez[1], v_sez[2]))

  return StateFlightElements(α, δ, r_mag, v_mag, φ_fpa, β, el.t)
end

function show_elements(io::IO, x::StateElements)
  t = typeof(x)
  println(io, "$t")
  for field in fieldnames(x)
    name = convert(String, field)
    value = getfield(x, field)
    if field in [:i, :Ω, :ω, :ν, :ω_true, :λ_true, :M, :φ_fpa, :β, :α, :δ, :φ_gc, :λ]
      value = rad2deg(value)
    end
    println(io, "$name: $value")
  end
end
Base.show(io::IO, x::StateElements) = show_elements(io, x)

function stateToElements(state::StateVector)
  t = state.t

  h = cross(state.r, state.v)
  n = cross(base_k, h)
  h_mag = norm(h)
  n_mag = norm(n)
  v_mag = norm(state.v)
  r_mag = norm(state.r)

  r_dot_v = dot(state.r, state.v)

  e = ((v_mag^2 - μ / r_mag) * state.r - r_dot_v * state.v) / μ
  e_mag = norm(e)

  ξ = v_mag^2 / 2 - μ / r_mag

  if abs(e_mag - 1) > 1e-12
    p = -μ / 2ξ * (1 - e_mag^2)
  else
    p = h_mag^2 / μ
  end

  i = acos(h[3] / h_mag)

  if e_mag > e_min

    if i > i_min # Normal Kelperian Elements

      Ω = acos(n[1] / n_mag)
      if n[2] < 0
        Ω = 2π - Ω
      end

      ω = acos(dot(n, e) / (n_mag * e_mag))
      if e[3] < 0
        2π - ω
      end

      ν = acos(dot(e, state.r) / (e_mag * r_mag))
      if r_dot_v < 0
        ν = 2π - ν
      end

      return StateKepler(p, e_mag, i, Ω, ω, ν, t)

    else # Elliptical Equatorial

      ω_true = acos(e[1] / e_mag)
      if e[2] < 0
        ω_true = 2π - ω_true
      end

      ν = acos(dot(e, state.r) / (e_mag * r_mag))
      if r_dot_v < 0
        ν = 2π - ν
      end

      return StateKeplerEquatorial(p, e_mag, i, ω_true, ν, t)
    end

  else

    if i > i_min # Circular Inclined

      Ω = acos(n[1] / n_mag)
      if n[2] < 0
        Ω = 2π - Ω
      end

      u = acos(dot(n, state.r) / (n_mag * r_mag))
      if state.r[3] < 0
        u = 2π - u
      end

      return StateKeplerCircular(p, e_mag, i, Ω, u, t)

    else # Circular Equatorial

      λ_true = acos(state.r[1] / r_mag)
      if state.r[2] < 0
        λ_true = 2π - λ_true
      end

      return StateKeplerCircularEquatorial(p, e_mag, i, λ_true, t)
    end
  end
end

function elementsToState(elements::StateElements; μ = 398600.4418)
  p = elements.p
  e = elements.e
  i = elements.i
  t = elements.t

  if isa(elements, StateKepler)
    Ω = elements.Ω
    ω = elements.ω
    ν = elements.ν
  elseif isa(elements, StateKeplerCircularEquatorial)
    Ω = 0.0
    ω = 0.0
    ν = elements.λ_true
  elseif isa(elements, StateKeplerCircular)
    Ω = elements.Ω
    ω = 0.0
    ν = elements.u
  elseif isa(elements, StateKeplerEquatorial)
    Ω = 0.0
    ω = elements.ω_true
    ν = elements.ν
  end

  r_PQW = @SVector [
    p * cos(ν) / (1 + e * cos(ν)),
    p * sin(ν) / (1 + e * cos(ν)),
    0
  ]

  v_PQW = @SVector [
    -sqrt(μ / p) * sin(ν),
    sqrt(μ / p) * (e + cos(ν)),
    0
  ]

  mat_PQW_to_IJK = rot3(-Ω) * rot1(-i) * rot3(-ω)

  r_IJK = mat_PQW_to_IJK * r_PQW
  v_IJK = mat_PQW_to_IJK * v_PQW

  return StateVector(r_IJK, v_IJK, t)
end

function kepler_equation_ell(M, e, max_iter = 50, ϵ = 1e-12)
  @assert 0.0 <= e < 1.0
  E = M + e * sin(M) + e^2 / 2.0 * sin(2M)
  for i = 1:max_iter
    ΔE = (M - E + e * sin(E)) / (1.0 - e * cos(E))
    E += ΔE
    if abs(ΔE) < ϵ
      return E
    end
  end
  return E
end

function kepler_equation_ell_series(M, f, max_ord = 50, ϵ = 1e-12)
  @assert 0.0 <= e < 1.0
  E = M
  for k in 1:max_ord
    ΔE = (2/k * besselj(k, k*f) * sin(k * M))
    E += ΔE
    if abs(ΔE) < ϵ
      return E
    end
  end
  return E
end

function kepler_equation_par(Δt, p)
  n_p = 2 * sqrt(μ / p^3)
  s = 0.5 * acot(3/2 * n_p * Δt)
  w = atan(cbrt(tan(s)))
  B = 2 * cot(2 * w)
end

function kepler_equation_par_cub(Δt, p)
  b = -3 * sqrt(μ / p^3) * Δt
  Δ = sqrt(1 + b^2)
  B = cbrt(-b + Δ) + cbrt(-b - Δ)
end

function kepler_equation_hyp(M, e, max_iter = 20, ϵ = 1e-12)
  @assert e > 1.0

  if e < 1.6
    if -π < M < 0.0 || M > π
      H = M - e
    else
      H = M + e
    end
  elseif e < 3.6 && abs(M) > π
    H = M - sign(M) * e
  else
    H = M / (e - 1)
  end

  for i = 1:max_iter
    ΔH = (M + H - e * sinh(H)) / (e * cosh(H) - 1)
    H += ΔH
    if abs(ΔH) < ϵ
      return H
    end
  end
  return H
end

function ν_to_EBH(ν, e)
  denom = 1 + e * cos(ν)
  if e < 1.0
    return E = atan2(sin(ν) * sqrt(1 - e^2) / denom, (e + cos(ν)) / denom)
  elseif e > 1.0
    ν_max = π - acos(1/e)
    @assert -ν_max <= ν <= ν_max
    return H = asinh(sin(ν) * sqrt(e^2 - 1) / (1 + e * cos(ν)))
  else
    return B = tan(ν/2)
  end
end

function EBH_to_ν(E, e)
  if e < 1.0
    denom = 1 - e * cos(E)
    return ν = atan2(sin(E) * sqrt(1 - e^2) / denom, (cos(E) - e) / denom)
  elseif e > 1.0
    denom = 1 - e * cosh(E)
    return ν = atan2(-sinh(E) * sqrt(e^2 - 1) / denom, (cosh(E) - e) / denom)
  else
    return 2 * atan(E)
  end
end

function ψ_to_c2c3(ψ)
  if ψ > 1e-6
    return (1.0 - cos(sqrt(ψ))) / ψ, (sqrt(ψ) - sin(sqrt(ψ))) / sqrt(ψ)^3
  elseif ψ < -1e-6
    return (1.0 - cosh(sqrt(-ψ))) / ψ, (sinh(sqrt(-ψ)) - sqrt(-ψ)) / sqrt(-ψ)^3
  else
    return 1/2, 1/6
  end
end

function ψ_to_c2c3_series(ψ, k_max = 10)
  ψk = 1.0
  c_2 = ψk / 2
  c_3 = ψk / 6
  g = 6.0
  for k in 1:k_max
    ψk *= (-ψ)
    g *= 2k + 2
    c_2 += ψk / g
    g *= 2k + 3
    c_3 += ψk / g
  end
  return c_2, c_3
end

propagate_kepler(state, epoch::Dates.TimeType) = propagate_kepler(state, Dates.value(epoch) / 1000.0)

function propagate_kepler(state, epoch::Real; ϵ = 1e-12, check_tol = 1e-12, max_iter = 20)
  Δt = (epoch - state.t)

  if abs(Δt) ≈ 0.
    return state
  end

  v_mag_2 = norm(state.v)^2
  r_mag = norm(state.r)

  r_dot_v = dot(state.r, state.v)

  ξ = v_mag_2 / 2 - μ / r_mag
  α = - v_mag_2 / μ + 2 / r_mag

  sqrt_μ = sqrt(μ)
  local Χ

  if α > 1e-6 # Elliptic
    Χ = sqrt_μ * Δt * α
    # Very close to solution, might fail to converge. Fix:
    if abs(α - 1) < 1e-6
      Χ += 0.05
    end

  elseif α < -1e-6 # Hyperpolic
    sqrt_a = sqrt(-1/α)
    Χ = sign(Δt) * sqrt_a * log( (-2μ * α * Δt) / (r_dot_v + sign(Δt) * sqrt_a * sqrt_μ * (1 - r_mag * α) ) )

  else # Parabolic
    h = norm(cross(state.r, state.v))
    p = h^2 / μ
    B = kepler_equation_par_cub(Δt, p)
    Χ = sqrt(p) * B
  end

  local ψ
  local r
  local c_2
  local c_3

  for i in 1:max_iter
    ψ = Χ^2 * α
    c_2, c_3 = ψ_to_c2c3(ψ)
    Χ_c_3 = Χ * (1 - ψ * c_3)
    r = Χ^2 * c_2 + r_dot_v / sqrt_μ * Χ_c_3 + r_mag * (1 - ψ * c_2)
    ΔΧ = (sqrt_μ * Δt - Χ^3 * c_3 - r_dot_v / sqrt_μ * Χ^2 * c_2 - r_mag * Χ_c_3) / r
    Χ += ΔΧ

    if abs(ΔΧ) < ϵ
      break
    end

    if iter === max_iter
      error("Did not converge after $iter iterations")
    end
  end

  Χ2_c_2 = Χ^2 * c_2
  f = 1.0 - Χ2_c_2 / r_mag
  g = Δt - Χ^3 / sqrt_μ * c_3
  ḟ = sqrt_μ / (r * r_mag) * Χ * (ψ * c_3 - 1.0)
  ġ = 1.0 - Χ2_c_2 / r

  check = f * ġ - ḟ * g - 1

  if !isapprox(check, 0; atol = check_tol)
    warn("check: ", check)
  end

  return StateVector(f * state.r + g * state.v, ḟ * state.r + ġ * state.v, epoch)
end

function find_tof(r0, r1, p)
  r0_mag = norm(r0)
  r1_mag = norm(r1)

  cosΔν = dot(r0, r1) / (r0_mag * r1_mag)
  Δν = acos(cosΔν)

  println(dot(r0, r1))

  k = r0_mag * r1_mag * (1 - cosΔν)
  l = r0_mag + r1_mag
  m = r0_mag * r1_mag * (1 + cosΔν)

  a = (m * k * p) / ((2m - l^2) * p^2 + 2 * k * l * p - k^2)

  f = 1 - r1_mag / p * (1 - cosΔν)
  g = r0_mag * r1_mag * sin(Δν) / sqrt(μ * p)

  if a > 1e-6
    ḟ = sqrt(μ / p) * tan(Δν / 2) * ((1 - cosΔν) / p - 1/r0_mag - 1/r1_mag)

    cosΔE = 1 - r0_mag / a * (1 - f)
    sinΔE = -r0_mag * r1_mag * ḟ / sqrt(a * μ)
    ΔE = atan2(sinΔE, cosΔE)

    return g + sqrt(a^3 / μ) * (ΔE - sinΔE)

  elseif a < -1e-6
    ΔH = acosh(1 + (f - 1) * r0_mag / a)
    return g + sqrt((-a)^3 / μ) * (sinh(ΔH) - ΔH)

  else
    c = sqrt(r0_mag^2 + r1_mag^2 - 2 * r0_mag * r1_mag * cosΔν)
    s = (r0_mag + r1_mag + c) / 2

    return 2/3 * sqrt(s^3 / 2μ) * (1 - ((s - c) / s)^(3/2))
  end
end

# Test

function test_find_tof()
  t0 = now()
  t1 = t0 + Dates.Millisecond(10 * 60 * 1e3)

  p = 7_500.000
  el0 = StateKepler(p, 0.1, deg2rad(15.0), deg2rad(0.0), deg2rad(0.0), deg2rad(0.0), t0)
  rv0 = elementsToState(el0)

  rv1 = propagate_kepler(rv0, t1)
  el1 = stateToElements(rv1)
  rv1b = elementsToState(el1)

  find_tof(rv0.r, rv1.r, p)
  find_tof(rv0.r, rv1.r, p*0.9)
  find_tof(rv0.r, rv1.r, p*1.1)

  println(el0)
  println(el1)

  println(rv1)
  println(rv1b)
end

function test_elementsToState()
  rv = StateVector([6524.834, 6862.875, 6448.296], [4.901327, 5.533756, -1.976341], now())
  el = stateToElements(rv)
  print(el)

  elements = StateKepler(200_000.000, 2.0, deg2rad(0.0), deg2rad(0.0), deg2rad(0.0), deg2rad(119.937), now())
  rv = elementsToState(elements)
  println(rv)

  el = stateToElements(rv)
  print(el)

  a = el.p / (1 - el.e^2)
  print(a *(1 - el.e), ", ", a *(1 + el.e))

  ϕ = 2 * asin(1 / el.e)
  println(rad2deg(ϕ))

  νc_d = rad2deg(acos(1 / el.e))
  println(-180 + νc_d, " < ", rad2deg(el.ν), " < ", 180 - νc_d)
end

function test_prop_kep()
  t0 = 0.
  t1 = 100 * 60

  el0 = StateKepler(14_000_000.000, 1.0, deg2rad(80.0), deg2rad(45.0), deg2rad(30.0), deg2rad(-20.0), t0)
  rv0 = elementsToState(el0)

  rv1 = propagate_kepler(rv0, t1)
  el1 = stateToElements(rv1)
  rv1b = elementsToState(el1)

  println(el0)
  println(el1)

  println(rv1)
  println(rv1b)
end

function test_ν_to_EBH()
  for e in linspace(0.0, 21.0, 1000)
    if e >= 1.0
      ν_max = π - acos(1/e) - 2 * eps(Float64)
    else
      ν_max = 2π
    end

    for ν in linspace(-ν_max, ν_max, 1000)
      E = ν_to_EBH(ν, e)
      err_ν = rem2pi(EBH_to_ν(E, e) - ν, RoundNearest)
      if abs(err_ν) > 1e-12
        println("e: $e, ν: $ν, E: $E, ν_max: $ν_max, error: $err_ν")
      end
    end
  end
end
