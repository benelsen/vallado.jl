include("state.jl")
include("nutation.jl")

function hohmann_circular(r_initial, r_final; μ = 398600.4418)
  a_trans = (r_initial + r_final) / 2

  v_initial = sqrt(μ / r_initial)
  v_final = sqrt(μ / r_final)

  v_trans_a = sqrt(2μ / r_initial - μ / a_trans)
  v_trans_b = sqrt(2μ / r_final - μ / a_trans)

  Δv_a = v_trans_a - v_initial
  Δv_b = v_final - v_trans_b
  Δv = abs(Δv_a) + abs(Δv_b)

  τ = π * sqrt(a_trans^3 / μ)

  return Δv, [Δv_a, Δv_b], τ
end

function bielliptic_circular(r_initial, r_intermediate, r_final; μ = 398600.4418)
  a_trans1 = (r_initial + r_intermediate) / 2
  a_trans2 = (r_intermediate + r_final) / 2

  v_initial = sqrt(μ / r_initial)
  v_final = sqrt(μ / r_final)

  v_trans1_a = sqrt(2μ / r_initial - μ / a_trans1)
  v_trans1_b = sqrt(2μ / r_intermediate - μ / a_trans1)
  v_trans2_b = sqrt(2μ / r_intermediate - μ / a_trans2)
  v_trans2_c = sqrt(2μ / r_final - μ / a_trans2)

  Δv_a = v_trans1_a - v_initial
  Δv_b = v_trans2_b - v_trans1_b
  Δv_c = v_final - v_trans2_c
  Δv = abs(Δv_a) + abs(Δv_b) + abs(Δv_c)

  τ = π * (sqrt(a_trans1^3 / μ) + sqrt(a_trans2^3 / μ))

  return Δv, [Δv_a, Δv_b, Δv_c], τ
end

function one_tangent_circular(r_initial, r_final, ν_trans_b, E_0 = 0.0, k = 0; μ = 398600.4418, periapsis = true)
  R_inv = r_initial / r_final

  if periapsis
    e_trans = (R_inv - 1) / (cos(ν_trans_b) - R_inv)
    a_trans = r_initial / (1 - e_trans)
  else
    e_trans = (R_inv - 1) / (cos(ν_trans_b) + R_inv)
    a_trans = r_initial / (1 + e_trans)
  end

  v_initial = sqrt(μ / r_initial)
  v_final = sqrt(μ / r_final)

  v_trans_a = sqrt(2μ / r_initial - μ / a_trans)
  v_trans_b = sqrt(2μ / r_final - μ / a_trans)

  ϕ_trans_b = atan2(e_trans * sin(ν_trans_b), (1 + e_trans * cos(ν_trans_b)))

  Δv_a = v_trans_a - v_initial
  Δv_b = sqrt(v_trans_b^2 + v_final^2 - 2 * v_trans_b * v_final * cos(ϕ_trans_b))
  Δv = abs(Δv_a) + abs(Δv_b)

  E = ν_to_EBH(ν_trans_b, e_trans)

  τ = sqrt(a_trans^3 / μ) * (k * 2π + (E - e_trans * sin(E)) - (E_0 - e_trans * sin(E_0)))
  return Δv, [Δv_a, Δv_b], τ
end

function launch_window(ut1_jd, λ, Ω, β, i; ω_E = 7.292_115_9e-5)
  ut1_jc = julianYears(floor(ut1_jd + 0.5) - 0.5)
  GMST0 = GMST_IAU_1980(ut1_jc)
  println(GMST0)

  λ_u = acos(cos(β) / sin(i))
  GMST = (Ω + λ_u - λ)
  println(GMST)
  return (GMST - GMST0) / ω_E / 86400
end

function inclination_change(Δi, ϕ_fpa, v_initial)
  return 2 * v_initial * cos(ϕ_fpa) * sin(Δi / 2)
end

function ascending_node_change_circular(ΔΩ, i_initial, v_initial)
  θ = acos(cos(i_initial)^2 + sin(i_initial)^2 * cos(ΔΩ))
  Δv = 2 * v_initial * sin(θ/2)
  u_initial = acos(tan(i_initial) * (cos(ΔΩ) - cos(θ)) / sin(θ))
  u_final = acos(cos(i_initial) * sin(i_initial) * (1 - cos(ΔΩ)) / sin(θ))
  return Δv, u_initial, u_final
end

function ascending_node_inclination_change_circular(i_initial, i_final, ΔΩ, v_initial)
  θ = acos(cos(i_initial) * cos(i_final) + sin(i_initial) * sin(i_final) * cos(ΔΩ))
  Δv = 2 * v_initial * sin(θ/2)
  u_initial = acos((sin(i_final) * cos(ΔΩ) - cos(θ) * sin(i_initial)) / (sin(θ) * cos(i_initial)))
  u_final = acos((cos(i_initial) * sin(i_final) - sin(i_initial) * cos(i_final) * cos(ΔΩ)) / sin(θ))
  return Δv, u_initial, u_final
end

function combined_plane_change(Δi, v_initial, v_trans_a, v_trans_b, v_final)
  # Δi = i_final - i_initial

  Δv_a = v_trans_a - v_initial
  Δv_b = v_final - v_trans_b

  R = v_initial / v_final * v_trans_a / v_trans_b

  s = atan2(sin(Δi), (R + cos(Δi))) / Δi

  for iter in 1:50
    s_new = asin((Δv_a * v_final * v_trans_b * sin((1 - s) * Δi)) / (Δv_b * v_initial * v_trans_a)) / Δi
    Δs = s_new - s
    s = s_new
    if abs(Δs) < 1e-8
      break
    end
  end

  Δi_initial = s * Δi
  Δi_final = (1 - s) * Δi

  Δv_initial = sqrt(v_initial^2 + v_trans_a^2 - 2 * v_initial * v_trans_a * cos(Δi_initial))
  Δv_final = sqrt(v_final^2 + v_trans_b^2 - 2 * v_final * v_trans_b * cos(Δi_final))

  Δv = Δv_initial + Δv_final

  return Δi_initial, Δi_final, Δv, Δv_initial, Δv_final
end

function combined_fixed_dv_circular(v_initial, v_final, Δv; increase = true)
  γ = acos(- (v_initial^2 + Δv^2 - v_final^2) / (2 * v_initial * Δv))

  if (increase && γ ≥ 0.0) || (!increase && γ ≤ 0.0)
    γ *= -1
  end

  Δi = acos((v_initial^2 + v_final^2 - Δv^2) / (2 * v_initial * v_final))

  if !increase
    Δi *= -1
  end

  return γ, Δi
end

function coplanar_phasing_circular_same(a_tgt, θ, k_tgt, k_int; μ = 398600.4418)
  ω_tgt = sqrt(μ / a_tgt^3)

  τ_phase = (2π * k_tgt + θ) / ω_tgt
  a_phase = cbrt(μ * (τ_phase / (2π * k_int))^2)

  Δv = 2 * abs(sqrt(2μ / a_tgt - μ / a_phase) - sqrt(μ / a_tgt))
  return τ_phase, a_phase, Δv
end

function coplanar_phasing_circular_different(a_tgt, a_int, θ_i, k; μ = 398600.4418)
  ω_tgt = sqrt(μ / a_tgt^3)
  ω_int = sqrt(μ / a_int^3)

  a_trans = (a_int + a_tgt) / 2
  τ_trans = π * sqrt(a_trans^3 / μ)

  α_L = ω_tgt * τ_trans
  θ = π - α_L
  τ_wait = (θ - θ_i + 2π * k) / (ω_int - ω_tgt)

  Δv = abs(sqrt(2μ / a_int - μ / a_trans) - sqrt(μ / a_int)) + abs(sqrt(2μ / a_tgt - μ / a_trans) - sqrt(μ / a_tgt))
  return τ_wait, τ_trans, a_trans, Δv
end

function noncoplanar_phasing(a_tgt, a_int, θ_i, k_tgt, k_int, u_int, Ω_int, λ_true0, Δi; μ = 398600.4418)
  ω_tgt = sqrt(μ / a_tgt^3)
  ω_int = sqrt(μ / a_int^3)

  a_trans = (a_int + a_tgt) / 2
  τ_trans = π * sqrt(a_trans^3 / μ)

  α_L = ω_tgt * τ_trans

  Δθ_int = π - u_int
  if Δθ_int < 0
    Δθ_int += π
  end

  Δt_node = Δθ_int / ω_int
  λ_true_tgt_1 = λ_true0 + ω_tgt * Δt_node
  λ_true_int_1 = Ω_int + π

  θ_new = λ_true_int_1 - λ_true_tgt_1
  α_new = π + θ_new

  p_phase = (α_new - α_L + 2π * k_tgt) / ω_tgt
  a_phase = cbrt(μ * (p_phase / (k_int * 2π))^2)

  Δv_phase = abs(sqrt(2μ / a_int - μ / a_phase) - sqrt(μ / a_int))
  Δv_trans_1 = abs(sqrt(2μ / a_int - μ / a_trans) - sqrt(2μ / a_int - μ / a_phase))
  Δv_trans_2 = sqrt((2μ / a_tgt - μ / a_trans) + μ / a_tgt - 2 * sqrt(2μ / a_tgt - μ / a_trans) * sqrt(μ / a_tgt) * cos(Δi))
  Δv = Δv_phase + Δv_trans_1 + Δv_trans_2

  τ_phase = 2π * sqrt(a_phase^3 / μ)
  τ = Δt_node + τ_phase + τ_trans

  return τ, Δt_node, τ_phase, τ_trans, Δv, Δv_phase, Δv_trans_1, Δv_trans_2, a_phase
end
