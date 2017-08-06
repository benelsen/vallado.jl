include("../maneuvers.jl")

# 6.1
hohmann_circular(6569.4781, 42159.4856)

# 6.2
bielliptic_circular(191.34411 + 6378.137, 503873 + 6378.137, 376310 + 6378.137)

# 6.3
one_tangent_circular(191.34411 + 6378.137, 35781.34857 + 6378.137, deg2rad(160.0))

# Table 6.3
μ = 398600.4418
r_initial = 6671.53
v_initial = sqrt(μ / r_initial)

r_final = 42163.95

a_trans = (r_initial + r_final) / 2
v_final = sqrt(μ / r_final)

v_trans_a = sqrt(2μ / r_initial - μ / a_trans)
v_trans_b = sqrt(2μ / r_final - μ / a_trans)

Δi_initial, Δi_final, Δv, Δv_initial, Δv_final, s, iter = combined_plane_change(deg2rad(10.0), v_initial, v_trans_a, v_trans_b, v_final)
rad2deg(Δi_initial)
rad2deg(Δi_final)

Δi_initial, Δi_final, Δv, Δv_initial, Δv_final, s, iter = combined_plane_change(deg2rad(45.0), v_initial, v_trans_a, v_trans_b, v_final)
rad2deg(Δi_initial)
rad2deg(Δi_final)

r_final = 26558.56

a_trans = (r_initial + r_final) / 2
v_final = sqrt(μ / r_final)

v_trans_a = sqrt(2μ / r_initial - μ / a_trans)
v_trans_b = sqrt(2μ / r_final - μ / a_trans)

Δi_initial, Δi_final, Δv, Δv_initial, Δv_final, s, iter = combined_plane_change(deg2rad(28.5), v_initial, v_trans_a, v_trans_b, v_final)
rad2deg(Δi_initial)
rad2deg(Δi_final)

# 6.4
inclination_change(deg2rad(15.0), 0.0, 5.892311)

e = 0.3
p = 17858.7836
ω = deg2rad(30.0)
a = p / (1 - e^2)

ν = deg2rad(330.0)
r = p / (1 + e * cos(ν))
v = sqrt(2μ / r - μ / a)
ϕ_fpa = atan2(e * sin(ν), 1 + e * cos(ν))

ν = deg2rad(330.0 - 180)
r = p / (1 + e * cos(ν))
v = sqrt(2μ / r - μ / a)
ϕ_fpa = atan2(e * sin(ν), 1 + e * cos(ν))


inclination_change(deg2rad(15.0), ϕ_fpa, v)

# 6.5
Δv, u_i, u_f = ascending_node_change_circular(deg2rad(45.0), deg2rad(55.0), 5.892311)
rad2deg(u_i), rad2deg(u_f)

# 6.6
Δv, u_i, u_f = ascending_node_inclination_change_circular(deg2rad(55.0), deg2rad(40.0), deg2rad(45.0), 5.892311)
rad2deg(u_i), rad2deg(u_f)

# 6.7
r_initial = 6378.137 + 191
r_final = 6378.137 + 35780
a_trans = (r_initial + r_final) / 2
v_initial = sqrt(μ / r_initial)
v_final = sqrt(μ / r_final)
v_trans_a = sqrt(2μ / r_initial - μ / a_trans)
v_trans_b = sqrt(2μ / r_final - μ / a_trans)
Δi_initial, Δi_final, Δv, Δv_initial, Δv_final = combined_plane_change(deg2rad(-28.5), v_initial, v_trans_a, v_trans_b, v_final)
rad2deg(Δi_initial), rad2deg(Δi_final)

γ_a, = combined_fixed_dv_circular(v_initial, v_trans_a, Δv_initial; increase = false)
γ_b, = combined_fixed_dv_circular(v_trans_b, v_final, Δv_final; increase = false)
rad2deg(γ_a), rad2deg(γ_b)

# 6.8
μ = 398600.4418

a_tgt = a_int = 12756.274
θ = deg2rad(-20.0)

τ_phase, a_phase, Δv = coplanar_phasing_circular_same(a_tgt, θ, 1, 1)
τ_phase, a_phase, Δv = coplanar_phasing_circular_same(a_tgt, θ, 2, 2)

# 6.9
μ = 398600.4418

a_int = 12756.274
a_tgt = 42159.48
θ_i = deg2rad(-20.0)

τ_wait, τ_trans, a_trans, Δv = coplanar_phasing_circular_different(a_tgt, a_int, θ_i, 1)
τ_wait / 60
τ_trans / 60

# 6.10
a_int = 7143.51
i_int = deg2rad(28.5)
Ω_int = deg2rad(45.0)
u_int = deg2rad(15.0)

a_tgt = 42159.4855
i_tgt = deg2rad(0.0)
λ_true0 = deg2rad(200.0)

Δi = i_tgt - i_int

τ, Δt_node, τ_phase, τ_trans, Δv, Δv_phase, Δv_trans_1, Δv_trans_2, a_phase = noncoplanar_phasing(a_tgt, a_int, θ_i, 0, 1, u_int, Ω_int, λ_true0, Δi)
println(τ / 60, " (", Δt_node / 60, ", ", τ_trans / 60, ", ", τ_phase / 60, ")")
println(Δv, " (", Δv_phase, ", ", Δv_trans_1, ", ", Δv_trans_2, ")")
