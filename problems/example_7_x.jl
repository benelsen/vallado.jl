include("../initialorbitdetermination.jl")

using StaticArrays

using Plots
# plotlyjs()
pyplot()
# gr()

r_i1 = SVector(-6518.1083, -2403.8479, -22.1722)
v_i1 = SVector(2.604057, -7.105717, -0.263218)
println(stateToElements(StateVector(r_i1, v_i1, 0.)))

r_t1 = SVector(6697.4756, 1794.5832, 0.)
v_t1 = SVector(-1.962373, 7.323674, 0.)
println(stateToElements(StateVector(r_t1, v_t1, 0.)))

Δv1, Δv2, = target(r_i1, r_t1, v_i1, v_t1, 46.68 * 60; way = :short, prop = false)
println(@sprintf("Δv1: %.3f m/s, Δv2: %.3f m/s", norm(Δv1) * 1e3, norm(Δv2) * 1e3))

begin
  t0, t1 = 0., 500.
  dt = 0.1
  ts = t0:dt:t1

  vs_s = zeros(length(ts), 4)
  for (i, t) in enumerate(ts)
    Δv1, Δv2, Δv_total, hit, way = target(r_i1, r_t1, v_i1, v_t1, t * 60.; way = :short, prop = false)
    vs_s[i, 1] = norm(Δv1)
    vs_s[i, 2] = norm(Δv2)
    vs_s[i, 3] = norm(Δv1) + norm(Δv2)
    vs_s[i, 4] = hit * 10.
  end

  vs_l = zeros(length(ts), 4)
  for (i, t) in enumerate(ts)
    Δv1, Δv2, Δv_total, hit, way = target(r_i1, r_t1, v_i1, v_t1, t * 60.; way = :long, prop = false)
    vs_l[i, 1] = norm(Δv1)
    vs_l[i, 2] = norm(Δv2)
    vs_l[i, 3] = norm(Δv1) + norm(Δv2)
    vs_l[i, 4] = hit * 10.
  end

  l = @layout [a;b]
  plot(
    plot(ts, vs_s, ylims=(0, 35), xlims=(t0, t1)),
    plot(ts, vs_l, ylims=(0, 35), xlims=(t0, t1)),
    layout = l)
  # gui()
end


r_i2 = SVector(5328.7862, 4436.1273, 101.4720)
v_i2 = SVector(-4.864779, 5.816486, 0.240163)
println(stateToElements(StateVector(r_i2, v_i2, 0.)))

r_t2 = SVector(6697.4756, 1794.5831, 0.)
v_t2 = SVector(-1.962372, 7.323674, 0.)
println(stateToElements(StateVector(r_t2, v_t2, 0.)))

Δv1, Δv2, = target(r_i2, r_t2, v_i2, v_t2, 40. * 60; t_m = +1.)
println(@sprintf("Δv1: %.3f km/s, Δv2: %.3f km/s", norm(Δv1) * 1, norm(Δv2) * 1))

Δv1, Δv2, = target(r_i2, r_t2, v_i2, v_t2, 40. * 60; t_m = -1.)
println(@sprintf("Δv1: %.3f km/s, Δv2: %.3f km/s", norm(Δv1) * 1, norm(Δv2) * 1))


begin
  t0, t1 = 0., 500.
  dt = 0.1
  ts = t0:dt:t1

  vs_s = zeros(length(ts), 4)
  for (i, t) in enumerate(ts)
    Δv1, Δv2, Δv_total, hit, way = target(r_i2, r_t2, v_i2, v_t2, t * 60.; way = :short)
    vs_s[i, 1] = norm(Δv1)
    vs_s[i, 2] = norm(Δv2)
    vs_s[i, 3] = norm(Δv1) + norm(Δv2)
    vs_s[i, 4] = hit * 10.
  end

  vs_l = zeros(length(ts), 4)
  for (i, t) in enumerate(ts)
    Δv1, Δv2, Δv_total, hit, way = target(r_i2, r_t2, v_i2, v_t2, t * 60.; way = :long)
    vs_l[i, 1] = norm(Δv1)
    vs_l[i, 2] = norm(Δv2)
    vs_l[i, 3] = norm(Δv1) + norm(Δv2)
    vs_l[i, 4] = hit * 10.
  end

  l = @layout [a;b]
  plot(
    plot(ts, vs_s, ylims=(0, 35), xlims=(t0, t1)),
    plot(ts, vs_l, ylims=(0, 35), xlims=(t0, t1)),
    layout = l)
end


begin
  t0, t1 = 0., 250.
  dt = 0.1
  ts = t0:dt:t1

  vs = zeros(length(ts), 5)
  for (i, t) in enumerate(ts)

    Δv1, Δv2, Δv_total, hit, way = target(r_i2, r_t2, v_i2, v_t2, t * 60., 0.; way = :both)

    vs[i, 1] = norm(Δv1)
    vs[i, 2] = norm(Δv2)
    vs[i, 3] = Δv_total
    vs[i, 4] = hit * 10.
    vs[i, 5] = way === :short ? 0 : 1
  end

  plot(ts, vs, ylims=(0, 25), xlims=(t0, t1))
end

begin
  ts_wait = 0. : 1. : 750.
  ts_transfer = 0. : 1. : 300.

  vs = zeros(length(ts_wait), length(ts_transfer))
  for (i, t_wait) in enumerate(ts_wait)
    for (j, t_transfer) in enumerate(ts_transfer)
      Δv1, Δv2, Δv_total, hit, way = target(r_i2, r_t2, v_i2, v_t2, t_transfer * 60., t_wait * 60.; way = :both)
      vs[i, j] = hit ? 100. : Δv_total
    end
  end

  i, j = ind2sub(vs, indmin(vs))
  println("Wait: ", ts_wait[i], "m Transfer: ", ts_transfer[j], "m Δv: ", vs[i, j], "km/s")
  contour(ts_transfer, ts_wait, vs, levels = 0.:.1:25., fill=true, xlabel = "transfer", ylabel = "wait", c=:viridis)
end

# Minimum Case
begin
  r_i3 = SVector(4965.2226, -1504.1795, 4617.2886)
  v_i3 = SVector(3.760556, -3.882606, -5.310687)
  r_t3 = SVector(5055.2385, -1599.4069, 4483.9878)
  v_t3 = SVector(3.612200, -3.836100, -5.443000)

  ts_wait = 0. : 1. : 250.
  ts_transfer = 0. : 1. : 250.

  vs = zeros(length(ts_wait), length(ts_transfer))
  for (i, t_wait) in enumerate(ts_wait)
    for (j, t_transfer) in enumerate(ts_transfer)
      Δv1, Δv2, Δv_total, hit, way = target(r_i3, r_t3, v_i3, v_t3, t_transfer * 60., t_wait * 60.; way = :both)
      vs[i, j] = hit ? 100. : Δv_total
    end
  end

  i, j = ind2sub(vs, indmin(vs))
  println("Wait: ", ts_wait[i], "m Transfer: ", ts_transfer[j], "m Δv: ", vs[i, j], "km/s")
  contour(ts_transfer, ts_wait, vs, levels = 0.:.1:25., fill=true, xlabel = "transfer", ylabel = "wait", c=:viridis)
end

# Minimum Case (Zoom)
begin
  r_i4 = SVector(4965.2226, -1504.1795, 4617.2886)
  v_i4 = SVector(3.760556, -3.882606, -5.310687)
  r_t4 = SVector(5055.2385, -1599.4069, 4483.9878)
  v_t4 = SVector(3.612200, -3.836100, -5.443000)

  ts_wait = 0. : 0.5 : 250.
  ts_transfer = 0. : 0.5 : 100.

  vs4 = zeros(length(ts_wait), length(ts_transfer))
  for (i, t_wait) in enumerate(ts_wait)
    for (j, t_transfer) in enumerate(ts_transfer)
      Δv1, Δv2, Δv_total, hit, way = target(r_i4, r_t4, v_i4, v_t4, t_transfer * 60., t_wait * 60.; way = :both)
      vs4[i, j] = hit ? 100. : Δv_total
    end
  end

  i, j = ind2sub(vs4, indmin(vs4))
  println("Wait: ", ts_wait[i], "m Transfer: ", ts_transfer[j], "m Δv: ", vs4[i, j], "km/s")
  contour(ts_transfer, ts_wait, vs4, levels = 0.:.05:2.0, fill=true, xlabel = "transfer", ylabel = "wait", c=:viridis, title="Min Case (Zoom)")
end

# Maximum Case (Zoom)
begin
  r_i5 = SVector(-6175.1034, 2757.0706, 1626.6556)
  v_i5 = SVector(2.376641, 1.139677, 7.078097)
  r_t5 = SVector(-6078.007289, 2796.641859, 1890.7135)
  v_t5 = SVector(2.654700, 1.018600, 7.015400)

  ts_wait = 0. : 0.5 : 250.
  ts_transfer = 0. : 0.5 : 100.

  vs5 = zeros(length(ts_wait), length(ts_transfer))
  for (i, t_wait) in enumerate(ts_wait)
    for (j, t_transfer) in enumerate(ts_transfer)
      Δv1, Δv2, Δv_total, hit, way = target(r_i5, r_t5, v_i5, v_t5, t_transfer * 60., t_wait * 60.; way = :both)
      vs5[i, j] = hit ? 100. : Δv_total
    end
  end

  i, j = ind2sub(vs5, indmin(vs5))
  println("Wait: ", ts_wait[i], "m Transfer: ", ts_transfer[j], "m Δv: ", vs5[i, j], "km/s")
  contour(ts_transfer, ts_wait, vs5, levels = 0.:.05:3.0, fill=true, xlabel = "transfer", ylabel = "wait", c=:viridis, title="Max Case (Zoom)")
end
