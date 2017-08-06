
include("state.jl")

# 1)
M = deg2rad(235.5)
e = 0.4

E = kepler_equation_ell(M, e)
ν = EBH_to_ν(E, e)

# 2)
elements = StateKepler(11_067.790, 0.832853, deg2rad(87.870), deg2rad(227.898), deg2rad(53.38), deg2rad(92.335), now())
convert(StateDelaunay, elements)
convert(StatePoincare, elements)
convert(StateEquinoctial, elements)
convert(StateFlightElements, elements)

# 3)
a = -2_797.425_069
e = 2.8
i = deg2rad(23.0)
ν = deg2rad(249.27)
r = 2_105_124.388

v = sqrt(2μ / r - μ/a)

p = a * (1 - e^2)
r_p = a * (1 - e)

t0 = now()

el0 = StateKepler(p, e, i, deg2rad(0), deg2rad(0), ν, t0)
rv0 = elementsToState(el0)

el1 = StateKepler(p, e, i, deg2rad(0), deg2rad(0), 0, t0)
rv1 = elementsToState(el1)

Δt = find_tof(rv0.r, rv1.r, p)

el1b = StateKepler(p, e, i, deg2rad(0), deg2rad(0), 0, t0 + Dates.Millisecond(round(Δt * 1e3)))
rv1b = elementsToState(el1b)

# 4)

# 5)
t0 = now()
t1 = t0 + Dates.Minute(55)

a = 6678.137
e = 0.000096
p = a * (1 - e^2)

n = sqrt(μ / a^3)

M1 = deg2rad(278.946_88)
M0 = M1 + (-55 * 60) * n

ν0 = EBH_to_ν(kepler_equation_ell(M0, e), e)
ν1 = EBH_to_ν(kepler_equation_ell(M1, e), e) + 2π

el0 = StateKepler(p, e, deg2rad(28.5), 0, 0, ν0, t0)
rv0 = elementsToState(el0)

el1 = StateKepler(p, e, deg2rad(28.5), 0, 0, ν1, t1)
rv1 = elementsToState(el1)

# 6)
Δλ = 360 - 25.0
P = (360 - Δλ) / 15 * 3600
a = cbrt(μ * P^2 / (2π)^2)

# 7)
find_tof([-2574.9533, 4267.0671, 4431.5026], [2700.6738, -4304.5378, -4358.2499], 6681.571)

# 8)

# 9)
# All of them

# 10)
