include("../celestial_phenomena.jl")

ut1_d = DateTime(2006, 4, 2, 0, 0, 0)
ut1_jd = dateToJD(ut1_d)

s = sun_vector(ut1_jd)
println(@sprintf("%.8f %.8f %.8f", s...))

r, α, δ = sun_radec(ut1_jd)
println("α: ", rad2deg(α), " δ: ", rad2deg(δ))
