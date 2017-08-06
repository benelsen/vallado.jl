include("../celestial_phenomena.jl")

ut1_d = DateTime(1994, 4, 28, 0, 0, 0)
ut1_jd = dateToJD(ut1_d)

r, λ_ecl, ϕ_ecl = moon_ecliptic(ut1_jd)
println("r: ", (r), "λ: ", rad2deg(λ_ecl), " ϕ: ", rad2deg(ϕ_ecl))

m = moon_vector(ut1_jd)
println(@sprintf("%.8f %.8f %.8f", m...))

r, α, δ = moon_radec(ut1_jd)
println("α: ", rad2deg(α), " δ: ", rad2deg(δ))
