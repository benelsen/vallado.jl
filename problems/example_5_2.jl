include("../celestial_phenomena.jl")

ut1_d = DateTime(1996, 3, 23, 12, 0, 0)
ut1_jd = dateToJD(ut1_d)

λ = deg2rad(0.0)
ϕ = deg2rad(85.0)
ut_sr, ut_ss = sunrise_sunset(ut1_jd, λ, ϕ)

d_to_dms(d_to_h(rad2deg(ut_sr)))
d_to_dms(d_to_h(rad2deg(ut_ss)))
