include("../celestial_phenomena.jl")

ut1_d = DateTime(1994, 5, 20, 20, 0, 0)
ut1_jd = dateToJD(ut1_d)

state = planet_ecliptic(ut1_jd)

AU = 149597870.700 # [km]
μ = 1.32712440018e20 / 1e9 # [km/s^2]
TU = sqrt(AU^3 / μ)

r_J2000 = rot1(-obliquity_IAU_1980(julianYears(ut1_jd))) * state.r # [AU]
v_J2000 = rot1(-obliquity_IAU_1980(julianYears(ut1_jd))) * state.v # [AU/TU]

r_J2000 * AU # [km]
v_J2000 * AU / TU # [km/s]
