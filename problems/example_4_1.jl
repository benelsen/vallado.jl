include("../coordinate_systems.jl")
include("../utils.jl")
include("../nutation.jl")

utc_dt = DateTime(1994, 5, 14, 13, 11, 20, 599)

# IAU2000A/2006
ΔAT = 28.0
ΔUT = -0.127940

x_p = as_to_rad(0.189443)
y_p = as_to_rad(0.306064)

ΔX = as_to_rad(0.000187)
ΔY = as_to_rad(0.000039)

# IAU1980
x_p = as_to_rad(0.189015)
y_p = as_to_rad(0.305910)

ΔUT = -0.127915
LOD = 0.0021632

Δψ = as_to_rad(-0.015955)
Δϵ = as_to_rad(-0.008602)


ut1_dt = UTCtoUT1(utc_dt, ΔUT)
tai_dt = UTCtoTAI(utc_dt, ΔAT)
tt_dt = TAItoTT(tai_dt)

ut1_jd = dateToJD(ut1_dt)
tt_jd = dateToJD(tt_dt)

λ_gd = deg2rad(-104.883)
ϕ_gd = deg2rad(39.007)
h_ell = 2.19456

gd = GeodeticCoordinate(λ_gd, ϕ_gd, h_ell)

r_site_ITRF = geodetic_to_cartesian(gd)
v_site_ITRF = SVector{3}(0.0, 0.0, 0.0)
println(@sprintf("%.7f %.7f %.7f", r_site_ITRF...))
println(@sprintf("%.4f %.4f %.4f", v_site_ITRF...))

T_IG = ITRF_to_GCRF(tt_jd, ut1_jd; x_p = x_p, y_p = y_p, ΔX = ΔX, ΔY = ΔY)
r_site_GCRF, v_site_GCRF = T_IG(r_site_ITRF, v_site_ITRF)
println(@sprintf("%.7f %.7f %.7f", r_site_GCRF...))
println(@sprintf("%.7f %.7f %.7f", v_site_GCRF...))


r = 4_437_725_220.51
α = deg2rad(h_to_d(dms_to_d(19, 39, 57.395)))
δ = deg2rad(dms_to_d(-20, 49, 24.58))

ṙ = -25.530_330_94
α̇ = deg2rad(-122.44e-9)
δ̇ = deg2rad( -17.94e-9)

r_ECI, v_ECI = radec_to_state(r, α, δ, ṙ, α̇, δ̇)
println(@sprintf("%.1f %.1f %.1f", r_ECI...))
println(@sprintf("%.4f %.4f %.4f", v_ECI...))

T_GI = GCRF_to_ITRF(tt_jd, ut1_jd; x_p = x_p, y_p = y_p, ΔX = ΔX, ΔY = ΔY)
r_ITRF, v_ITRF = T_GI(r_ECI, v_ECI)
println(@sprintf("%.7f %.7f %.7f", r_ITRF...))
println(@sprintf("%.7f %.7f %.7f", v_ITRF...))


geocentric = state_to_radec(r_ECI, v_ECI)
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [geocentric[1], rad2deg(mod2pi(geocentric[2])), rad2deg(geocentric[3]), geocentric[4], rad2deg.(geocentric[5:6])...]...))

r_ECI_t, v_ECI_t = geocentric_to_topocentric(r_site_GCRF, v_site_GCRF, r_ECI, v_ECI)
topocentric = state_to_radec(r_ECI_t, v_ECI_t)
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [topocentric[1], rad2deg(mod2pi(topocentric[2])), rad2deg(topocentric[3]), topocentric[4], rad2deg.(topocentric[5:6])...]...))

r_SEZ, v_SEZ = ECEF_to_SEZ(gd, r_site_ITRF, r_ITRF, v_ITRF)
horizon = SEZ_to_azel(r_SEZ, v_SEZ)
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [horizon[1], rad2deg(mod2pi(horizon[2])), rad2deg(horizon[3]), horizon[4], rad2deg.(horizon[5:6])...]...))

ecliptic = radec_to_ecliptic(geocentric[2], geocentric[3], geocentric[5], geocentric[6])
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [geocentric[1], rad2deg(mod2pi(ecliptic[1])), rad2deg(ecliptic[2]), geocentric[4], rad2deg.(ecliptic[3:4])...]...))

radec = ecliptic_to_radec(ecliptic...)
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [geocentric[1], rad2deg(mod2pi(radec[1])), rad2deg(radec[2]), geocentric[4], rad2deg.(radec[3:4])...]...))


r_ECI = SVector{3}(5_036.736_529, -10_806.660_797, -4_534.633_784)
v_ECI = SVector{3}(2.684_385_5, -5.759_592_0, -2.416_809_3)
println(@sprintf("%.1f %.1f %.1f", r_ECI...))
println(@sprintf("%.4f %.4f %.4f", v_ECI...))

T_GI = GCRF_to_ITRF(tt_jd, ut1_jd; x_p = x_p, y_p = y_p, ΔX = ΔX, ΔY = ΔY)
r_ITRF, v_ITRF = T_GI(r_ECI, v_ECI)
println(@sprintf("%.7f %.7f %.7f", r_ITRF...))
println(@sprintf("%.7f %.7f %.7f", v_ITRF...))


geocentric = state_to_radec(r_ECI, v_ECI)
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [geocentric[1], rad2deg(mod2pi(geocentric[2])), rad2deg(geocentric[3]), geocentric[4], rad2deg.(geocentric[5:6])...]...))

r_ECI_t, v_ECI_t = geocentric_to_topocentric(r_site_GCRF, v_site_GCRF, r_ECI, v_ECI)
topocentric = state_to_radec(r_ECI_t, v_ECI_t)
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [topocentric[1], rad2deg(mod2pi(topocentric[2])), rad2deg(topocentric[3]), topocentric[4], rad2deg.(topocentric[5:6])...]...))

r_SEZ, v_SEZ = ECEF_to_SEZ(gd, r_site_ITRF, r_ITRF, v_ITRF)
horizon = SEZ_to_azel(r_SEZ, v_SEZ)
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [horizon[1], rad2deg(mod2pi(horizon[2])), rad2deg(horizon[3]), horizon[4], rad2deg.(horizon[5:6])...]...))

ecliptic = radec_to_ecliptic(geocentric[2], geocentric[3], geocentric[5], geocentric[6])
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [geocentric[1], rad2deg(mod2pi(ecliptic[1])), rad2deg(ecliptic[2]), geocentric[4], rad2deg.(ecliptic[3:4])...]...))

radec = ecliptic_to_radec(ecliptic...)
println(@sprintf("%.3f %.7f %.7f | %.7f %.12f %.12f", [geocentric[1], rad2deg(mod2pi(radec[1])), rad2deg(radec[2]), geocentric[4], rad2deg.(radec[3:4])...]...))
