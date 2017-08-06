include("../initialorbitdetermination.jl")

obs_raw = [
   1 2012 08 20 11 32 28.00 333.028738 -2.022317
   2 2012 08 20 11 36 28.00 345.235515  7.648921
   3 2012 08 20 11 40 28.00   0.939913 18.667717
   4 2012 08 20 11 44 28.00  21.235600 29.086871
   5 2012 08 20 11 48 28.00  45.025748 35.664741
   6 2012 08 20 11 52 28.00  67.886655 36.996583
   7 2012 08 20 11 56 28.00  86.208078 34.719667
   8 2012 08 20 12 00 28.00  99.845522 30.928387
   9 2012 08 20 12 04 28.00 110.078585 26.767438
  10 2012 08 20 12 08 28.00 118.058822 22.680214
  11 2012 08 20 12 12 28.00 124.552101 18.801899
  12 2012 08 20 12 16 28.00 130.038810 15.154469
  13 2012 08 20 12 20 28.00 134.823380 11.722965
  14 2012 08 20 12 24 28.00 139.104912  8.483073
  15 2012 08 20 12 26 28.00 141.101028  6.927305]

obs_all = [
  begin
    eop = EOP(0, as_to_rad(0.137495), as_to_rad(0.342416), -0.609641, 0, 0, 0)

    gd = GeodeticCoordinate(deg2rad(-110.), deg2rad(40.), 2.)
    r_site_ECEF = geodetic_to_cartesian(gd)

    utc_dt = DateTime(obs_raw[i, 2], obs_raw[i, 3], obs_raw[i, 4], obs_raw[i, 5], obs_raw[i, 6], obs_raw[i, 7])
    ut1_dt = UTCtoUT1(utc_dt, eop.ΔUT1)
    ut1_jd = dateToJD(ut1_dt)
    tt_jd = dateToJD(TAItoTT(UTCtoTAI(utc_dt, 35)))

    r_site = ITRF_to_GCRF(tt_jd, ut1_jd; x_p = eop.x, y_p = eop.y, ΔX = eop.ΔX, ΔY = eop.ΔY, LOD = eop.ΔLOD)(r_site_ECEF)

    ObservationsRhoRaDec(
      utc_dt,
      r_site,
      deg2rad(obs_raw[i, 8]), deg2rad(obs_raw[i, 9])
    )
  end for i in 1:size(obs_raw, 1)]


# Laplace method
state = laplace_method(obs_all[[3,5,6]])
println(state)
elements = stateToElements(state)
println(elements)

# Gauss method
r = gauss_method(obs_all[[3,5,6]])

# Double-r method
state = double_r(obs_all[[3,5,6]]; ϵ = 1e-9)
println(state)
elements = stateToElements(state)
println(elements)
