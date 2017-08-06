include("../initialorbitdetermination.jl")

x_ECI, v_ECI = site_track(
  DateTime(1995, 5, 20, 3, 17, 2),
  GeodeticCoordinate(deg2rad(-104.883), deg2rad(39.007), 2.187),
  604.68, deg2rad(205.6), deg2rad(30.7),
  2.08, deg2rad(0.15), deg2rad(0.17)
)
