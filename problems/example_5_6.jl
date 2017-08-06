include("../celestial_phenomena.jl")

ut1_d = DateTime(1995, 2, 15, 12, 0, 0)
ut1_jd = dateToJD(ut1_d)

r1 = SVector(0.0, -4464.696, -5102.509)
r2 = SVector(0.0, -5740.323,  3189.068)

sight(r1, r2)

light(r1, ut1_jd)
