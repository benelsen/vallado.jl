include("time.jl")
include("earth_motion.jl")

println( Time.dateToJD(now()) )

### Problem 1
λ = deg2rad(-121.53)
d_lcl = DateTime(2004, 5, 4, 14, 26, 0)

d_utc = d_lcl + Dates.Hour(8)
println(d_utc)

ΔUT1 = 0.046_94
ΔAT = 25.0
d_tai = Time.UTCtoTAI(d_utc, ΔAT)
println(d_tai)

d_ut1 = Time.UTCtoUT1(d_utc, ΔUT1)
println(d_ut1)

jd_ut1 = Time.dateToJD(d_ut1)
println(jd_ut1)

### Problem 2
Φgc = deg2rad(35.05)
λ = deg2rad(-106.40)
h_ell = 1866.03
d = DateTime(1996, 5, 2)

# LHA_eq, GHA_eq, LST, GST
GHA = -λ
LHA = 0

jd = Time.dateToJD(d)

GST = Time.JDtoGST(jd)
LST = GST + λ

### Problem 3
RJ = 71_492
μJ = 1.268e8

TU = sqrt(RJ^3 / μJ)

println(TU)

### Problem 4
# 0

### Problem 5
Φgd = deg2rad(40.0)
λ = deg2rad(-100.0)

ΔUT1 = 0.40233
ΔAT = 26.0

d_ut1 = DateTime(1991, 4, 6, 7, 51, 28, 789)
r_j2000 = [5102.5096; 6123.01152; 6378.1363] # km
v_j2000 = [-4.743_219_6; 0.790_536_6; 5.533_756_19] # km/s

d_utc = Time.UT1toUTC(d_ut1, ΔUT1)
d_tai = Time.UTCtoTAI(d_utc, ΔAT)
d_tdt = Time.TAItoTDT(d_tai)

jd_tdt = Time.dateToJD(d_tdt)

jd_tdb = Time.TDTtoTDB(jd_tdt)
# jd_tdb = 2_448_352.828_085_322_4

P = EarthMotion.precessionMatrix(jd_tdb)
r_ECI_mod = P * r_j2000
@sprintf "%.8f %.8f %.8f" r_ECI_mod...

N = EarthMotion.nutationMatrix(jd_tdb)
r_ECI_tod = N * r_ECI_mod
@sprintf "%.8f %.8f %.8f" r_ECI_tod...

# PQW, SEZ, RSW, NTW, EQW
