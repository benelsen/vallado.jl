include("../nutation.jl")

tt_jc = 0.042_623_631_9
ut1_jd = 2_453_101.827_406_783

xp = as_to_rad(-0.140_682)
yp = as_to_rad( 0.333_309)

r_ITRF = [-1033.4793830, 7901.2952754, 6380.3565958]

(X, Y, s) = precession_nutation_IAU_2000_2006_CIO_XYs(tt_jc)
dX = as_to_rad(-0.000205)
dY = as_to_rad(-0.000136)

rad_to_as(X) - 80.531880 # 2.318754084740249e-6
rad_to_as(Y) - 7.273921 # -2.2267244936813313e-5
rad_to_as(s) # OK

Q = precession_nutation_IAU_2000_2006_CIO(tt_jc; dX = dX, dY = dY)

R_CIO = CIRS_to_TIRS_CIO(ut1_jd)

W_CIO = polar_motion(tt_jc, xp, yp)

T_CIO = Q * transpose(R_CIO) * transpose(W_CIO)

r_GCRF_CIO = T_CIO * r_ITRF
println(r_GCRF_CIO)


N, P, B, R_EQU = precession_nutation_IAU_2000_2006_EQU(tt_jc, ut1_jd)

W_EQU = polar_motion(tt_jc, xp, yp)

T_EQU = transpose(B) * transpose(P) * transpose(N) * transpose(R_EQU) * transpose(W_EQU)

r_GCRF_EQU = T_EQU * r_ITRF
println(r_GCRF_EQU)
println(r_GCRF_EQU - r_GCRF_CIO)


N, P_FW, B, R = precession_nutation_IAU_2000_2006_EQU_FW(tt_jc, ut1_jd)

T_EQU_FW = transpose(P_FW) * transpose(N) * transpose(R_EQU) * transpose(W_EQU)

r_GCRF_EQU_FW = T_EQU_FW * r_ITRF
println(r_GCRF_EQU_FW)
println(r_GCRF_EQU_FW - r_GCRF_CIO)
