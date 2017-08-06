include("../nutation.jl")

tt_jc = 0.042_623_631_9
ut1_jd = 2_453_101.827_406_783
ut1_jc = julianYears(ut1_jd)

xp = as_to_rad(-0.140_682)
yp = as_to_rad( 0.333_309)

δ∆ψ_1980 = as_to_rad(-0.052_195)
δ∆ϵ_1980 = as_to_rad(-0.003_875)

r_ITRF = [-1033.479_383_00, 7901.295_275_40, 6380.356_595_80]

W = polar_motion_IAU_1980(xp, yp)

r_PEF = transpose(W) * r_ITRF
println("PEF: ", @sprintf("%.8f %.8f %.8f", r_PEF...))

rad2deg(GMST_IAU_1980(ut1_jc))

Δψ, Δϵ = nutation_IAU_1980_EQU_ang(tt_jc)

rad2deg(Δψ)
rad2deg(Δϵ)

ϵ_A = obliquity_IAU_1980(tt_jc)
rad2deg(GAST_IAU_1982(ut1_jc, Δψ, ϵ_A, tt_jc))

N, P, R = precession_nutation_IAU_1980(tt_jc, ut1_jc)

r_TOD = transpose(R) * r_PEF
println("TOD: ", @sprintf("%.8f %.8f %.8f", r_TOD...))

r_MOD = transpose(N) * r_TOD
println("MOD: ", @sprintf("%.8f %.8f %.8f", r_MOD...))

r_J2000 = transpose(P) * r_MOD
println("J2000: ", @sprintf("%.8f %.8f %.8f", r_J2000...))

T = transpose(P) * transpose(N) * transpose(R) * transpose(W)

r_J2000 = T * r_ITRF
println("J2000: ", @sprintf("%.8f %.8f %.8f", r_J2000...))
