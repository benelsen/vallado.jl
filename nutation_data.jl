using StaticArrays

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.5.1 - Table 5.1a
# γ  l  l′ F  D  Ω Period [d]   ΔXp s    c ΔYp s     c
const polar_motion_libration_coeffs = readdlm("data/polar_motion_libration.dlm")
const polar_motion_libration_secular_coeffs = [-3.8, -4.3]

# IERS Conventions 2010 - IERS Technical Note No. 36 - Gérard Petit, Brian Luzum
# 2012-08-10
# Section 5.5.3 - Table 5.1b
# γ  l l′  F  D  Ω  Period [d]  UT1 s     c LOD s     c
const UT1_LOD_libration_coeffs = readdlm("data/UT1_LOD_libration.dlm")

# nutation_2000A_X.dlm (tab5.2a.txt)
# 2012-08-10
# Section 5.x.y - Table 5.2a
#   X = polynomial part + non-polynomial part
# - 16617. + 2004191898. t - 429782.9 t^2 - 198618.34 t^3 + 7.578 t^4 + 5.9285 t^5
#   Sum_i[a_{s,0})_i * sin(ARG) + a_{c,0})_i * cos(ARG)]
# + Sum_i,j=1,4 [a_{s,j})_i * sin(ARG) + a_{c,j})_i * cos(ARG)] * t^j
# i    a_{s,j})_i      a_{c,j})_i    l    l'   F    D   Om L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
const nutation_2000A_X_coeffs = readdlm("data/nutation_2000A_X.dlm")

# nutation_2000A_Y.dlm (tab5.2b.txt)
# 2012-08-10
# Section 5.x.y - Table 5.2b
#   Y = polynomial part + non-polynomial part
# - 6951. - 25896. t - 22407274.7 t^2 + 1900.59 t^3 + 1112.526 t^4 + 0.1358 t^5
#   Sum_i[b_{c,0})_i * cos(ARG) + b_{s,0})_i * sin(ARG)]
# + Sum_i,j=1,4 [b_{c,j})_i * cos(ARG) + b_{s,j})_i * sin(ARG)] * t^j
# i    b_{s,j})_i      b_{c,j})_i    l    l'   F    D   Om L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
const nutation_2000A_Y_coeffs = readdlm("data/nutation_2000A_Y.dlm")

# nutation_2000A_s.dlm (tab5.2d.txt)
# 2012-08-10
# Section 5.x.y - Table 5.2d
#   s + XY/2 = polynomial part + non-polynomial part
# 94.0 + 3808.65 t - 122.68 t^2 - 72574.11 t^3 + 27.98 t^4 + 15.62 t^5
#   Sum_i[C_{s,0})_i * sin(ARG) + C_{c,0})_i * cos(ARG)]
# + Sum_i,j=1,4 [C_{s,j})_i * sin(ARG) + C_{c,j})_i * cos(ARG)] * t^j
# i    C_{s,j})_i      C_{c,j})_i    l    l'   F    D   Om L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
const nutation_2000A_s_coeffs = readdlm("data/nutation_2000A_s.dlm")

# nutation_2000A_GST.dlm (tab5.2e.txt)
# 2012-08-10
# Section 5.x.y - Table 5.2e
# GST = ERA(UT1) + polynomial part  + DeltaPsi*cos(epsilon_A) + non-polynomial additional part
#
#   ERA(UT1) = 2*Pi*(0.7790572732640 + 1.00273781191135448.T_u)
#   where T_u = Julian UT1 date - 2451545.0, and UT1 = UTC + (UT1 - UTC)
#
#   0.014506 + 4612.156534 t + 1.3915818 t^2 - 0.00000044 t^3 - 0.000029956 t^4 - 0.0000000368 t^5
#
#   DeltaPsi*cos(epsilon_A)
#   IAU 200A_R06 expression for the classical "equation of the equinoxes"
#
#   Sum_i[C'_{s,0})_i * sin(ARG) + C'_{c,0})_i * cos(ARG)]
# + Sum_i[C'_{s,1})_i * sin(ARG) + C'_{c,1})_i * cos(ARG)] * t
# i   C'_{s,j})_i     C'_{c,j})_i    l    l'   F    D   Om L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
const nutation_2000A_GST_coeffs = readdlm("data/nutation_2000A_GST.dlm")

# nutation_2006_2000A_diff_X.dlm (tab5.2f.txt)
# 2012-08-10
# Section 5.x.y - Table 5.2f
#   X_IAU 2006 - X_IAU2000A = polynomial part + non-polynomial part
# 155. t - 2564.0 t^2 + 2.20 t^3 + 53.63 t^4
# i      D_{s,j})_i     D_{c,j})_i   l    l'   F    D   Om
const nutation_2006_2000A_diff_X_coeffs = readdlm("data/nutation_2006_2000A_diff_X.dlm")

# nutation_2006_2000A_diff_Y.dlm (tab5.2f.txt)
# 2012-08-10
# Section 5.x.y - Table 5.2f
#   Y_IAU 2006 - Y_IAU2000A = polynomial part + non-polynomial part
# - 0.2 - 514.  t - 23.7  t^2 + 58.31 t^3 - 0.53 t^4 - 0.85 t^5
# i      D_{s,j})_i     D_{c,j})_i   l    l'   F    D   Om
const nutation_2006_2000A_diff_Y_coeffs = readdlm("data/nutation_2006_2000A_diff_Y.dlm")

# nutation_2006_2000A_diff_Dpsi.dlm (tab5.2f.txt)
# 2012-08-10
# Section 5.x.y - Table 5.2f
#   Dpsi(IAU 2000A_R06 - IAU 2000A)
# i    D_{s,j})_i      D_{c,j})_i    l    l'   F    D   Om
const nutation_2006_2000A_diff_Dpsi_coeffs = readdlm("data/nutation_2006_2000A_diff_Dpsi.dlm")

# nutation_2006_2000A_diff_Deps.dlm (tab5.2f.txt)
# 2012-08-10
# Section 5.x.y - Table 5.2f
#   Deps(IAU 2000A_R06 - IAU 2000A)
# i    D_{s,j})_i      D_{c,j})_i    l    l'   F    D   Om
const nutation_2006_2000A_diff_Deps_coeffs = readdlm("data/nutation_2006_2000A_diff_Deps.dlm")

# nutation_2000_R06_Dpsi.dlm (tab5.3a.txt)
# 2012-08-10
# Section 5.x.y - Table 5.3a
#   Sum_i[A_i * sin(ARG) + A"_i * cos(ARG)]
# + Sum_i[A'_i * sin(ARG) + A"'_i * cos(ARG)] * t           (see Chapter 5, Eq. (35))
# i        A_i             A"_i     l    l'   F    D    Om  L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
# i        A'_i            A"'_i    l    l'   F    D    Om  L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
const nutation_2000_R06_Dpsi_coeffs = readdlm("data/nutation_2000_R06_Dpsi.dlm")


# nutation_2000_R06_Deps.dlm (tab5.3b.txt)
# 2012-08-10
# Section 5.x.y - Table 5.3b
#   Sum_i[B_i * cos(ARG) + B"_i * sin(ARG)]
# + Sum_i[B'_i * cos(ARG) + B"'_i * sin(ARG)] * t      (see Chapter 5, Eq. (35))
# i         B"_i            B_i      l    l'   F    D   Om  L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
# i         B"'_i           B'_i     l    l'   F    D   Om  L_Me L_Ve  L_E L_Ma  L_J L_Sa  L_U L_Ne  p_A
const nutation_2000_R06_Deps_coeffs = readdlm("data/nutation_2000_R06_Deps.dlm")


nutation_values = readdlm("data/nut_IAU1980.dat")

const nutation_coeffs_arguments = SMatrix{106, 5}(nutation_values[:, 1:5])
const nutation_coeffs_longitude = SMatrix{106, 2}(nutation_values[:, 7:8])
const nutation_coeffs_obliquity = SMatrix{106, 2}(nutation_values[:, 9:10])
