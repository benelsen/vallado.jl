include("../initialorbitdetermination.jl")

r0 = SVector(15_945.34, 0., 0.)
r1 = SVector(12_214.838_99, 10_249.467_31, 0.)

v0, = lambert_minumum_energy(r0, r1)
println(@sprintf("Min E: %10.7f %10.7f %10.7f |", v0...))

v0, v1, = lambert_gauss(r0, r1, 76. * 60)
println(@sprintf("Gauss: %10.7f %10.7f %10.7f | %10.7f %10.7f %10.7f", v0..., v1...))

v0, v1 = lambert_univeral(r0, r1, 76. * 60; ψ_low = 0.)
println(@sprintf("Univ:  %10.7f %10.7f %10.7f | %10.7f %10.7f %10.7f", v0..., v1...))

v0, v1 = lambert_batin(r0, r1, 76. * 60)
println(@sprintf("Batin: %10.7f %10.7f %10.7f | %10.7f %10.7f %10.7f", v0..., v1...))
