include("../initialorbitdetermination.jl")

# Herrick Gibbs Method
r = [SVector(3419.85564, 6019.82602, 2784.60022), SVector(2935.91195, 6326.18324, 2660.59584), SVector(2434.95202, 6597.38674, 2521.52311)]
t = [0, (1 + 16.48/60) / 60 / 24, (2 + 33.04/60) / 60 / 24]
state = herrick_gibbs_method(r, t)
println(state)
elements = stateToElements(state)
println(elements)
