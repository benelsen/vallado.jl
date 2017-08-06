include("../initialorbitdetermination.jl")

# Gibbs method
r = [SVector(0., 0., 6378.137), SVector(0., -4464.696, -5102.509), SVector(0., 5740.323, 3189.068)]
state = gibbs_method(r, now())
println(state)
elements = stateToElements(state)
println(elements)
