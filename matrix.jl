
using StaticArrays

const base_i = SVector{3, Float64}([1, 0, 0])
const base_j = SVector{3, Float64}([0, 1, 0])
const base_k = SVector{3, Float64}([0, 0, 1])

function rot1(α)
  sinα = sin(α)
  cosα = cos(α)

  return @SMatrix [
    1 0 0
    0 cosα sinα
    0 -sinα cosα
  ]
end

function rot2(α)
  sinα = sin(α)
  cosα = cos(α)

  return @SMatrix [
    cosα 0 -sinα
    0 1 0
    sinα 0 cosα
  ]
end

function rot3(α)
  sinα = sin(α)
  cosα = cos(α)

  return @SMatrix [
    cosα sinα  0
    -sinα cosα  0
    0 0 1
  ]
end
