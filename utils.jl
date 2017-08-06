import Base.Math.rem2pi

struct NotImplementedException <: Exception end

rem2pi(x) = rem2pi(x, RoundToZero)

function dms_to_d(d, m, s = 0)
  return sign(d) * (abs(d) + m / 60 + s / 3600)
end

function d_to_dms(x)
  d_ = abs(x)
  d = floor(d_)
  m_ = (d_ - d) * 60
  m = floor(m_)
  s = (m_ - m) * 60
  return sign(x) * d, m, s
end

h_to_d(x) = x * 15
d_to_h(x) = x / 15

h_to_rad(x) = x |> h_to_d |> deg2rad
rad_to_h(x) = x |> rad2deg |> d_to_h

hms_to_rad(x) = dms_to_d(x...) |> h_to_d |> deg2rad
rad_to_hms(x) = x |> rad2deg |> d_to_h |> d_to_dms
