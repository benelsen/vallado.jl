function lagrange(ys, xs, t)
  n = length(xs)

  L = 0.0
  for i in 1:n

    a = 1.0
    for k in 1:n
      if k === i
        continue
      end
      a *= (t - xs[k]) / (xs[i] - xs[k])
    end

    L += ys[i] * a
  end

  return L
end

function lagrange_d1(ys, xs, t)
  n = length(xs)

  L = 0.0
  for i in 1:n

    a = 0.
    b = 1.
    for k in 1:n
      if k === i
        continue
      end

      c = 1.
      for l in 1:n
        if l === i || l === k
          continue
        end
        c *= t - xs[l]
      end

      a += c
      b *= xs[i] - xs[k]
    end

    L += ys[i] * (a / b)
  end

  return L
end

function lagrange_d2(ys, xs, t)
  n = length(xs)

  L = 0.0
  for i in 1:n

    a = 0.0
    for l in 1:n
      if l === i
        continue
      end

      b = 0.0
      for m in 1:n
        if m === l || m === i
          continue
        end

        c = 1.0
        for k in 1:n
          if k === m || k === l || k === i
            continue
          end

          c *= (t - xs[k]) / (xs[i] - xs[k])
        end

        b += 1 / (xs[i] - xs[m]) * c
      end

      a += 1 / (xs[i] - xs[l]) * b
    end

    L += ys[i] * a
  end

  return L
end
