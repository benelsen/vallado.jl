using Interpolations

include("time.jl")

eop_c04_itrf14_iau2000_data, eop_c04_itrf14_iau2000_header  = readdlm("data/eopc04_14_IAU2000.62-now.csv", ';'; header = true)
eop_c04_itrf14_iau1980_data, eop_c04_itrf14_iau1980_header  = readdlm("data/eopc04_14_IAU1980.62-now.csv", ';'; header = true)
eop_c04_itrf08_iau2000_data, eop_c04_itrf08_iau2000_header  = readdlm("data/eopc04_08_IAU2000.62-now.csv", ';'; header = true)
eop_c04_itrf08_iau1980_data, eop_c04_itrf08_iau1980_header  = readdlm("data/eopc04_08_IAU1980.62-now.csv", ';'; header = true)

eop_c04_itrf14_iau2000_header = eop_c04_itrf14_iau2000_header[[collect(1:4); collect(6:9); collect(11:14); collect(20:23)]]
eop_c04_itrf14_iau2000_data = convert(Array{Float64}, eop_c04_itrf14_iau2000_data[:, [collect(1:4); collect(6:9); collect(11:14); collect(20:23)]])
eop_c04_itrf14_iau2000 = interpolate(eop_c04_itrf14_iau2000_data[:, collect(5:16)], (BSpline(Linear()), NoInterp()), OnGrid())

eop_c04_itrf08_iau1980_header = eop_c04_itrf08_iau1980_header[[collect(1:4); collect(6:9); collect(11:14); collect(16:19)]]
eop_c04_itrf08_iau1980_data = convert(Array{Float64}, eop_c04_itrf08_iau1980_data[:, [collect(1:4); collect(6:9); collect(11:14); collect(16:19)]])
eop_c04_itrf08_iau1980 = interpolate(eop_c04_itrf08_iau1980_data[:, collect(5:16)], (BSpline(Linear()), NoInterp()), OnGrid())

struct EOP
    jd :: Real
    x :: Real
    σ_x :: Real
    y :: Real
    σ_y :: Real
    ΔUT1 :: Real
    σ_ΔUT1 :: Real
    ΔLOD :: Real
    σ_ΔLOD :: Real
    ΔX :: Real
    σ_ΔX :: Real
    ΔY :: Real
    σ_ΔY :: Real
end

EOP(jd, x, y, ΔUT1, ΔLOD, ΔX, ΔY) = EOP(jd, x, 0, y, 0, ΔUT1, 0, ΔLOD, 0, ΔX, 0, ΔY, 0)

struct EOP_iau1980
    jd :: Real
    x :: Real
    σ_x :: Real
    y :: Real
    σ_y :: Real
    ΔUT1 :: Real
    σ_ΔUT1 :: Real
    ΔLOD :: Real
    σ_ΔLOD :: Real
    Δψ :: Real
    σ_Δψ :: Real
    Δϵ :: Real
    σ_Δϵ :: Real
end

EOP_iau1980(jd, x, y, ΔUT1, ΔLOD, Δψ, Δϵ) = EOP_iau1980(jd, x, 0, y, 0, ΔUT1, 0, ΔLOD, 0, Δψ, 0, Δϵ, 0)


function get_c04_eop(jd; frame = 2014, theory = 2000)
    mjd = jd - 2_400_000.5
    idx = mjd - eop_c04_itrf14_iau2000_data[1, 1] + 1

    return EOP(
    jd,
    eop_c04_itrf14_iau2000[idx,  1],
    eop_c04_itrf14_iau2000[idx,  2],
    eop_c04_itrf14_iau2000[idx,  3],
    eop_c04_itrf14_iau2000[idx,  4],
    eop_c04_itrf14_iau2000[idx,  5],
    eop_c04_itrf14_iau2000[idx,  6],
    eop_c04_itrf14_iau2000[idx,  7],
    eop_c04_itrf14_iau2000[idx,  8],
    eop_c04_itrf14_iau2000[idx,  9],
    eop_c04_itrf14_iau2000[idx, 10],
    eop_c04_itrf14_iau2000[idx, 11],
    eop_c04_itrf14_iau2000[idx, 12]
    )
end

function get_c04_eop_1980(jd; frame = 2008, theory = 1980)
    mjd = jd - 2_400_000.5
    idx = mjd - eop_c04_itrf08_iau1980_data[1, 1] + 1

    return EOP_iau1980(
    jd,
    eop_c04_itrf08_iau1980[idx,  1],
    eop_c04_itrf08_iau1980[idx,  2],
    eop_c04_itrf08_iau1980[idx,  3],
    eop_c04_itrf08_iau1980[idx,  4],
    eop_c04_itrf08_iau1980[idx,  5],
    eop_c04_itrf08_iau1980[idx,  6],
    eop_c04_itrf08_iau1980[idx,  7],
    eop_c04_itrf08_iau1980[idx,  8],
    eop_c04_itrf08_iau1980[idx,  9],
    eop_c04_itrf08_iau1980[idx, 10],
    eop_c04_itrf08_iau1980[idx, 11],
    eop_c04_itrf08_iau1980[idx, 12]
    )
end

# get_c04_eop(dateToJD(DateTime(2017, 4, 3, 12)))
