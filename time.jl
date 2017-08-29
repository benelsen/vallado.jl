const ωE = 0.000_072_921_158_553_0

function JDtoGMST(jd)
    T = julianYears(jd)
    θGMST = mod(@evalpoly(T, 67_310.548_41, (876_600 * 60^2 + 8_640_184.812_866), 0.093_104, -6.2e-6), 86400) # [s]
    return θGMST / 86400 * 2π
    # θGMST = @evalpoly(T, 4.894961212823059, 230121.67531542323, 6.770713944903336e-6, -4.5087672343186846e-10)
    # return mod2pi(θGMST)
end

function GMSTtoLST(θGMST, λ)
    θLST = θGMST + λ
end

# Meeus (p 61ff)
function dateToJD(date::DateTime, julian = false)
    month = Dates.month(date)
    da = fld(14 - month, 12)
    y = Dates.year(date) - da
    m = month + 12 * da
    a = fld(y, 100)
    b = (2 - a + fld(a, 4)) * !julian
    t = (((Dates.millisecond(date) * 1e-3 + Dates.second(date)) / 60 + Dates.minute(date)) / 60 + Dates.hour(date)) / 24
    return floor(365.25 * (y + 4716)) + floor(30.6001 * (m + 1)) + Dates.day(date) + b - 1524.5 + t
end

# Meeus (p 63ff)
function JDtoDate(jd)
    jd += 0.5
    a, f = fldmod(jd, 1)

    # First day of the Gregorian Calendar (1582-10-15)
    if a >= 2299161
        α = fld(a - 1867216.25, 36524.25)
        a += 1 + α - fld(α, 4)
    end

    b = a + 1524
    c = fld(b - 122.1, 365.25)
    d = floor(365.25 * c)
    e = fld(b - d, 30.6001)

    day = b - d - floor(Integer, 30.6001 * e)

    if e < 14
        month = e - 1
    else
        month = e - 13
    end

    if month > 2
        year = c - 4716
    else
        year = c - 4715
    end

    hour, f = divrem(f * 24, 1)
    minute, f = divrem(f * 60, 1)
    second, f = divrem(f * 60, 1)
    millisecond = floor(f * 1e3)

    # Consider returning Date and Time objects for ns-precision
    return DateTime(year, month, day, hour, minute, second, millisecond)
end

function JDtoMJD(jd)
    jd - 2_400_000.5
end

function MJDtoJD(mjd)
    mjd + 2_400_000.5
end

function julianYears(jd, since = 2_451_545.0)
    (jd - since) / 36_525
end

function TAItoTT(tai::Dates.TimeType)
    tai + Dates.Millisecond(32_184)
end

function TTtoTAI(tdt::Dates.TimeType)
    tdt - Dates.Millisecond(32_184)
end

function UT1toUTC(ut1::Dates.TimeType, ΔUT1)
    ut1 - Dates.Microsecond(ΔUT1 * 1e6)
end

function UTCtoUT1(utc::Dates.TimeType, ΔUT1)
    utc + Dates.Microsecond(ΔUT1 * 1e6)
end

function UTCtoTAI(utc::Dates.TimeType, ΔAT)
    utc + Dates.Second(ΔAT)
end

function TAItoUTC(tai::Dates.TimeType, ΔAT)
    tai - Dates.Second(ΔAT)
end

# (3-54) - Vallado et al. 2013 page 193 (220)
# Astronomical Almanac 2012:B7
# ±30 μs for 1980-2050
function TT_to_TDB(tt)
    Ttt = julianYears(tt)

    # M = 357.527_723_3 + 35_999.050_34 * Ttt # [º]
    # Δλ_mean = 246.11 + 0.902_517_92 * (tt - 2_451_545.0) # [º]
    # tdb = tt + 0.001_657 * sin(M) + 0.000_022 * sin(Δλ_mean) # [s]

    M = 6.240_035_938_744_247 + 628.301_956_024_184_2 * Ttt # [rad]
    Δλ_mean = 4.2954298220832445 + 0.01575190926225078 * (tt - 2_451_545.0) # [rad]
    tdb = tt + 1.917824074074074e-8 * sin(M) + 2.546296296296296e-10 * sin(Δλ_mean) # [d]
end

# (3-53) - Vallado et al. 2013 page 193 (220)
# Fairhead et al. 1990
# ±10 μs for 1600-2200
function TT_to_TDB_fairhead_trunc(tt)
    Ttt = julianYears(tt)
    tdb = tt + (1657e-6 * sin( 628.3076 * Ttt + 6.2401) +
    22e-6 * sin( 575.3385 * Ttt + 4.2970) +
    14e-6 * sin(1256.6152 * Ttt + 6.1969) +
    5e-6 * sin( 606.9777 * Ttt + 4.0212) +
    5e-6 * sin(  52.9691 * Ttt + 0.4444) +
    2e-6 * sin(  21.3299 * Ttt + 5.5431) +
    10e-6 * Ttt * sin( 628.3076 * Ttt + 4.2490)) / 86400
end

# Fairhead et al. 1990
function TT_to_TDB_fairhead(tt)
    throw(error("Not implemented"))
end

# Harada et al. 2003
function TT_to_TDB_harada(tt)
    throw(error("Not implemented"))
end

function TDB_to_TCB(tdb, t0 = 2_443_144.500_372_5)
    tdb0 = -6.55e-5 / 86400
    return tdb + 1.550_519_767_72e-8 * (tdb - t0) - tdb0
end

function TT_to_TCB(tt, t0 = 2_443_144.500_372_5)
    return tt + 1.550_519_767_72e-8 * (tt - t0)
end

function TCG_to_TT(tcg, t0 = 2_443_144.500_372_5)
    L_G = 6.969_290_134e-10
    tcg - L_G * (tcg - t0)
end

function TT_to_TCG(tt, t0 = 2_443_144.500_372_5)
    L_G = 6.969_290_134e-10
    tt + L_G/(1 - L_G) * (tt - t0)
end

#=
t = DateTime(2004, 5, 14, 16, 43)
ut1 = UTCtoUT1(t, -0.463326)
tai = UTCtoTAI(t, 32)
gps = tai - Dates.Second(19)
tt = TAItoTT(tai)

tt_jd = dateToJD(tt)

tdb_jd = TT_to_TDB_fairhead_trunc(tt_jd)
julianYears(tdb_jd)
JDtoDate(tdb_jd)


tcb_jd = TDB_to_TCB(tdb_jd)
JDtoDate(tcb_jd)

tcb_jd = TT_to_TCB(tt_jd)
JDtoDate(tcb_jd)


tcg_jd = TT_to_TCG(tt_jd)
JDtoDate(tcg_jd)
=#
