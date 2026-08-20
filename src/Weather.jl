"Calendar-aligned weather and management inputs supplied from a unit-aware DataFrame."
@system Weather begin
    "Simulation calendar configured with the first weather timestamp."
    calendar(context) ~ ::Calendar

    "Unit-aware weather table aligned to calendar time."
    WDATA: weather_data ~ provide(index = :time, init = calendar.time, parameter)

    "Calendar date associated with the current weather record."
    DATE: calendar_date ~ drive::date(from = WDATA, by = :date)

    "Calendar day of year, numbered one on 1 January and 365 or 366 on the final day of the year."
    DOY(DATE): day_of_year => Cropbox.Dates.dayofyear(DATE) ~ track::int

    """
    Initial value of the monotonically increasing FST simulation day.

    The value is derived once from the calendar start date. It can therefore
    begin at any day of year without adding a separate time column to the
    weather table.
    """
    TIMEI(t0 = calendar.init): initial_source_simulation_day => begin
        Cropbox.Dates.dayofyear(Cropbox.Dates.Date(t0))
    end ~ preserve

    "Rate that advances FST simulation time by one day per elapsed day."
    RTIME: source_simulation_day_rate => 1 ~ preserve(u"d^-1")

    """
    Monotonically increasing FST simulation day used by source timing rules.

    Unlike calendar day of year, this accumulated state does not reset at a
    year boundary. Its increment follows the configured `Clock.step`, including
    subdaily simulations.
    """
    TIME(RTIME): source_simulation_day ~ accumulate(init = TIMEI)

    "Daily minimum air temperature."
    TMIN: minimum_temperature ~ drive(from = WDATA, by = :tmin, u"°C")

    "Daily maximum air temperature."
    TMAX: maximum_temperature ~ drive(from = WDATA, by = :tmax, u"°C")

    "Daily global solar radiation incident above the canopy."
    RAD: solar_radiation ~ drive(from = WDATA, by = :rad, u"MJ/m^2/d")

    "Daily rainfall depth."
    RAIN: rainfall ~ drive(from = WDATA, by = :rain, u"mm")

    "Daily fertilizer nitrogen input rate supplied at the management boundary."
    NFERT: fertilizer_nitrogen_input ~ drive(from = WDATA, by = :nfert, u"g/m^2/d")

    "Daily mean vapour pressure."
    VP: vapour_pressure ~ drive(from = WDATA, by = :vp, u"kPa")

    "Daily mean wind speed."
    WIND: wind_speed ~ drive(from = WDATA, by = :wind, u"m/s")

    "Daily mean air temperature calculated after scale-safe conversion to kelvin."
    TEMP(TMIN, TMAX): mean_temperature => (0.5 * (u"K"(TMIN) + u"K"(TMAX))) ~ track(u"°C")
end
