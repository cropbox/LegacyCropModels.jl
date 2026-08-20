const REQUIRED_WEATHER_COLUMNS = (:doy, :tmin, :tmax, :rad, :rain, :vp, :wind)

function weather_day(value)
    day = Cropbox.deunitfy(value, u"d")
    isinteger(day) || throw(ArgumentError("weather `doy` values must be whole days"))
    Int(day)
end

function weather_year(value)
    year = Cropbox.deunitfy(value)
    isinteger(year) || throw(ArgumentError("weather `year` values must be integers"))
    Int(year)
end

"""
Prepare a daily weather table indexed by calendar time.

Column names containing units, such as `doy (d)` and `tmin (°C)`, are parsed
with `unitfy()`. Bare numeric columns are retained and interpreted in the units
declared by the corresponding `Weather` drives. A `date` column is used when
available. Otherwise, dates are constructed from `year` and `doy`, or from the
explicit `year` keyword for a yearless synthetic table. Dates must be
consecutive and their calendar day of year must agree with `doy`. Monotonic
FST simulation time is derived inside the `Weather` system from the calendar
start and the active clock step, so the input table needs no separate time
counter.
"""
function prepare_weather_data(
    weather::DataFrame;
    year::Union{Nothing,Integer} = nothing,
    timezone = tz"UTC",
)
    nrow(weather) > 0 || throw(ArgumentError("weather data are empty"))

    prepared = unitfy(copy(weather))
    rename!(prepared, lowercase.(names(prepared)))

    available = propertynames(prepared)
    missing_columns = filter(column -> column ∉ available, REQUIRED_WEATHER_COLUMNS)
    isempty(missing_columns) || throw(ArgumentError(
        "weather data are missing required columns: $(join(missing_columns, ", "))"
    ))

    :nfert in available || (prepared[!, :nfert] = zeros(nrow(prepared)))

    days = weather_day.(prepared.doy)
    dates = if :date in available
        Date.(prepared.date)
    elseif :year in available
        years = weather_year.(prepared.year)
        [Date(y, 1, 1) + Day(d - 1) for (y, d) in zip(years, days)]
    elseif !isnothing(year)
        [Date(year, 1, 1) + Day(d - 1) for d in days]
    else
        throw(ArgumentError(
            "weather data without a date or year column require the `year` keyword"
        ))
    end

    all(diff(dates) .== Day(1)) || throw(ArgumentError(
        "weather dates must be consecutive daily values"
    ))

    calendar_days = dayofyear.(dates)
    days == calendar_days || throw(ArgumentError(
        "weather `doy` values do not agree with their calendar dates"
    ))

    if :year in available
        weather_year.(prepared.year) == Dates.year.(dates) || throw(ArgumentError(
            "weather `year` values do not agree with their calendar dates"
        ))
    end

    prepared[!, :date] = dates
    prepared[!, :time] = [ZonedDateTime(DateTime(date), timezone) for date in dates]

    select!(prepared, :time, :date, Not([:time, :date]))
    prepared
end

"""
Expand daily weather to the interval used by a subdaily simulation.

Daily values are retained at every interval within the day. This preserves the
legacy interpretation of radiation and meteorological means and distributes a
daily rainfall amount uniformly through the model's daily rainfall rate.
"""
function resample_weather_data(weather::DataFrame, step)
    daily = prepare_weather_data(weather)
    interval = unitfy(step, u"hr")
    intervals_per_day = Cropbox.deunitfy(24u"hr" / interval)

    intervals_per_day > 0 || throw(ArgumentError("weather interval must be positive"))
    isinteger(intervals_per_day) || throw(ArgumentError(
        "weather interval must divide one day exactly"
    ))

    count = Int(round(intervals_per_day))
    count == 1 && return daily

    rows = repeat(collect(1:nrow(daily)), inner = count)
    expanded = daily[rows, :]
    period = convert(Dates.Millisecond, interval)
    expanded[!, :time] = [
        daily.time[i] + (j - 1) * period
        for i in 1:nrow(daily) for j in 1:count
    ]
    expanded[!, :date] = Date.(expanded.time)
    expanded
end

"Calendar start represented by the first timestamp of a prepared weather table."
weather_start(weather::DataFrame) = first(weather.time)

"""
Create the standard Cropbox clock, calendar, and weather configuration.

The calendar starts at the first timestamp of the prepared weather table. The
declared clock step must match the table interval checked by `provide`.
"""
function weather_config(weather::DataFrame; step = 1u"d")
    @config(
        Clock => :step => step,
        Calendar => :init => weather_start(weather),
        Weather => :weather_data => weather,
    )
end

"Load and prepare a daily weather CSV for `Weather`."
function load_weather_data(
    filepath::AbstractString;
    year::Union{Nothing,Integer} = nothing,
    timezone = tz"UTC",
)
    weather = CSV.read(filepath, DataFrame)
    prepare_weather_data(weather; year, timezone)
end
