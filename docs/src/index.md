# LegacyCropModels.jl

[LegacyCropModels.jl](https://github.com/cropbox/LegacyCropModels.jl) provides
Cropbox implementations of selected LINTUL crop models. This guide shows how
to load the included daily weather data, configure each model, run a
simulation, and inspect the output.

## Installation

Install the package directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/cropbox/LegacyCropModels.jl")
```

## Getting Started

Load Cropbox and the model package:

```@example lintul
using Cropbox
using LegacyCropModels
```

The included weather file contains daily temperature, solar radiation,
rainfall, vapor pressure, and wind speed. Units are written in the column names
and are attached when the file is loaded. Because the example file contains
day of year rather than calendar dates, a year is supplied explicitly.

```@example lintul
weather_file = pkgdir(LegacyCropModels, "data", "weather.csv")
weather = load_weather_data(weather_file; year = 2001)
base = weather_config(weather)
```

`weather_config()` initializes the simulation calendar from the first weather
record and connects the table to the shared `Weather` system. Calendar day of
year is retained for seasonal calculations. The FST-style simulation day is
initialized from that start date and then accumulated with the configured
clock, so it remains monotonic across year boundaries without an extra input
column. A custom weather table can be used when it has the same required
columns and unit annotations.

## LINTUL 1

`Lintul1Model` represents potential production without water or nitrogen
limitation. The weather configuration is sufficient for a standard run.

```@example lintul
lintul1 = simulate(Lintul1Model;
    config = base,
    stop = :FINISH,
    verbose = false,
)

lintul1.WTR[end]
```

## LINTUL 2

`Lintul2Model` adds the root zone water balance. The following configuration
sets the initial volumetric water content and disables automatic irrigation.

```@example lintul
lintul2_config = @config(
    base,
    Lintul2System => (
        :WCI => 0.36,
        :IRRIGF => 0,
    ),
)

lintul2 = simulate(Lintul2Model;
    config = lintul2_config,
    stop = :FINISH,
    verbose = false,
)

lintul2.WTR[end]
```

## LINTUL 3

`Lintul3Model` adds soil and crop nitrogen dynamics. This example enables
automatic irrigation and supplies the initial inorganic soil nitrogen pool.
Daily fertilizer input is zero unless an `nfert` column is included in the
weather table.

```@example lintul
lintul3_config = @config(
    base,
    Lintul3System => (
        :WCI => 0.36,
        :IRRIGF => 1,
        :FERTNPI => 10u"g/m^2",
    ),
)

lintul3 = simulate(Lintul3Model;
    config = lintul3_config,
    stop = :FINISH,
    verbose = false,
)

lintul3.WTR[end]
```

The package also provides `Lintul3SavedStateModel`, which retains the saved
preanthesis development state of the selected FST source program. It can be run
with a configuration for `Lintul3SavedStateSystem`.

```@example lintul
saved_state_config = @config(
    base,
    Lintul3SavedStateSystem => (
        :WCI => 0.36,
        :IRRIGF => 1,
        :FERTNPI => 10u"g/m^2",
    ),
)

saved_state = simulate(Lintul3SavedStateModel;
    config = saved_state_config,
    stop = :FINISH,
    verbose = false,
)

saved_state.WTR[end]
```

## Working with Results

`simulate()` returns a data frame containing the model trajectory. Columns can
be selected for analysis or passed to Cropbox visualization functions.

```@example lintul
first(lintul3[:, [:time, :DVS, :LAI, :WTR]], 5)
```

Source variable names, descriptive aliases, units, and documentation can be
inspected with `look()`.

```julia
look(Lintul1System, :GLAE)
look(Lintul3System, :NNI)
```

See the [Cropbox documentation](https://cropbox.github.io/Cropbox.jl/stable/)
for additional configuration, simulation, visualization, and evaluation
workflows.
