using Test
using Cropbox
using LegacyCropModels

weather = load_weather_data(joinpath(@__DIR__, "..", "data", "weather.csv"); year = 2001)
base = weather_config(weather)

lintul2 = @config(
    base,
    Lintul2System => (
        :WCI => 0.36,
        :IRRIGF => 0,
    ),
)

lintul3 = @config(
    base,
    Lintul3System => (
        :WCI => 0.36,
        :IRRIGF => 1,
        :FERTNPI => 10u"g/m^2",
    ),
)

lintul3_saved_state = @config(
    base,
    Lintul3SavedStateSystem => (
        :WCI => 0.36,
        :IRRIGF => 1,
        :FERTNPI => 10u"g/m^2",
    ),
)

cases = (
    (
        name = "LINTUL 1",
        model = Lintul1Model,
        config = base,
        expected = 1530.72u"g/m^2",
    ),
    (
        name = "LINTUL 2",
        model = Lintul2Model,
        config = lintul2,
        expected = 875.49u"g/m^2",
    ),
    (
        name = "LINTUL 3",
        model = Lintul3Model,
        config = lintul3,
        expected = 1459.98u"g/m^2",
    ),
    (
        name = "LINTUL 3 saved state",
        model = Lintul3SavedStateModel,
        config = lintul3_saved_state,
        expected = 1472.96u"g/m^2",
    ),
)

@testset "Public LINTUL package" begin
    for case in cases
        @testset "$(case.name)" begin
            result = simulate(case.model;
                config = case.config,
                stop = :FINISH,
                verbose = false,
            )

            @test length(result.WTR) > 1
            @test isapprox(result.WTR[end], case.expected; atol = 0.01u"g/m^2")
        end
    end
end
