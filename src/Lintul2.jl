"""
Crop processes represented by the selected LINTUL-2 spring wheat program.

This system inherits the potential production processes of `Lintul1System` and
adds calendar and soil water controls on emergence, root zone expansion, a
single layer water balance, Penman evaporation and transpiration, drainage,
runoff, irrigation, and water stress effects on growth, allocation, and leaf
area development.
"""
@system Lintul2System(Lintul1System) begin
    "Calendar day of year supplied by the shared weather system."
    DOY  ~ hold

    "Monotonically accumulated simulation day corresponding to FST `TIME`."
    TIME ~ hold

    "Daily precipitation amount supplied by the shared weather system."
    RAIN ~ hold

    "Daily mean vapour pressure supplied by the shared weather system."
    VP   ~ hold

    "Daily mean wind speed supplied by the shared weather system."
    WIND ~ hold

    "Switch enabling automatic irrigation toward field capacity."
    IRRIGF: irrigation_switch                => 1     ~ preserve(parameter, min = 0, max = 1)

    "Volumetric soil water content at the air dryness limit."
    WCAD: air_dry_water_content              => 0.08  ~ preserve(parameter)

    "Volumetric soil water content at the permanent wilting point."
    WCWP: wilting_point_water_content        => 0.23  ~ preserve(parameter)

    "Volumetric soil water content at field capacity."
    WCFC: field_capacity_water_content       => 0.36  ~ preserve(parameter)

    "Volumetric water content above which oxygen stress begins."
    WCWET: anaerobic_threshold_water_content => 0.48  ~ preserve(parameter)

    "Volumetric soil water content at saturation."
    WCST: saturation_water_content           => 0.55  ~ preserve(parameter)

    "Initial volumetric water content in the rooted soil layer."
    WCI: initial_water_content               => 0.36  ~ preserve(parameter)

    "Initial depth of the rooted soil layer."
    ROOTDI: initial_rooting_depth            => 100   ~ preserve(parameter, u"mm")

    "Maximum depth available for root zone expansion."
    ROOTDM: maximum_rooting_depth            => 1200  ~ preserve(parameter, u"mm")

    "Maximum daily increase in rooting depth."
    RRDMAX: maximum_rooting_depth_rate       => 12    ~ preserve(parameter, u"mm/d")

    "Upper limit on daily drainage from the rooted soil layer."
    DRATE: maximum_drainage_rate             => 50    ~ preserve(parameter, u"mm/d")

    "Atmospheric demand constant used to derive the critical water content."
    TRANCO: transpiration_constant           => 8     ~ preserve(parameter, u"mm/d")

    "Latent heat required to evaporate one millimeter of water per square meter."
    LHVAP: latent_heat_of_vaporization       => 2.4   ~ preserve(parameter, u"MJ/m^2/mm")

    "Psychrometric constant used by the selected Penman calculation."
    PSYCH: psychrometric_constant            => 0.067 ~ preserve(parameter, u"kPa/K")

    "FST simulation day on which seedling emergence can begin."
    DOYEM: emergence_day                    => 60             ~ preserve(parameter)

    "Flag raised when monotonic FST simulation time reaches `DOYEM`."
    EDAY(TIME, DOYEM): emergence_day_reached => (TIME >= DOYEM) ~ flag

    "Whether leaf area has already been established."
    LAIEST(LAI): leaf_area_established      => (LAI > 0)      ~ flag

    """
    LINTUL-2 emergence switch used by thermal time and root extension.

    The selected source combines the calendar and water test with the
    established LAI test through `MAX`. The switch therefore remains active
    after leaf area has been established, even if water content later reaches
    the wilting point.
    """
    EMERG(EDAY, WC, WCWP, LAIEST): emergence_switch => begin
        (EDAY && WC > WCWP) || LAIEST
    end ~ flag

    "Thermal time rate gated by the selected LINTUL-2 emergence switch."
    RTSUM(DTEFF): thermal_time_rate => DTEFF ~ track(u"K", when = EMERG)

    "Flag indicating that root zone water content supports root extension."
    RWOK(WC, WCWP): root_zone_at_or_above_wilting_point => (WC >= WCWP)     ~ flag

    "Flag indicating that rooting depth remains below its maximum."
    RSPACE(ROOTD, ROOTDM): rootable_depth_available     => (ROOTD < ROOTDM) ~ flag

    """
    Condition for root zone deepening in the selected LINTUL-2 source.

    Rooting depth increases only before anthesis, while root zone water remains
    at or above wilting point, and while deeper soil is available.
    """
    RGROW(EMERG, RWOK, PRE_ANTHESIS, RSPACE): root_growth_condition => begin
        EMERG && RWOK && PRE_ANTHESIS && RSPACE
    end ~ flag

    "Rate of root depth extension while the source growth condition is active."
    RROOTD(RRDMAX): rooting_depth_growth_rate => RRDMAX ~ track(u"mm/d", when = RGROW)

    "Current rooting depth integrated only while `RGROW` is true."
    ROOTD(RROOTD): rooting_depth ~ accumulate(init = ROOTDI, u"mm")

    "Daily precipitation represented as the uniform rate used by the FST model."
    RAIN_rate(RAIN): rainfall_rate                          => (RAIN / 1u"d") ~ track(u"mm/d")

    "Maximum daily canopy interception proportional to leaf area index."
    RNINTC_MAX(LAI): maximum_canopy_interception_rate       => (0.25LAI)      ~ track(u"mm/d")

    """
    Rainfall intercepted by the canopy.

    The selected source limits daily interception to `0.25 * LAI` millimeters
    and cannot intercept more water than the incoming rainfall rate.
    """
    RNINTC(RAIN_rate, RNINTC_MAX): canopy_interception_rate => RAIN_rate      ~ track(u"mm/d", max = RNINTC_MAX)

    "Temperature interval above 0 degrees Celsius used in the Penman relations."
    TC(TEMP): celsius_temperature_interval          => (TEMP - 0u"°C")                                ~ track(u"K")

    "Absolute temperature reconstructed with the 273 K offset used by the source routine."
    TABS(TC): source_absolute_temperature           => (TC + 273u"K")                                  ~ track(u"K")

    "Daily blackbody radiation calculated from the source Stefan-Boltzmann relation."
    BBRAD(TABS): black_body_radiation               => (5.668e-8u"W/m^2/K^4" * TABS^4)                 ~ track(u"J/m^2/d")

    "Saturation vapour pressure at the daily mean temperature."
    SVP(TC): saturation_vapor_pressure              => (0.611u"kPa" * exp(17.4 * TC / (TC + 239u"K"))) ~ track(u"kPa")

    "Local slope of the saturation vapour pressure curve."
    SLOPE(SVP, TC): saturation_vapor_pressure_slope => (4158.6u"K" * SVP / (TC + 239u"K")^2)           ~ track(u"kPa/K")

    "Humidity-dependent factor applied to outgoing longwave radiation."
    RLWNF(VP, SVP): net_longwave_radiation_factor => (0.55 * (1 - VP / SVP)) ~ track(min = 0)

    "Daily outgoing net longwave radiation estimated by the source routine."
    RLWN(BBRAD, RLWNF): net_longwave_radiation    => (BBRAD * RLWNF)         ~ track(u"J/m^2/d")

    "Net radiation available to the soil surface after albedo and longwave losses."
    NRADS(RAD, RLWN): net_soil_radiation => (RAD * (1 - 0.15) - RLWN) ~ track(u"J/m^2/d")

    "Net radiation available to the crop canopy after albedo and longwave losses."
    NRADC(RAD, RLWN): net_crop_radiation => (RAD * (1 - 0.25) - RLWN) ~ track(u"J/m^2/d")

    "Wind speed retained under the FST Penman identifier `WN`."
    WN(WIND): source_wind_speed         => WIND                                      ~ track(u"m/s")

    "Aerodynamic vapor transfer coefficient derived from wind speed."
    WDF(WN): vapor_transfer_coefficient => (2.63u"mm/d/kPa" * (1 + 0.54u"s/m" * WN)) ~ track(u"mm/d/kPa")

    "Radiative contribution to the Penman soil evaporation demand."
    PENMRS(NRADS, SLOPE, PSYCH): soil_radiative_penman_term => (NRADS * SLOPE / (SLOPE + PSYCH)) ~ track(u"J/m^2/d")

    "Radiative contribution to the Penman canopy transpiration demand."
    PENMRC(NRADC, SLOPE, PSYCH): crop_radiative_penman_term => (NRADC * SLOPE / (SLOPE + PSYCH)) ~ track(u"J/m^2/d")

    "Aerodynamic contribution shared by the soil and canopy Penman calculations."
    PENMD(LHVAP, WDF, SVP, VP, PSYCH, SLOPE): aerodynamic_penman_term => begin
        LHVAP * WDF * (SVP - VP) * PSYCH / (SLOPE + PSYCH)
    end ~ track(u"J/m^2/d")

    """
    Potential soil evaporation from the selected LINTUL-2 Penman routine.

    The exposed soil fraction is `exp(-0.5 * LAI)`. Dividing the radiative and
    aerodynamic energy terms by `LHVAP` converts available energy to water
    depth without removing physical units.
    """
    PEVAP(LAI, PENMRS, PENMD, LHVAP): potential_evaporation => (exp(-0.5LAI) * (PENMRS + PENMD) / LHVAP) ~ track(u"mm/d")

    """
    Potential crop transpiration from the selected LINTUL-2 Penman routine.

    The canopy fraction is `1 - exp(-0.5 * LAI)`. The source subtracts one half
    of intercepted rainfall before constraining the result to nonnegative
    transpiration.
    """
    PTRAN(LAI, PENMRC, PENMD, LHVAP, RNINTC): potential_transpiration => begin
        (1 - exp(-0.5LAI)) * (PENMRC + PENMD) / LHVAP - 0.5RNINTC
    end ~ track(u"mm/d", min = 0)

    """
    Increment from wilting point to the critical water content threshold.

    The source makes the threshold depend on atmospheric demand through
    `PTRAN / (PTRAN + TRANCO)` and retains a minimum increment of 0.01.
    """
    WCCRI(PTRAN, TRANCO, WCFC, WCWP): critical_water_content_increment => begin
        PTRAN / (PTRAN + TRANCO) * (WCFC - WCWP)
    end ~ track(min = 0.01)

    "Critical root zone water content above which stress in the dry range is absent."
    WCCR(WCWP, WCCRI): critical_water_content => (WCWP + WCCRI) ~ track

    "Flag selecting the water stress relation for the wet range."
    WET(WC, WCCR): wet_stress_branch          => (WC > WCCR)    ~ flag

    "Stress factor for the wet range, reduced as water content approaches saturation."
    FRW(WCST, WC, WCWET): wet_side_water_stress_factor => begin
        (WCST - WC) / (WCST - WCWET)
    end ~ track(min = 0, max = 1, when = WET)

    "Stress factor for the dry range, reduced as water content approaches wilting point."
    FRD(WC, WCWP, WCCR): dry_side_water_stress_factor  => begin
        (WC - WCWP) / (WCCR - WCWP)
    end ~ track(min = 0, max = 1, when = !WET)

    """
    Water stress factor combining dry conditions and oxygen limitation under wet conditions.

    The mutually exclusive `FRD` and `FRW` branches give one under favorable
    water content and approach zero near wilting point or saturation.
    """
    FR(FRW, FRD): water_stress_factor         => (FRW + FRD)    ~ track

    "Reduction in potential soil evaporation as the surface dries below field capacity."
    EVAPF(WC, WCAD, WCFC): soil_evaporation_reduction_factor   => ((WC - WCAD) / (WCFC - WCAD)) ~ track(min = 0, max = 1)

    "Soil evaporation demand before the root zone storage limit is applied."
    EVAPP(PEVAP, EVAPF): evaporation_before_availability_limit => (PEVAP * EVAPF)               ~ track(u"mm/d")

    "Canopy transpiration demand after atmospheric water stress but before storage limitation."
    TRANP(PTRAN, FR): transpiration_before_availability_limit  => (PTRAN * FR)                  ~ track(u"mm/d")

    "Water amount remaining in the root zone at the air dryness limit."
    WAAD(WCAD, ROOTD): air_dry_water_amount                    => (WCAD * ROOTD)                ~ track(u"mm")

    "Combined evaporation and transpiration demand before the storage limit."
    EVTRAN(EVAPP, TRANP): combined_evapotranspiration_demand   => (EVAPP + TRANP)               ~ track(u"mm/d")

    "Flag indicating a nonzero combined evapotranspiration demand."
    EVNZ(EVTRAN): nonzero_evapotranspiration_demand            => (EVTRAN != 0u"mm/d")          ~ flag

    """
    Safe denominator corresponding to the source `NOTNUL(EVTRAN)` call.

    A value of one is used only when total demand is zero. The accompanying
    branch then contributes no evaporation or transpiration, avoiding division
    by zero without changing an active water flux.
    """
    EVNN(EVTRAN): source_notnul_evapotranspiration_demand      => EVTRAN                        ~ track(u"mm/d", init = 1u"mm/d", when = EVNZ)

    """
    Fraction of potential evapotranspiration supported by available soil water.

    The storage limit over one integration step prevents evaporation and
    transpiration from drawing the root zone below its air dryness limit.
    """
    AVAILF(WA, WAAD, DELT, EVNN): water_availability_factor => begin
        ((WA - WAAD) / DELT) / EVNN
    end ~ track(max = 1)

    "Actual soil evaporation after applying the root zone availability limit."
    EVAP(EVAPP, AVAILF): actual_evaporation                    => (EVAPP * AVAILF)      ~ track(u"mm/d")

    "Actual crop transpiration after applying the root zone availability limit."
    TRAN(TRANP, AVAILF): actual_transpiration                  => (TRANP * AVAILF)      ~ track(u"mm/d")

    "Flag indicating nonzero potential transpiration."
    PTNZ(PTRAN): nonzero_potential_transpiration               => (PTRAN != 0u"mm/d")   ~ flag

    """
    Safe denominator corresponding to `NOTNUL(PTRAN)` in the source.

    The fallback is active only when potential transpiration is zero, so the
    resulting reduction factor remains defined without altering nonzero demand.
    """
    PTNN(PTRAN): source_notnul_potential_transpiration         => PTRAN                 ~ track(u"mm/d", init = 1u"mm/d", when = PTNZ)

    "Ratio of actual to potential transpiration used to reduce crop growth."
    TRANRF(TRAN, PTNN): transpiration_reduction_factor => (TRAN / PTNN)       ~ track

    "Root zone water amount at field capacity."
    WAFC(WCFC, ROOTD): field_capacity_water_amount     => (WCFC * ROOTD)      ~ track(u"mm")

    "Root zone water amount at saturation."
    WAST(WCST, ROOTD): saturation_water_amount         => (WCST * ROOTD)      ~ track(u"mm")

    "Initial root zone water amount from initial content and depth."
    WAI(WCI, ROOTDI): initial_water_amount             => (WCI * ROOTDI)      ~ preserve(u"mm")

    "Net water input after interception, evaporation, and transpiration."
    NETRAIN(RAIN_rate, RNINTC, EVAP, TRAN): post_interception_water_input => (RAIN_rate - RNINTC - EVAP - TRAN) ~ track(u"mm/d")

    """
    Drainage that removes water above field capacity.

    The rate includes current net input and is limited by the source parameter
    `DRATE`.
    """
    DRAIN(WA, WAFC, NETRAIN, DRATE, DELT): drainage_rate => begin
        (WA - WAFC) / DELT + NETRAIN
    end ~ track(u"mm/d", min = 0, max = DRATE)

    "Runoff required to keep root zone storage at or below saturation."
    RUNOFF(WA, WAST, NETRAIN, DRAIN, DELT): runoff_rate => begin
        (WA - WAST) / DELT + NETRAIN - DRAIN
    end ~ track(u"mm/d", min = 0)

    "Automatic irrigation rate that refills the root zone to field capacity."
    IRRIG(IRRIGF, WAFC, WA, DELT, NETRAIN, DRAIN, RUNOFF): irrigation_rate => begin
        IRRIGF * ((WAFC - WA) / DELT - (NETRAIN - DRAIN - RUNOFF))
    end ~ track(u"mm/d", min = 0)

    "Rainfall plus irrigation remaining after canopy interception."
    INFLOW(RAIN_rate, IRRIG, RNINTC): net_surface_inflow => (RAIN_rate + IRRIG - RNINTC) ~ track(u"mm/d")

    "Water added to the modeled root zone as roots enter soil at field capacity."
    EXPLOR(RROOTD, WCFC): root_exploration_water_rate    => (RROOTD * WCFC)              ~ track(u"mm/d")

    """
    Net root zone water balance rate.

    Inputs are rainfall, irrigation, and water encountered by root extension.
    Interception, runoff, transpiration, evaporation, and drainage are losses.
    """
    RWA(RAIN_rate, EXPLOR, IRRIG, RNINTC, RUNOFF, TRAN, EVAP, DRAIN): water_accumulation_rate => begin
        RAIN_rate + EXPLOR + IRRIG - RNINTC - RUNOFF - TRAN - EVAP - DRAIN
    end ~ track(u"mm/d")

    "Root zone water amount integrated from the complete water balance."
    WA(RWA): root_zone_water_amount                    ~ accumulate(init = WAI, u"mm")

    "Volumetric root zone water content obtained from water amount and rooting depth."
    WC(WA, ROOTD): root_zone_water_content            => (WA / ROOTD) ~ track

    "Cumulative surface inflow after canopy interception."
    CINFLOW(INFLOW): cumulative_surface_inflow         ~ accumulate(u"mm")

    "Cumulative precipitation supplied to the model."
    TRAIN(RAIN_rate): cumulative_rainfall              ~ accumulate(u"mm")

    "Cumulative water encountered through root zone expansion."
    CEXPLOR(EXPLOR): cumulative_root_exploration_water ~ accumulate(u"mm")

    "Cumulative actual soil evaporation."
    CEVAP(EVAP): cumulative_evaporation                ~ accumulate(u"mm")

    "Cumulative actual crop transpiration."
    CTRAN(TRAN): cumulative_transpiration              ~ accumulate(u"mm")

    "Cumulative drainage from the rooted soil layer."
    CDRAIN(DRAIN): cumulative_drainage                 ~ accumulate(u"mm")

    "Cumulative surface runoff."
    CRUNOFF(RUNOFF): cumulative_runoff                 ~ accumulate(u"mm")

    """
    Residual of the integrated root zone water balance.

    Values near zero indicate closure among initial storage, cumulative inputs,
    final storage, and cumulative losses.
    """
    WBAL(WAI, CINFLOW, CEXPLOR, WA, CEVAP, CTRAN, CDRAIN, CRUNOFF): water_balance_error => begin
        WAI + CINFLOW + CEXPLOR - WA - CEVAP - CTRAN - CDRAIN - CRUNOFF
    end ~ track(u"mm")

    "Daily growth limited by water in the LINTUL radiation use efficiency relation."
    GTOTAL(PARINT, LUE, TRANRF): daily_growth_rate => (LUE * PARINT * TRANRF) ~ track(u"g/m^2/d")

    "Condition that establishes initial LAI only after the emergence day and with available soil water."
    SEEDLAI(EDAY, LAI, WC, WCWP): source_lai_seeding_condition => begin
        EDAY && LAI == 0 && WC > WCWP
    end ~ flag

    "Root allocation fraction before modification by water stress."
    FRTWET(FRTTB, TSUM_AFGEN): unstressed_root_fraction  => FRTTB(TSUM_AFGEN) ~ track

    "Leaf allocation fraction before modification by water stress."
    FLVT(FLVTB, TSUM_AFGEN): unstressed_leaf_fraction    => FLVTB(TSUM_AFGEN) ~ track

    "Stem allocation fraction before modification by water stress."
    FSTT(FSTTB, TSUM_AFGEN): unstressed_stem_fraction    => FSTTB(TSUM_AFGEN) ~ track

    "Storage organ allocation fraction before modification by water stress."
    FSOT(FSOTB, TSUM_AFGEN): unstressed_storage_fraction => FSOTB(TSUM_AFGEN) ~ track

    """
    Water stress modifier that preferentially increases root allocation.

    `1 / (TRANRF + 0.5)` is constrained to be at least one. The companion
    `FSHMOD` rescales shoot fractions so that total allocation remains closed.
    """
    FRTMOD(TRANRF): root_partitioning_modifier           => (1 / (TRANRF + 0.5))             ~ track(min = 1)

    "Root allocation fraction after applying the water stress modifier."
    FRT(FRTWET, FRTMOD): root_allocation_fraction        => (FRTWET * FRTMOD)                ~ track

    "Shoot rescaling factor that preserves total allocation after root modification."
    FSHMOD(FRT, FRTMOD): shoot_partitioning_modifier     => ((1 - FRT) / (1 - FRT / FRTMOD)) ~ track

    "Fraction of growth allocated to leaves after adjustment for water stress."
    FLV(FLVT, FSHMOD): leaf_allocation_fraction          => (FLVT * FSHMOD)                  ~ track

    "Fraction of growth allocated to stems after adjustment for water stress."
    FST(FSTT, FSHMOD): stem_allocation_fraction          => (FSTT * FSHMOD)                  ~ track

    "Fraction of growth allocated to storage organs after adjustment for water stress."
    FSO(FSOT, FSHMOD): storage_organ_allocation_fraction => (FSOT * FSHMOD)                  ~ track

    "Fraction of growth allocated to all aboveground organs."
    FSH(FRT): shoot_allocation_fraction                  => (1 - FRT)                        ~ track

    "Sum of organ allocation fractions after adjustment for water stress."
    FRACT(FRT, FLV, FST, FSO): allocation_sum            => (FRT + FLV + FST + FSO)          ~ track

    """
    Early exponential LAI growth reduced by transpiration stress.

    The source scales exponential canopy expansion by `TRANRF`, so water stress
    reduces early leaf area growth.
    """
    GLAE(LAI, RGRL, DTEFF, TRANRF, DELT): exponential_lai_growth_rate => begin
        LAI * expm1(RGRL * DTEFF * DELT) / DELT * TRANRF
    end ~ track(u"d^-1")
end

"""
LINTUL-2 water limited production model in Cropbox.

LINTUL (Light INTerception and UtiLization) was developed within the Wageningen
crop modeling tradition. LINTUL-2 extends the potential production model
LINTUL-1 with a simple root zone water balance for studying drought effects.
The selected FST (Fortran Simulation Translator) program is a spring wheat
example in which water availability affects transpiration, dry matter
production, allocation, leaf area development, and rooting depth.

# References

- Spitters, C. J. T. and Schapendonk, A. H. C. M. (1990). Evaluation of
  breeding strategies for drought tolerance in potato by means of crop growth
  simulation. *Plant and Soil*, 123, 193-203.
  <https://doi.org/10.1007/BF00011268>
- van Ittersum, M. K. et al. (2003). On approaches and applications of the
  Wageningen crop models. *European Journal of Agronomy*, 18, 201-234.
  <https://doi.org/10.1016/S1161-0301(02)00106-5>
"""
@system Lintul2Model(Weather, Lintul2System, Controller)
