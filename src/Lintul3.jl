"""
Crop processes represented by the selected LINTUL-3 spring wheat program.

This system retains the shared LINTUL crop structure while replacing the
phenology, allocation, water, initialization, and turnover relations specific
to the source. It represents development modified by daylength, root zone water,
staged soil evaporation, inorganic soil nitrogen, fertilizer recovery,
nitrogen demand and uptake, translocation, nitrogen stress, organ growth and
death, and carbon, water, and nitrogen balance diagnostics. Daily fertilizer
input enters through `NFERT`. The default development relation uses explicit
anchors of zero at emergence, one at anthesis, and two at maturity.
"""
@system Lintul3System(Lintul2System) begin
    """
    Daily fertilizer nitrogen rate supplied through the management boundary.

    This external input replaces the selected FST program's `FERTAB(TIME)`
    schedule, which is indexed by simulation time. Fertilizer recovery and soil
    nitrogen accounting remain model processes.
    """
    NFERT ~ hold

    "Source day of seedling emergence for the selected spring wheat program."
    DOYEM: emergence_day => 90 ~ preserve(parameter)

    """
    LINTUL-3 emergence switch recomputed at every integration step.

    Unlike LINTUL-2, the selected spring wheat source multiplies the calendar,
    water, and established LAI conditions. The switch can therefore return to
    false during drought after emergence.
    """
    EMERG(EDAY, WC, WCWP, LAIEST): emergence_switch => begin
        EDAY && WC > WCWP && LAIEST
    end ~ flag

    "Thermal time rate gated by the selected LINTUL-3 emergence switch."
    RTSUM(DTEFF): thermal_time_rate => DTEFF ~ track(u"K", when = EMERG)

    "Thermal time required for the preanthesis development component to reach one."
    TSUMAN: thermal_time_to_anthesis               => 800              ~ preserve(parameter, u"K*d")

    "Thermal time required for the postanthesis development component to reach one."
    TSUMMT: thermal_time_anthesis_to_maturity      => 1030             ~ preserve(parameter, u"K*d")

    "Thermal time threshold after which temperature-dependent leaf ageing begins."
    TSUMAG: thermal_time_to_leaf_ageing            => 800              ~ preserve(parameter, u"K*d")

    "Dry matter produced per unit of intercepted photosynthetically active radiation."
    LUE: light_use_efficiency                      => 2.8              ~ preserve(parameter, u"g/MJ")

    "Reference green leaf area per unit leaf dry mass before stage and nitrogen corrections."
    SLAC: reference_specific_leaf_area             => 0.022            ~ preserve(parameter, u"m^2/g")

    "Green leaf dry mass specified by the FST starting state."
    WLVGI: fst_initial_green_leaf_weight           => 2.4              ~ preserve(parameter, u"g/m^2")

    "Living root dry mass specified by the FST starting state."
    WRTLI: fst_initial_root_weight                 => 3.6              ~ preserve(parameter, u"g/m^2")

    "Effective thermal time at the source simulation boundary."
    TSUMI: initial_thermal_time                    => 0                ~ preserve(parameter, u"K*d")

    "Initial development stage derived from the preanthesis thermal-time scale."
    DVSI(TSUMI, TSUMAN): initial_development_stage => (TSUMI / TSUMAN) ~ preserve

    """
    Development stage multiplier applied to reference specific leaf area.

    The selected spring wheat parameterization keeps the multiplier at one,
    while retaining the table preserves the source interface and permits a
    stage-dependent replacement through configuration.
    """
    SLACF                                          => [
        0    1
        2    1
        2.1  1
    ] ~ interpolate(parameter)

    "Development-stage multiplier applied to reference specific leaf area at initialization."
    SLACFI(SLACF, DVSI): initial_sla_correction    => SLACF(DVSI)       ~ preserve

    "Specific leaf area at initialization after applying the development-stage multiplier."
    ISLA(SLAC, SLACFI): initial_specific_leaf_area => (SLAC * SLACFI)   ~ preserve(u"m^2/g")

    "Leaf area index at initialization, calculated from initial leaf mass and specific leaf area."
    LAII(WLVGI, ISLA): initial_leaf_area_index     => (WLVGI * ISLA)    ~ preserve

    "Green leaf dry mass used to initialize the accumulated living leaf pool."
    WLVI(WLVGI): initial_green_leaf_weight         => WLVGI             ~ preserve(u"g/m^2")

    "Root dry mass used to initialize the accumulated living root pool."
    WRTI(WRTLI): initial_root_weight               => WRTLI             ~ preserve(u"g/m^2")

    "Stem dry mass at the source simulation boundary."
    WSTI: initial_stem_weight                      => 0                 ~ preserve(parameter, u"g/m^2")

    "Storage organ dry mass at the source simulation boundary."
    WSOI: initial_storage_organ_weight             => 0                 ~ preserve(parameter, u"g/m^2")

    "Volumetric soil water content at the air dryness limit."
    WCAD: air_dry_water_content              => 0.1  ~ preserve(parameter)

    "Volumetric soil water content below which water uptake and root extension cease."
    WCWP: wilting_point_water_content        => 0.2  ~ preserve(parameter)

    "Volumetric soil water content used as the drainage and refill reference."
    WCFC: field_capacity_water_content       => 0.4  ~ preserve(parameter)

    "Volumetric water content above which oxygen stress begins."
    WCWET: anaerobic_threshold_water_content => 0.45 ~ preserve(parameter)

    "Volumetric soil water content at saturation, used to limit root-zone storage."
    WCST: saturation_water_content           => 0.5  ~ preserve(parameter)

    "Volumetric water content used to initialize root-zone water storage."
    WCI: initial_water_content               => 0.4  ~ preserve(parameter)

    "Rooting depth used to initialize the root-zone water balance."
    ROOTDI: initial_rooting_depth            => 100  ~ preserve(parameter, u"mm")

    "Maximum depth that the root zone can reach in the selected soil profile."
    ROOTDM: maximum_rooting_depth            => 1200 ~ preserve(parameter, u"mm")

    "Upper limit on daily root-zone depth extension."
    RRDMAX: maximum_rooting_depth_rate       => 12   ~ preserve(parameter, u"mm/d")

    "Upper limit on daily drainage below the root zone."
    DRATE: maximum_drainage_rate             => 30   ~ preserve(parameter, u"mm/d")

    "Source coefficient relating evaporative demand to water-stress onset."
    TRANCO: transpiration_constant           => 8    ~ preserve(parameter, u"mm/d")

    "Switch enabling automatic refill toward the source irrigation target."
    IRRIGF: irrigation_switch                => 1    ~ preserve(parameter, min = 0, max = 1)

    "Switch selecting flooded or nonflooded water-stress and irrigation relations."
    WMFAC: flooded_water_management_switch   => 0    ~ preserve(parameter)

    "Volumetric water content assigned to soil newly incorporated by root extension."
    WCSUBS: subsoil_water_content            => 0.3  ~ preserve(parameter)

    "Latitude used by the source astronomical daylength calculation."
    LAT: latitude       => 52        ~ preserve(parameter)

    "Value of pi retained from the FST source for astronomical calculations."
    PI: fst_pi_constant => 3.1415926 ~ preserve

    """
    Preanthesis development multiplier as a function of astronomical daylength.

    Values below eight hours reduce thermal time accumulation. The response is
    one from eight hours onward in the selected spring wheat parameterization.
    """
    PHOTTB                                                => [
         0  0
         8  1
        10  1
        12  1
        18  1
    ] ~ interpolate(parameter, knotunit = u"hr")

    "Maximum solar declination sine."
    SINDCM(PI): maximum_solar_declination_sine            => sin(PI * 23.45 / 180)                   ~ track

    "Solar declination sine."
    SINDEC(DOY, PI, SINDCM): solar_declination_sine       => (-SINDCM * cos(2PI * (DOY + 10) / 365)) ~ track

    "Solar declination cosine."
    COSDEC(SINDEC): solar_declination_cosine              => sqrt(1 - SINDEC^2)                      ~ track

    "Latitude sine."
    SINLAT(LAT, PI): latitude_sine                        => sin(PI * LAT / 180)                     ~ track

    "Latitude cosine."
    COSLAT(LAT, PI): latitude_cosine                      => cos(PI * LAT / 180)                     ~ track

    "Astronomical daylength numerator."
    A(SINLAT, SINDEC): astronomical_daylength_numerator   => (SINLAT * SINDEC)                       ~ track

    "Astronomical daylength denominator."
    B(COSLAT, COSDEC): astronomical_daylength_denominator => (COSLAT * COSDEC)                       ~ track

    "Ratio of the solar declination and latitude terms used in daylength."
    DAYLRATIO(A, B): astronomical_daylength_ratio         => (A / B)                                 ~ track

    """
    Astronomical day length used by the source photoperiod response.

    Solar declination and latitude are combined to estimate the interval
    between sunrise and sunset in hours.
    """
    DAYL(DAYLRATIO, PI): astronomical_daylength      => (12 * (1 + 2 / PI * asin(DAYLRATIO))) ~ track(u"hr")

    "Flag selecting extrapolation below the first photoperiod table knot."
    PHOT_LO(DAYL): daylength_below_photoperiod_table => (DAYL < 0u"hr")                       ~ flag

    "Flag selecting extrapolation above the last photoperiod table knot."
    PHOT_HI(DAYL): daylength_above_photoperiod_table => (DAYL > 18u"hr")                      ~ flag

    "Flag selecting interpolation within the photoperiod table domain."
    PHOT_IN(PHOT_LO, PHOT_HI): daylength_inside_photoperiod_table => begin
        !PHOT_LO && !PHOT_HI
    end ~ flag

    "Preanthesis condition that uses the tabulated photoperiod response."
    PHOT_TAB(PRE_ANTHESIS, PHOT_IN): tabulated_photoperiod_condition => begin
        PRE_ANTHESIS && PHOT_IN
    end ~ flag

    "Preanthesis condition that uses low-daylength extrapolation."
    PHOT_LOW(PRE_ANTHESIS, PHOT_LO): low_photoperiod_branch   => begin
        PRE_ANTHESIS && PHOT_LO
    end ~ flag

    "Preanthesis condition that uses the upper photoperiod endpoint."
    PHOT_HIGH(PRE_ANTHESIS, PHOT_HI): high_photoperiod_branch => begin
        PRE_ANTHESIS && PHOT_HI
    end ~ flag

    "Photoperiod multiplier interpolated within the tabulated daylength domain."
    PHOTPF_TABLE(PHOTTB, DAYL): interpolated_photoperiod_factor => PHOTTB(DAYL)    ~ track(when = PHOT_TAB)

    "Linear photoperiod multiplier below the first table knot at eight hours."
    PHOTPF_LOW(DAYL): low_daylength_photoperiod_extrapolation => (DAYL / 8u"hr") ~ track(when = PHOT_LOW)

    "Constant photoperiod multiplier above the last table knot."
    PHOTPF_HIGH: high_daylength_photoperiod_extrapolation     => 1               ~ track(when = PHOT_HIGH)

    "Preanthesis photoperiod multiplier assembled from the three table domains."
    PHOTPF1(PHOTPF_TABLE, PHOTPF_LOW, PHOTPF_HIGH): pre_anthesis_photoperiod_factor => begin
        PHOTPF_TABLE + PHOTPF_LOW + PHOTPF_HIGH
    end ~ track

    "Postanthesis multiplier of one, which removes photoperiod limitation."
    PHOTPF2: post_anthesis_photoperiod_factor                     => 1                   ~ track(when = !PRE_ANTHESIS)

    """
    Photoperiod multiplier applied to thermal time accumulation.

    Before anthesis it follows `PHOTTB(DAYL)` with extrapolation consistent
    with the source. After anthesis it is one.
    """
    PHOTPF(PHOTPF1, PHOTPF2): photoperiod_factor                  => (PHOTPF1 + PHOTPF2) ~ track

    "Effective thermal time rate after applying the photoperiod multiplier."
    RTSUMP(RTSUM, PHOTPF): photoperiod_adjusted_thermal_time_rate => (RTSUM * PHOTPF)    ~ track(u"K")

    "Thermal time adjusted for photoperiod and accumulated from emergence."
    TSUM(RTSUMP): thermal_time ~ accumulate(init = TSUMI, u"K*d")

    "Intercepted PAR set to zero before the source emergence day."
    PARINT(PAR, FINT): intercepted_par => (PAR * FINT) ~ track(u"MJ/m^2/d", when = EDAY)

    "Preanthesis development rate derived from photoperiod-adjusted thermal time."
    DVR1(RTSUMP, TSUMAN): pre_anthesis_development_rate  => (RTSUMP / TSUMAN) ~ track(u"d^-1", when = PRE_ANTHESIS)

    "Postanthesis development rate derived from photoperiod adjusted thermal time."
    DVR2(RTSUMP, TSUMMT): post_anthesis_development_rate => (RTSUMP / TSUMMT) ~ track(u"d^-1", when = !PRE_ANTHESIS)

    "Thermal time from emergence through anthesis to maturity."
    TTSUM(TSUMAN, TSUMMT): total_thermal_time         => (TSUMAN + TSUMMT)            ~ track(u"K*d")

    "Source stop condition reached when development or total thermal time passes maturity."
    FINISH(DVS, TSUM, TTSUM): source_finish_condition => (DVS > 2.01 || TSUM > TTSUM) ~ flag

    """
    Unstressed root allocation fraction over development stage.

    The values form the baseline partitioning schedule before water or nitrogen
    stress modifies the balance between root and shoot growth.
    """
    FRTTB => [
        0     0.6
        0.33  0.58
        0.4   0.55
        0.8   0.1
        1     0
        2     0
    ] ~ interpolate(parameter)

    """
    Unstressed green leaf allocation fraction over development stage.

    This schedule is combined with the stem and storage schedules after root
    allocation has been separated from total daily growth.
    """
    FLVTB => [
        0     0.4
        0.33  0.42
        0.4   0.405
        0.8   0.36
        1     0.1
        1.01  0
        2     0
    ] ~ interpolate(parameter)

    """
    Unstressed stem allocation fraction over development stage.

    The sharp change around anthesis transfers an increasing share of shoot
    growth from stems to storage organs.
    """
    FSTTB => [
        0     0
        0.33  0
        0.4   0.045
        0.8   0.54
        1     0.9
        1.01  0.25
        2     0
    ] ~ interpolate(parameter)

    """
    Unstressed storage organ allocation fraction over development stage.

    Allocation is zero before anthesis and becomes the dominant shoot sink
    during reproductive development in the selected spring wheat source.
    """
    FSOTB => [
        0     0
        0.33  0
        0.4   0
        0.8   0
        1     0
        1.01  0.75
        2     1
    ] ~ interpolate(parameter)

    """
    Temperature response of the relative leaf senescence rate.

    The interpolation knots define the daily fraction of green leaf area and
    biomass lost through temperature-dependent ageing. Separate declarations
    reproduce the source terminal-segment extrapolation outside the knot range.
    """
    RDRT  => [
        -10  0
         10  0.02
         15  0.03
         30  0.05
         50  0.09
    ] ~ interpolate(parameter, u"d^-1", knotunit = u"°C")

    "Leaf death rate below the first temperature knot using the source terminal slope."
    RDRTMP_LOW(TEMP): low_temperature_leaf_death_extrapolation   => begin
        0.001u"K^-1*d^-1" * (TEMP - (-10u"°C"))
    end ~ track(u"d^-1", when = RDRT_LO)

    "Leaf death rate above the final temperature knot using the source terminal slope."
    RDRTMP_HIGH(TEMP): high_temperature_leaf_death_extrapolation => begin
        0.09u"d^-1" + 0.002u"K^-1*d^-1" * (TEMP - 50u"°C")
    end ~ track(u"d^-1", when = RDRT_HI)

    """
    Root growth condition in the selected spring wheat source.

    Unlike inherited LINTUL-2 behavior, root extension is not restricted to
    the preanthesis period.
    """
    RGROW(EMERG, RSPACE): root_growth_condition => begin
        EMERG && RSPACE
    end ~ flag

    """
    Remaining rootable depth expressed as an equivalent rate.

    Dividing the remaining depth by `DELT` maps the source increment limit to
    the Cropbox rate state and prevents an integration step from exceeding
    `ROOTDM`.
    """
    RROOTD_LIMIT(ROOTDM, ROOTD, DELT): remaining_root_depth_rate => ((ROOTDM - ROOTD) / DELT) ~ track(u"mm/d")

    "Actual rooting depth increase, limited by both `RRDMAX` and remaining rootable depth."
    RROOTD(RRDMAX): rooting_depth_growth_rate                    => RRDMAX                    ~ track(u"mm/d", max = RROOTD_LIMIT, when = RGROW)

    "Current rooting depth integrated from the limited extension rate."
    ROOTD(RROOTD): rooting_depth ~ accumulate(init = ROOTDI, u"mm")

    "Water entering the modeled root zone as roots extend into subsoil."
    EXPLOR(RROOTD, WCSUBS): root_exploration_water_rate => (RROOTD * WCSUBS) ~ track(u"mm/d")

    "Flag indicating that rooting depth can be used as a denominator."
    RDNZ(ROOTD): nonzero_rooting_depth                  => (ROOTD != 0u"mm") ~ flag

    """
    Safe rooting-depth denominator corresponding to source `NOTNUL(ROOTD)`.

    The fallback is used only at zero depth. With a positive root zone, the
    actual rooting depth is retained for the water-content calculation.
    """
    RDNN(ROOTD): source_notnul_rooting_depth            => ROOTD             ~ track(u"mm", init = 1u"m", when = RDNZ)

    "Volumetric root zone water content calculated with a source compatible safe denominator."
    WC(WA, RDNN): root_zone_water_content               => (WA / RDNN)       ~ track

    "Zero interception capacity specified by the selected spring wheat source."
    RNINTC_MAX: maximum_canopy_interception_rate => 0 ~ preserve(u"mm/d")

    "Canopy interception fixed at zero for this source parameterization."
    RNINTC: canopy_interception_rate             => 0 ~ preserve(u"mm/d")

    "Potential soil evaporation from the source Penman terms and exposed soil fraction."
    PEVAP(LAI, PENMRS, PENMD, LHVAP): potential_evaporation   => (exp(-0.5LAI) * (PENMRS + PENMD) / LHVAP)       ~ track(u"mm/d", min = 0)

    "Potential crop transpiration from the source Penman terms and canopy fraction."
    PTRAN(LAI, PENMRC, PENMD, LHVAP): potential_transpiration => ((1 - exp(-0.5LAI)) * (PENMRC + PENMD) / LHVAP) ~ track(u"mm/d", min = 0)

    "Flag selecting flooded water management in the source `EVAPTR` routine."
    FLOODED(WMFAC): flooded_water_management => (WMFAC >= 1) ~ flag

    "Wet-range condition under flooded management."
    FWET(FLOODED, WET): flooded_above_critical_water     => begin
        FLOODED && WET
    end ~ flag

    "Wet-range condition under nonflooded management."
    NFWET(FLOODED, WET): nonflooded_above_critical_water => begin
        !FLOODED && WET
    end ~ flag

    "Water status factor fixed at one in the flooded wet branch."
    FRF: flooded_wet_water_stress_factor => 1 ~ track(when = FWET)

    "Oxygen stress factor that declines as nonflooded soil approaches saturation."
    FRW(WCST, WC, WCWET): wet_side_water_stress_factor => begin
        (WCST - WC) / (WCST - WCWET)
    end ~ track(min = 0, max = 1, when = NFWET)

    "Drought stress factor that declines as soil water approaches the wilting point."
    FRD(WC, WCWP, WCCR): dry_side_water_stress_factor  => begin
        (WC - WCWP) / (WCCR - WCWP)
    end ~ track(min = 0, max = 1, when = !WET)

    "Water status factor assembled from mutually exclusive flooded, wet, and dry branches."
    FR(FRF, FRW, FRD): water_stress_factor => (FRF + FRW + FRD) ~ track

    """
    Initial state for the source's staged soil evaporation routine.

    `DSLR` counts drying time after rainfall. Wet days reset the state, whereas
    consecutive dry days advance it by the active integration interval.
    """
    DSLRI: initial_dry_soil_stage                            => 0                           ~ preserve(parameter, u"d")

    "Flag raised when rainfall resets the staged soil evaporation calculation."
    RAINRST(RAIN_rate): evaporation_stage_reset              => (RAIN_rate >= 0.5u"mm/d")   ~ flag

    "Drying duration assigned on a rainfall reset interval."
    DSLR_RESET: reset_dry_soil_stage                         => 1u"d"                       ~ track(u"d", when = RAINRST)

    "Drying duration after one additional rain-free clock interval."
    DSLR_ADVANCE(DSLR, DELT): advanced_dry_soil_stage        => (DSLR + DELT)               ~ track(u"d", when = !RAINRST)

    "Next drying-duration state selected from reset or continued drying."
    DSLR_NEXT(DSLR_RESET, DSLR_ADVANCE): next_dry_soil_stage => (DSLR_RESET + DSLR_ADVANCE) ~ track(u"d")

    "Equivalent rate that moves the drying-duration state to `DSLR_NEXT`."
    RDSLR(DSLR, DSLR_NEXT, DELT): dry_soil_stage_change_rate => ((DSLR_NEXT - DSLR) / DELT) ~ track

    "Accumulated duration controlling second-stage soil evaporation."
    DSLR(RDSLR): dry_soil_evaporation_stage ~ accumulate(init = DSLRI, u"d")

    "Drying duration at the start of the active integration interval."
    DSLR0(DSLR_NEXT, DELT): previous_dry_soil_stage => (DSLR_NEXT - DELT) ~ track(u"d", min = 0u"d")

    """
    Decline in evaporation with the square root of drying duration.

    This is the increment over one integration step in the source `sqrt(DSLR)`
    relation, normalized by `DELT` so that it acts as a dimensionless rate
    modifier.
    """
    EVSMXI(DSLR_NEXT, DSLR0, DELT): dry_soil_evaporation_increment => begin
        (sqrt(DSLR_NEXT / 1u"d") - sqrt(DSLR0 / 1u"d")) / (DELT / 1u"d")
    end ~ track

    "Second-stage evaporation limit implied by the square-root drying relation."
    EVSMXT(PEVAP, EVSMXI): maximum_dry_soil_evaporation => (PEVAP * EVSMXI)     ~ track(u"mm/d")

    "Potential soil evaporation used on a rainfall reset interval."
    EVSW(PEVAP): wet_day_soil_evaporation               => PEVAP                ~ track(u"mm/d", when = RAINRST)

    "Soil evaporation on a dry interval, limited by potential evaporation."
    EVSD(EVSMXT, RAIN_rate): dry_day_soil_evaporation   => (EVSMXT + RAIN_rate) ~ track(u"mm/d", max = PEVAP, when = !RAINRST)

    "Soil evaporation demand assembled from wet- and dry-interval branches."
    EVS(EVSW, EVSD): staged_soil_evaporation            => (EVSW + EVSD)        ~ track(u"mm/d")

    "Staged soil evaporation demand before root-zone water availability is applied."
    EVAPP(EVS): evaporation_before_availability_limit   => EVS                  ~ track(u"mm/d")

    "Net water input after actual evaporation and transpiration, with interception omitted."
    NETRAIN(RAIN_rate, EVAP, TRAN): post_interception_water_input => (RAIN_rate - EVAP - TRAN) ~ track(u"mm/d")

    "Field capacity irrigation target."
    IRRFC(WAFC): field_capacity_irrigation_target                 => WAFC                      ~ track(u"mm", when = !FLOODED)

    "Saturation irrigation target."
    IRRST(WAST): saturation_irrigation_target                     => WAST                      ~ track(u"mm", when = FLOODED)

    "Root zone water amount targeted by automatic irrigation under the active water regime."
    IRRTGT(IRRFC, IRRST): source_irrigation_target => begin
        IRRFC + IRRST
    end ~ track(u"mm")

    """
    Automatic irrigation rate required to approach the active water target.

    The target is field capacity under nonflooded management and saturation
    under flooded management. Current net input, drainage, and runoff are
    included so the interval does not overfill the root zone.
    """
    IRRIG(IRRIGF, IRRTGT, WA, DELT, NETRAIN, DRAIN, RUNOFF): irrigation_rate => begin
        IRRIGF * ((IRRTGT - WA) / DELT - (NETRAIN - DRAIN - RUNOFF))
    end ~ track(u"mm/d", min = 0)

    "Surface water input to the root zone after combining rainfall and irrigation."
    INFLOW(RAIN_rate, IRRIG): net_surface_inflow => (RAIN_rate + IRRIG) ~ track(u"mm/d")

    "Net root zone water rate after rainfall, exploration, irrigation, runoff, transpiration, evaporation, and drainage."
    RWA(RAIN_rate, EXPLOR, IRRIG, RUNOFF, TRAN, EVAP, DRAIN): water_accumulation_rate => begin
        RAIN_rate + EXPLOR + IRRIG - RUNOFF - TRAN - EVAP - DRAIN
    end ~ track(u"mm/d")

    "Initial inorganic nitrogen available in the modeled soil pool."
    FERTNPI: initial_inorganic_soil_nitrogen       => 0      ~ preserve(parameter, u"g/m^2")

    "Potential daily addition of mineralized nitrogen to the inorganic soil pool."
    RTMINP: potential_soil_nitrogen_mineralization => 0.1    ~ preserve(parameter, u"g/m^2/d")

    "Fraction of applied fertilizer nitrogen recovered in the inorganic soil pool."
    NRF: fertilizer_n_recovery_fraction            => 0.7    ~ preserve(parameter, min = 0, max = 1)

    "Maximum nitrogen mass fraction used to calculate storage organ nitrogen demand."
    NMAXSO: maximum_storage_n_fraction             => 0.0165 ~ preserve(parameter)

    "Ratio used to derive maximum stem nitrogen concentration from the leaf value."
    LSNR: stem_to_leaf_maximum_n_ratio             => 0.5    ~ preserve(parameter)

    "Ratio used to derive maximum root nitrogen concentration from the leaf value."
    LRNR: root_to_leaf_maximum_n_ratio             => 0.5    ~ preserve(parameter)

    "Fraction of maximum organ nitrogen concentration treated as the optimum level."
    FRNX: optimum_to_maximum_n_ratio               => 0.5    ~ preserve(parameter)

    "Time coefficient controlling nitrogen transfer to storage organs."
    TCNT: nitrogen_translocation_time_coefficient  => 10     ~ preserve(parameter, u"d")

    "Development stage at which nitrogen translocation to storage organs begins."
    DVSNT: storage_n_translocation_stage           => 0.8    ~ preserve(parameter)

    "Development stage after which soil nitrogen uptake and mineralization stop."
    DVSNLT: soil_n_uptake_end_stage                => 1      ~ preserve(parameter)

    "Development stage at which root senescence begins."
    DVSDR: root_death_start_stage                  => 1      ~ preserve(parameter)

    "Minimum leaf nitrogen mass fraction retained after remobilization."
    RNFLV: residual_leaf_n_fraction                => 0.004  ~ preserve(parameter)

    "Minimum stem nitrogen mass fraction retained after remobilization."
    RNFST: residual_stem_n_fraction                => 0.002  ~ preserve(parameter)

    "Minimum root nitrogen mass fraction retained after remobilization."
    RNFRT: residual_root_n_fraction                => 0.002  ~ preserve(parameter)

    "Nitrogen mass fraction used to initialize the green leaf nitrogen pool."
    NFRLVI: initial_leaf_n_fraction                => 0.06   ~ preserve(parameter)

    "Nitrogen mass fraction used to initialize the stem nitrogen pool."
    NFRSTI: initial_stem_n_fraction                => 0.03   ~ preserve(parameter)

    "Nitrogen mass fraction used to initialize the living root nitrogen pool."
    NFRRTI: initial_root_n_fraction                => 0.03   ~ preserve(parameter)

    "Maximum root contribution relative to nitrogen available from leaves and stems."
    FNTRT: root_n_translocation_fraction           => 0.15   ~ preserve(parameter)

    "Response coefficient controlling how nitrogen stress reduces light use efficiency."
    NLUE: nitrogen_effect_on_lue                   => 0.2    ~ preserve(parameter)

    "Response coefficient controlling how nitrogen stress reduces juvenile leaf area expansion."
    NLAI: nitrogen_effect_on_juvenile_lai          => 1      ~ preserve(parameter)

    "Response coefficient controlling how nitrogen status modifies specific leaf area."
    NSLA: nitrogen_effect_on_specific_leaf_area    => 1      ~ preserve(parameter)

    "Response coefficient controlling how nitrogen stress shifts shoot partitioning."
    NPART: nitrogen_effect_on_partitioning         => 1      ~ preserve(parameter)

    "Maximum relative leaf death coefficient induced by nitrogen stress."
    RDRNS: relative_n_stress_leaf_death_rate       => 0.03   ~ preserve(parameter, u"d^-1")

    "Relative living root death rate after the configured development threshold."
    RDRRT: relative_root_death_rate                => 0.03   ~ preserve(parameter, u"d^-1")

    """
    Maximum green leaf nitrogen fraction over development stage.

    The declining schedule represents dilution and remobilization as the crop
    develops. Stem and root maxima are derived from this leaf reference.
    """
    NMXLV => [
        0    0.06
        0.4  0.04
        0.7  0.03
        1    0.02
        2    0.014
        2.1  0.014
    ] ~ interpolate(parameter)

    "Development stage bounded to the flat endpoint domain of nitrogen and SLA tables."
    DVSTAB(DVS): flat_endpoint_development_stage   => DVS             ~ track(min = 0, max = 2.1)

    "Maximum green leaf nitrogen mass fraction at the current development stage."
    NMAXLV(NMXLV, DVSTAB): maximum_leaf_n_fraction => NMXLV(DVSTAB)   ~ track

    "Maximum stem nitrogen mass fraction derived from the current leaf value."
    NMAXST(NMAXLV, LSNR): maximum_stem_n_fraction  => (LSNR * NMAXLV) ~ track

    "Maximum root nitrogen mass fraction derived from the current leaf value."
    NMAXRT(NMAXLV, LRNR): maximum_root_n_fraction  => (LRNR * NMAXLV) ~ track

    "Optimum leaf nitrogen mass fraction used in the nitrogen nutrition index."
    NOPTLV(NMAXLV, FRNX): optimum_leaf_n_fraction  => (FRNX * NMAXLV) ~ track

    "Optimum stem nitrogen mass fraction used in the nitrogen nutrition index."
    NOPTST(NMAXST, FRNX): optimum_stem_n_fraction  => (FRNX * NMAXST) ~ track

    "Combined green leaf and stem biomass used to calculate shoot nitrogen status."
    TBGMR(WLVG, WST): green_matter_biomass               => (WLVG + WST)         ~ track(u"g/m^2")

    "Flag indicating that green shoot biomass can be used as a denominator."
    TBGNZ(TBGMR): nonzero_green_matter_biomass           => (TBGMR != 0u"g/m^2") ~ flag

    "Safe green shoot biomass denominator corresponding to source `NOTNUL(TBGMR)`."
    TBGNN(TBGMR): source_notnul_green_matter_biomass     => TBGMR                ~ track(u"g/m^2", init = 1u"g/m^2", when = TBGNZ)

    "Nitrogen contained in green leaves and stems."
    NUPGMR(ANLV, ANST): green_matter_nitrogen            => (ANLV + ANST)        ~ track(u"g/m^2")

    "Actual nitrogen fraction of combined green leaf and stem biomass."
    NFGMR(NUPGMR, TBGNN): actual_green_matter_n_fraction => (NUPGMR / TBGNN)     ~ track

    "Residual nitrogen fraction weighted by current green leaf and stem biomass."
    NRMR(WLVG, WST, RNFLV, RNFST, TBGNN): residual_green_matter_n_fraction => begin
        (WLVG * RNFLV + WST * RNFST) / TBGNN
    end ~ track

    "Leaf nitrogen amount at the optimum concentration for current green leaf mass."
    NOPTL(WLVG, NOPTLV): optimum_leaf_n_amount                                 => (WLVG * NOPTLV)                  ~ track(u"g/m^2")

    "Stem nitrogen amount at the optimum concentration for current stem mass."
    NOPTS(WST, NOPTST): optimum_stem_n_amount                                  => (WST * NOPTST)                   ~ track(u"g/m^2")

    "Biomass weighted optimum nitrogen fraction of green leaves and stems."
    NOPTMR(NOPTL, NOPTS, TBGNN): optimum_green_matter_n_fraction               => ((NOPTL + NOPTS) / TBGNN)        ~ track

    "Difference between optimum and residual green shoot nitrogen fractions."
    NNIDEN(NOPTMR, NRMR): nitrogen_index_denominator                           => (NOPTMR - NRMR)                  ~ track

    "Flag indicating that the nitrogen nutrition index denominator is nonzero."
    NDENNZ(NNIDEN): nonzero_nitrogen_index_denominator                         => (NNIDEN != 0)                    ~ flag

    "Safe denominator corresponding to source `NOTNUL(NOPTMR - NRMR)`."
    NDENNN(NNIDEN): source_notnul_nitrogen_index_denominator                   => NNIDEN                           ~ track(init = 1, when = NDENNZ)

    """
    Unconstrained nitrogen nutrition index of green shoot biomass.

    The source relation is
    `(actual N fraction - residual N fraction) /
    (optimum N fraction - residual N fraction)`.
    """
    NNIRAW(NFGMR, NRMR, NDENNN): unconstrained_nitrogen_nutrition_index => ((NFGMR - NRMR) / NDENNN) ~ track

    "Nitrogen nutrition index constrained to the source range from 0.001 to 1 during an active emergence interval."
    NNI(NNIRAW): nitrogen_nutrition_index => NNIRAW ~ track(min = 0.001, max = 1, when = EMERG)

    "Multiplicative reduction of radiation use efficiency under nitrogen stress."
    NF_LUE(NNI, NLUE): nitrogen_lue_factor                                     => exp(-NLUE * (1 - NNI))           ~ track

    "Multiplicative reduction of juvenile LAI expansion under nitrogen stress."
    NF_LAI(NNI, NLAI): nitrogen_juvenile_lai_factor                            => exp(-NLAI * (1 - NNI))           ~ track

    "Development stage multiplier interpolated from the specific leaf area table."
    SLA_CORR(SLACF, DVSTAB): developmental_sla_correction => SLACF(DVSTAB) ~ track

    """
    Specific leaf area adjusted for development stage and nitrogen status.

    The source multiplies reference SLA by the developmental table and by
    `exp(-NSLA * (1 - NNI))`.
    """
    SLA(SLAC, SLA_CORR, NNI, NSLA): specific_leaf_area => begin
        SLAC * SLA_CORR * exp(-NSLA * (1 - NNI))
    end ~ track(u"m^2/g")

    "Nitrogen required to raise current green leaf biomass to its stage-dependent maximum concentration."
    NDEML(WLVG, NMAXLV, ANLV): leaf_n_demand               => (NMAXLV * WLVG - ANLV)  ~ track(u"g/m^2", min = 0)

    "Nitrogen required to raise current stem biomass to its stage-dependent maximum concentration."
    NDEMS(WST, NMAXST, ANST): stem_n_demand                => (NMAXST * WST - ANST)   ~ track(u"g/m^2", min = 0)

    "Nitrogen required to raise current root biomass to its stage-dependent maximum concentration."
    NDEMR(WRT, NMAXRT, ANRT): root_n_demand                => (NMAXRT * WRT - ANRT)   ~ track(u"g/m^2", min = 0)

    "Total nitrogen required to bring leaves, stems, and roots to their maximum concentrations."
    NDEMTO(NDEML, NDEMS, NDEMR): total_vegetative_n_demand => (NDEML + NDEMS + NDEMR) ~ track(u"g/m^2", min = 0)

    "Condition allowing soil nitrogen uptake before the configured development limit."
    NLIMIT(DVS, DVSNLT, WC, WCWP): soil_n_uptake_allowed => (DVS < DVSNLT && WC >= WCWP) ~ flag

    "Condition activating mineralization after emergence and before the uptake limit."
    NMINON(EMERG, NLIMIT): soil_n_mineralization_condition => begin
        EMERG && NLIMIT
    end ~ flag

    "Condition activating crop nitrogen uptake after the emergence day and before its stage limit."
    NUPON(EDAY, NLIMIT): soil_n_uptake_condition => begin
        EDAY && NLIMIT
    end ~ flag

    "Soil nitrogen mineralization rate while crop uptake is allowed."
    RTMIN(RTMINP): actual_soil_n_mineralization        => RTMINP         ~ track(u"g/m^2/d", when = NMINON)

    "Fertilizer amount entering during the active integration interval."
    FERTN(NFERT, DELT): fertilizer_n_amount            => (NFERT * DELT) ~ track(u"g/m^2")

    "Fertilizer nitrogen entering the inorganic soil pool after source recovery loss."
    FERTNS(FERTN, NRF): recovered_fertilizer_n_amount  => (FERTN * NRF)  ~ track(u"g/m^2")

    "Crop nitrogen uptake limited by both organ demand and inorganic soil supply."
    NUPSOIL(NDEMTO): soil_limited_crop_n_uptake_amount => NDEMTO         ~ track(u"g/m^2", max = TNSOIL)

    "Nonnegative nitrogen amount taken up from soil during the active interval."
    NUPTA(NUPSOIL): crop_n_uptake_amount               => NUPSOIL        ~ track(u"g/m^2", min = 0u"g/m^2")

    "Soil nitrogen uptake amount converted to a rate over the active clock interval."
    NUPTR(NUPTA, DELT): crop_n_uptake_rate             => (NUPTA / DELT) ~ track(u"g/m^2/d", when = NUPON)

    """
    Net rate of change in the inorganic soil nitrogen pool.

    Recovered fertilizer and mineralization add nitrogen. Crop uptake removes
    nitrogen during the active integration interval.
    """
    RNSOIL(FERTNS, DELT, NUPTR, RTMIN): inorganic_soil_n_change_rate => begin
        FERTNS / DELT - NUPTR + RTMIN
    end ~ track(u"g/m^2/d")

    """
    Inorganic nitrogen available to the plant in the represented soil layer.

    The state gains recovered fertilizer and mineralized nitrogen and loses
    crop uptake.
    """
    TNSOIL(RNSOIL): inorganic_soil_nitrogen ~ accumulate(init = FERTNPI, u"g/m^2")

    "Flag indicating that vegetative nitrogen demand can be used as a denominator."
    NDEMNZ(NDEMTO): nonzero_total_n_demand       => (NDEMTO != 0u"g/m^2") ~ flag

    "Safe total-demand denominator corresponding to source `NOTNUL(NDEMTO)`."
    NDEMNN(NDEMTO): source_notnul_total_n_demand => NDEMTO                ~ track(u"g/m^2", init = 1u"g/m^2", when = NDEMNZ)

    "Crop nitrogen uptake assigned to leaves in proportion to leaf demand."
    RNULV(NUPTR, NDEML, NDEMNN): leaf_n_uptake_rate => (NUPTR * NDEML / NDEMNN) ~ track(u"g/m^2/d", when = EMERG)

    "Crop nitrogen uptake assigned to stems in proportion to stem demand."
    RNUST(NUPTR, NDEMS, NDEMNN): stem_n_uptake_rate => (NUPTR * NDEMS / NDEMNN) ~ track(u"g/m^2/d", when = EMERG)

    "Crop nitrogen uptake assigned to roots in proportion to root demand."
    RNURT(NUPTR, NDEMR, NDEMNN): root_n_uptake_rate => (NUPTR * NDEMR / NDEMNN) ~ track(u"g/m^2/d", when = EMERG)

    "Leaf nitrogen available for translocation above the residual leaf concentration."
    ATNLV(ANLV, WLVG, RNFLV): translocatable_leaf_n                    => (ANLV - WLVG * RNFLV)          ~ track(u"g/m^2", min = 0u"g/m^2")

    "Stem nitrogen available for translocation above the residual stem concentration."
    ATNST(ANST, WST, RNFST): translocatable_stem_n                     => (ANST - WST * RNFST)           ~ track(u"g/m^2", min = 0u"g/m^2")

    "Root nitrogen available above the residual root concentration."
    ANRTAV(ANRT, WRT, RNFRT): translocatable_root_n_available          => (ANRT - WRT * RNFRT)           ~ track(u"g/m^2")

    "Root contribution to translocatable nitrogen, limited by the available root pool."
    ATNRT(ATNLV, ATNST, FNTRT): translocatable_root_n                  => ((ATNLV + ATNST) * FNTRT)      ~ track(u"g/m^2", max = ANRTAV)

    "Total nitrogen available for redistribution from leaves, stems, and roots."
    ATN(ATNLV, ATNST, ATNRT): total_translocatable_n                   => (ATNLV + ATNST + ATNRT)        ~ track(u"g/m^2")

    """
    Storage organ nitrogen demand represented by a proportional filling rate.

    The deficit from `NMAXSO * WSO` is distributed over the source
    translocation time coefficient `TCNT`.
    """
    NDEMSO(WSO, ANSO, NMAXSO, TCNT): storage_n_demand_rate             => ((WSO * NMAXSO - ANSO) / TCNT) ~ track(u"g/m^2/d", min = 0)

    "Condition enabling nitrogen transfer to storage organs after `DVSNT`."
    NSTORE(DVS, DVSNT): storage_n_translocation_condition              => (DVS >= DVSNT)                 ~ flag

    "Potential nitrogen supply rate from the total translocatable pool."
    NSUPSO(ATN, TCNT): storage_n_supply_rate                           => (ATN / TCNT)                   ~ track(u"g/m^2/d", when = NSTORE)

    "Storage organ nitrogen accumulation limited by demand and translocatable supply."
    RNSO(NDEMSO): storage_n_accumulation_rate                          => NDEMSO                         ~ track(u"g/m^2/d", max = NSUPSO)

    "Flag indicating that total translocatable nitrogen can be used as a denominator."
    ATNNZ(ATN): nonzero_translocatable_n                               => (ATN != 0u"g/m^2")             ~ flag

    "Safe translocatable-nitrogen denominator corresponding to source `NOTNUL(ATN)`."
    ATNNN(ATN): source_notnul_translocatable_n                         => ATN                            ~ track(u"g/m^2", init = 1u"g/m^2", when = ATNNZ)

    "Leaf contribution to storage-organ nitrogen transfer, proportional to its available pool."
    RNTLV(RNSO, ATNLV, ATNNN): leaf_n_translocation_rate               => (RNSO * ATNLV / ATNNN)         ~ track(u"g/m^2/d")

    "Stem contribution to storage-organ nitrogen transfer, proportional to its available pool."
    RNTST(RNSO, ATNST, ATNNN): stem_n_translocation_rate               => (RNSO * ATNST / ATNNN)         ~ track(u"g/m^2/d")

    "Root contribution to storage-organ nitrogen transfer, proportional to its available pool."
    RNTRT(RNSO, ATNRT, ATNNN): root_n_translocation_rate               => (RNSO * ATNRT / ATNNN)         ~ track(u"g/m^2/d")

    """
    Growth route selected by the more limiting resource.

    Water limitation scales the RUE relation by `TRANRF`. Otherwise the source
    applies the exponential nitrogen response `NF_LUE`.
    """
    WGROW(TRANRF, NNI): water_limited_growth_route  => (TRANRF <= NNI) ~ flag

    "Condition selecting nitrogen-stress rather than water-stress partitioning."
    NPSEL(TRANRF, NNI): nitrogen_partitioning_route => (TRANRF >= NNI) ~ flag

    "Daily dry matter production when water is the dominant growth limitation."
    GTOTALW(PARINT, LUE, TRANRF): water_limited_growth_rate     => (LUE * PARINT * TRANRF) ~ track(u"g/m^2/d", when = WGROW)

    "Daily dry matter production when nitrogen is the dominant growth limitation."
    GTOTALN(PARINT, LUE, NF_LUE): nitrogen_limited_growth_rate  => (LUE * PARINT * NF_LUE) ~ track(u"g/m^2/d", when = !WGROW)

    "Daily dry matter production assembled from the active resource-limitation route."
    GTOTAL(GTOTALW, GTOTALN): daily_growth_rate                 => (GTOTALW + GTOTALN)     ~ track(u"g/m^2/d")

    "Flag selecting linear extrapolation below the first allocation-table knot."
    PART_LO(DVS): dvs_below_partitioning_tables => (DVS < 0) ~ flag

    "Flag selecting linear extrapolation above the last allocation-table knot."
    PART_HI(DVS): dvs_above_partitioning_tables => (DVS > 2) ~ flag

    "Flag selecting interpolation within the allocation-table domain."
    PART_IN(PART_LO, PART_HI): dvs_inside_partitioning_tables => begin
        !PART_LO && !PART_HI
    end ~ flag

    "Unstressed root allocation interpolated within the tabulated development domain."
    FRTWET_TABLE(FRTTB, DVS): interpolated_unstressed_root_fraction => FRTTB(DVS)             ~ track(when = PART_IN)

    "Linear root-allocation extrapolation below the first development-stage knot."
    FRTWET_LOW(DVS): low_dvs_root_fraction_extrapolation     => (0.6 - 0.02DVS / 0.33) ~ track(when = PART_LO)

    "Root allocation beyond the final development-stage knot."
    FRTWET_HIGH: high_dvs_root_fraction_extrapolation        => 0                      ~ track(when = PART_HI)

    "Unstressed root allocation assembled across interpolation and endpoint domains."
    FRTWET(FRTWET_TABLE, FRTWET_LOW, FRTWET_HIGH): unstressed_root_fraction => begin
        FRTWET_TABLE + FRTWET_LOW + FRTWET_HIGH
    end ~ track

    "Unstressed leaf allocation interpolated within the tabulated development domain."
    FLVT_TABLE(FLVTB, DVS): interpolated_unstressed_leaf_fraction => FLVTB(DVS)             ~ track(when = PART_IN)

    "Linear leaf-allocation extrapolation below the first development-stage knot."
    FLVT_LOW(DVS): low_dvs_leaf_fraction_extrapolation     => (0.4 + 0.02DVS / 0.33) ~ track(when = PART_LO)

    "Leaf allocation beyond the final development-stage knot."
    FLVT_HIGH: high_dvs_leaf_fraction_extrapolation        => 0                      ~ track(when = PART_HI)

    "Unstressed leaf allocation assembled across interpolation and endpoint domains."
    FLVT(FLVT_TABLE, FLVT_LOW, FLVT_HIGH): unstressed_leaf_fraction => begin
        FLVT_TABLE + FLVT_LOW + FLVT_HIGH
    end ~ track

    "Unstressed stem allocation interpolated within the tabulated development domain."
    FSTT_TABLE(FSTTB, DVS): interpolated_unstressed_stem_fraction => FSTTB(DVS)              ~ track(when = PART_IN)

    "Stem allocation below the first development-stage knot."
    FSTT_LOW: low_dvs_stem_fraction_extrapolation          => 0                       ~ track(when = PART_LO)

    "Linear stem-allocation extrapolation above the final development-stage knot."
    FSTT_HIGH(DVS): high_dvs_stem_fraction_extrapolation   => (-0.25(DVS - 2) / 0.99) ~ track(when = PART_HI)

    "Unstressed stem allocation assembled across interpolation and endpoint domains."
    FSTT(FSTT_TABLE, FSTT_LOW, FSTT_HIGH): unstressed_stem_fraction => begin
        FSTT_TABLE + FSTT_LOW + FSTT_HIGH
    end ~ track

    "Unstressed storage allocation interpolated within the tabulated development domain."
    FSOT_TABLE(FSOTB, DVS): interpolated_unstressed_storage_fraction => FSOTB(DVS)                 ~ track(when = PART_IN)

    "Storage allocation below the first development-stage knot."
    FSOT_LOW: low_dvs_storage_fraction_extrapolation          => 0                          ~ track(when = PART_LO)

    "Linear storage-allocation extrapolation above the final development-stage knot."
    FSOT_HIGH(DVS): high_dvs_storage_fraction_extrapolation   => (1 + 0.25(DVS - 2) / 0.99) ~ track(when = PART_HI)

    "Unstressed storage allocation assembled across interpolation and endpoint domains."
    FSOT(FSOT_TABLE, FSOT_LOW, FSOT_HIGH): unstressed_storage_fraction => begin
        FSOT_TABLE + FSOT_LOW + FSOT_HIGH
    end ~ track

    """
    Water stress modifier that increases root allocation when water is more
    limiting than nitrogen.

    Allocation under nitrogen limitation follows the separate `FLVMOD` and `MODIF`
    branch below.
    """
    FRTMOD(TRANRF): root_partitioning_modifier          => (1 / (TRANRF + 0.5))               ~ track(min = 1, when = !NPSEL)

    "Root allocation fraction after the water stress modifier is applied."
    FRTW(FRTWET, FRTMOD): water_stressed_root_fraction  => (FRTWET * FRTMOD)                  ~ track(when = !NPSEL)

    "Shoot rescaling factor that preserves allocation closure after root adjustment."
    FSHMOD(FRTW, FRTMOD): shoot_partitioning_modifier   => ((1 - FRTW) / (1 - FRTW / FRTMOD)) ~ track(when = !NPSEL)

    "Nitrogen stress modifier that directly reduces the unstressed leaf fraction."
    FLVMOD(NNI, NPART): nitrogen_leaf_partitioning_modifier => begin
        exp(-NPART * (1 - NNI))
    end ~ track(when = NPSEL)

    "Leaf allocation fraction after the nitrogen stress modifier is applied."
    FLVN(FLVT, FLVMOD): nitrogen_limited_leaf_fraction => begin
        FLVT * FLVMOD
    end ~ track(when = NPSEL)

    "Modifier that rescales nonleaf shoot allocation after the leaf adjustment."
    MODIF(FLVN, FLVMOD): nitrogen_nonleaf_partitioning_modifier => begin
        (1 - FLVN) / (1 - FLVN / FLVMOD)
    end ~ track(when = NPSEL)

    "Root allocation contribution selected when water is the dominant limitation."
    FRTWATER(FRTW): water_limited_root_fraction            => FRTW                    ~ track(when = !NPSEL)

    "Root allocation contribution selected when nitrogen is the dominant limitation."
    FRTN(FRTWET, MODIF): nitrogen_limited_root_fraction    => (FRTWET * MODIF)        ~ track(when = NPSEL)

    "Active root allocation fraction assembled from the mutually exclusive stress routes."
    FRT(FRTWATER, FRTN): root_allocation_fraction          => (FRTWATER + FRTN)       ~ track

    "Leaf allocation contribution selected when water is the dominant limitation."
    FLVWATER(FLVT, FSHMOD): water_limited_leaf_fraction    => (FLVT * FSHMOD)         ~ track(when = !NPSEL)

    "Active leaf allocation fraction assembled from the mutually exclusive stress routes."
    FLV(FLVWATER, FLVN): leaf_allocation_fraction          => (FLVWATER + FLVN)       ~ track

    "Stem allocation contribution selected when water is the dominant limitation."
    FSTWATER(FSTT, FSHMOD): water_limited_stem_fraction    => (FSTT * FSHMOD)         ~ track(when = !NPSEL)

    "Stem allocation contribution selected when nitrogen is the dominant limitation."
    FSTN(FSTT, MODIF): nitrogen_limited_stem_fraction      => (FSTT * MODIF)          ~ track(when = NPSEL)

    "Active stem allocation fraction assembled from the mutually exclusive stress routes."
    FST(FSTWATER, FSTN): stem_allocation_fraction          => (FSTWATER + FSTN)       ~ track

    "Storage allocation contribution selected when water is the dominant limitation."
    FSOWATER(FSOT, FSHMOD): water_limited_storage_fraction => (FSOT * FSHMOD)         ~ track(when = !NPSEL)

    "Storage allocation contribution selected when nitrogen is the dominant limitation."
    FSON(FSOT, MODIF): nitrogen_limited_storage_fraction   => (FSOT * MODIF)          ~ track(when = NPSEL)

    "Active storage organ allocation fraction assembled from the stress routes."
    FSO(FSOWATER, FSON): storage_organ_allocation_fraction => (FSOWATER + FSON)       ~ track

    "Fraction of total daily growth allocated above ground after root allocation."
    FSH(FRT): shoot_allocation_fraction                    => (1 - FRT)               ~ track

    "Sum of active root, leaf, stem, and storage fractions used to check allocation closure."
    FRACT(FRT, FLV, FST, FSO): allocation_sum              => (FRT + FLV + FST + FSO) ~ track

    """
    Juvenile exponential LAI growth under combined water and nitrogen stress.

    The source scales exponential canopy expansion by both `TRANRF` and
    `NF_LAI`, so water and nitrogen stress jointly reduce early leaf area growth.
    """
    GLAE(LAI, RGRL, DTEFF, TRANRF, NF_LAI, DELT): exponential_lai_growth_rate => begin
        LAI * expm1(RGRL * DTEFF * DELT) / DELT * TRANRF * NF_LAI
    end ~ track(u"d^-1")

    "Leaf area growth obtained from new green leaf mass and current specific leaf area."
    GLAL(GLV, SLA): linear_lai_growth_rate => (GLV * SLA) ~ track(u"d^-1")

    "Condition selecting exponential canopy expansion at early development and low LAI."
    JUVLAI(DVS, LAI, LAIJ): juvenile_lai_condition => begin
        DVS < 0.2 && LAI < LAIJ
    end ~ flag

    "Exponential LAI growth contribution active during the juvenile canopy phase."
    GLAIJ(GLAE): juvenile_lai_growth_rate     => GLAE            ~ track(u"d^-1", when = JUVLAI)

    "Leaf mass based LAI growth contribution active after the juvenile canopy phase."
    GLAIM(GLAL): mature_lai_growth_rate       => GLAL            ~ track(u"d^-1", when = !JUVLAI)

    "Gross LAI increase assembled from the mutually exclusive canopy growth phases."
    GLAI(GLAIJ, GLAIM): gross_lai_growth_rate => (GLAIJ + GLAIM) ~ track(u"d^-1", when = EDAY)

    "Green leaf area per unit ground area integrated from growth and death rates."
    LAI(RLAI): leaf_area_index ~ accumulate(init = LAII)

    "Condition activating temperature dependent leaf ageing after `TSUMAG`."
    LAGE(TSUM, TSUMAG): leaf_ageing_condition => (TSUM >= TSUMAG) ~ flag

    "Relative leaf death rate from temperature dependent ageing."
    RDRDV(RDRTMP): aging_leaf_death_rate      => RDRTMP           ~ track(u"d^-1", when = LAGE)

    "Relative leaf death rate caused by canopy shading above `LAICR`."
    RDRSH(LAI, LAICR, RDRSHM): shading_leaf_death_rate => begin
        RDRSHM * (LAI - LAICR) / LAICR
    end ~ track(u"d^-1", min = 0)

    """
    Additional leaf death caused by nitrogen stress.

    The source scales green leaf biomass by `RDRNS * (1 - NNI)`, so this loss
    vanishes at optimal nitrogen status.
    """
    DLVNS(WLVG, RDRNS, NNI): nitrogen_stress_leaf_biomass_death_rate => begin
        WLVG * RDRNS * (1 - NNI)
    end ~ track(u"g/m^2/d", when = EDAY)

    "Leaf area loss associated with nitrogen-stress leaf biomass death."
    DLAINS(DLVNS, SLA): nitrogen_stress_lai_death_rate  => (DLVNS * SLA)    ~ track(u"d^-1")

    "Green leaf biomass lost through temperature- or shading-driven senescence."
    DLVS(WLVG, RDR): senescence_leaf_biomass_death_rate => (WLVG * RDR)     ~ track(u"g/m^2/d")

    "Leaf area lost through temperature- or shading-driven senescence."
    DLAIS(LAI, RDR): senescence_lai_death_rate          => (LAI * RDR)      ~ track(u"d^-1")

    "Total green leaf biomass loss from senescence and nitrogen stress."
    DLV(DLVS, DLVNS): leaf_biomass_death_rate           => (DLVS + DLVNS)   ~ track(u"g/m^2/d")

    "Total LAI loss from senescence and nitrogen stress."
    DLAI(DLAIS, DLAINS): leaf_area_death_rate           => (DLAIS + DLAINS) ~ track(u"d^-1")

    "Condition activating living root senescence at and after `DVSDR`."
    RTDEATH(DVS, DVSDR): root_death_condition => (DVS >= DVSDR) ~ flag

    "Living root biomass loss after the root senescence threshold."
    DRRT(WRT, RDRRT): root_death_rate         => (WRT * RDRRT)  ~ track(u"g/m^2/d", when = RTDEATH)

    "Net green leaf biomass rate after subtracting senescence and nitrogen-stress losses."
    RWLVG(GLV, DLV): net_green_leaf_growth_rate        => (GLV - DLV)           ~ track(u"g/m^2/d", when = EMERG)

    "Net living root growth after allocation and root death."
    RWRT(GTOTAL, FRT, DRRT): root_growth_rate          => (GTOTAL * FRT - DRRT) ~ track(u"g/m^2/d", when = EMERG)

    "Daily stem biomass growth allocated from total dry matter production."
    RWST(GTOTAL, FST): stem_growth_rate                => (GTOTAL * FST)        ~ track(u"g/m^2/d", when = EMERG)

    "Daily storage-organ biomass growth allocated from total dry matter production."
    RWSO(GTOTAL, FSO): storage_organ_growth_rate       => (GTOTAL * FSO)        ~ track(u"g/m^2/d", when = EMERG)

    "Cumulative root biomass removed from the living root pool."
    WDRT(DRRT): dead_root_biomass                  ~ accumulate(u"g/m^2")

    "Nitrogen removed with dead leaf biomass at the residual leaf concentration."
    RNLDLV(DLV, RNFLV): leaf_n_loss_rate  => (DLV * RNFLV)  ~ track(u"g/m^2/d")

    "Nitrogen removed with dead root biomass at the residual root concentration."
    RNLDRT(DRRT, RNFRT): root_n_loss_rate => (DRRT * RNFRT) ~ track(u"g/m^2/d")

    "Initial nitrogen in green leaves from source leaf mass and concentration."
    ANLVI(WLVGI, NFRLVI): initial_leaf_nitrogen    => (WLVGI * NFRLVI)         ~ preserve(u"g/m^2")

    "Initial nitrogen in stems from source stem mass and concentration."
    ANSTI(WSTI, NFRSTI): initial_stem_nitrogen     => (WSTI * NFRSTI)          ~ preserve(u"g/m^2")

    "Initial nitrogen in living roots from source root mass and concentration."
    ANRTI(WRTLI, NFRRTI): initial_root_nitrogen    => (WRTLI * NFRRTI)         ~ preserve(u"g/m^2")

    "Initial nitrogen amount in storage organs at the source boundary."
    ANSOI: initial_storage_nitrogen                => 0                        ~ preserve(u"g/m^2")

    "Net leaf nitrogen rate after uptake, translocation, and senescence loss."
    RNLV(RNULV, RNTLV, RNLDLV): leaf_n_change_rate => (RNULV - RNTLV - RNLDLV) ~ track(u"g/m^2/d")

    "Net stem nitrogen rate after uptake and translocation."
    RNST(RNUST, RNTST): stem_n_change_rate         => (RNUST - RNTST)          ~ track(u"g/m^2/d")

    "Net living root nitrogen rate after uptake, translocation, and death loss."
    RNRT(RNURT, RNTRT, RNLDRT): root_n_change_rate => (RNURT - RNTRT - RNLDRT) ~ track(u"g/m^2/d")

    "Nitrogen retained in living leaves after uptake, translocation, and loss."
    ANLV(RNLV): leaf_nitrogen_amount    ~ accumulate(init = ANLVI, u"g/m^2")

    "Nitrogen retained in stems after uptake and translocation."
    ANST(RNST): stem_nitrogen_amount    ~ accumulate(init = ANSTI, u"g/m^2")

    "Nitrogen retained in living roots after uptake, translocation, and death loss."
    ANRT(RNRT): root_nitrogen_amount    ~ accumulate(init = ANRTI, u"g/m^2")

    "Nitrogen accumulated in storage organs after the translocation stage begins."
    ANSO(RNSO): storage_nitrogen_amount ~ accumulate(init = ANSOI, u"g/m^2")

    "Nitrogen contained in living leaves, stems, roots, and storage organs."
    NCROP(ANLV, ANST, ANRT, ANSO): total_crop_nitrogen => (ANLV + ANST + ANRT + ANSO) ~ track(u"g/m^2")

    """
    Lower of the water and nitrogen status indices.

    This is retained with the other whole crop state and efficiency diagnostics
    named in the source.
    """
    RNW(TRANRF): combined_water_nitrogen_stress_factor   => TRANRF            ~ track(max = NNI)

    "Total aboveground dry mass in leaves, stems, and storage organs."
    TAGBM(WLV, WST, WSO): total_above_ground_biomass     => (WLV + WST + WSO) ~ track(u"g/m^2")

    "Aboveground dry mass retained under the inherited LINTUL identifier `WAD`."
    WAD(TAGBM): above_ground_biomass                     => TAGBM             ~ track(u"g/m^2")

    "Whole crop dry mass including living roots and aboveground organs."
    WTR(WRT, TAGBM): total_biomass                       => (WRT + TAGBM)     ~ track(u"g/m^2")

    "Daily dry matter production divided by incident photosynthetically active radiation."
    LUECAL(GTOTAL, PAR): calculated_light_use_efficiency => (GTOTAL / PAR)    ~ track(u"g/MJ")

    "Cumulative incident photosynthetically active radiation."
    CUMPAR(PAR): cumulative_par                         ~ accumulate(u"MJ/m^2")

    "Nitrogen contained in living leaves, stems, and storage organs."
    NTAG(ANLV, ANST, ANSO): above_ground_nitrogen     => (ANLV + ANST + ANSO) ~ track(u"g/m^2")

    "Flag allowing aboveground biomass to be used as a concentration denominator."
    TAGNZ(TAGBM): nonzero_above_ground_biomass        => (TAGBM != 0u"g/m^2") ~ flag

    "Safe aboveground biomass divisor corresponding to source `NOTNUL(TAGBM)`."
    TAGNN(TAGBM): source_notnul_above_ground_biomass  => TAGBM                ~ track(u"g/m^2", init = 1u"g/m^2", when = TAGNZ)

    "Nitrogen mass fraction of living aboveground biomass."
    NTAC(NTAG, TAGNN): above_ground_nitrogen_fraction => (NTAG / TAGNN)       ~ track

    "Ratio of living root biomass to total living aboveground biomass."
    RTSH(WRT, TAGBM): root_to_shoot_ratio             => (WRT / TAGBM)        ~ track

    "Nonzero green leaf biomass."
    WLVGNZ(WLVG): nonzero_green_leaf_biomass => (WLVG != 0u"g/m^2") ~ flag

    "Nonzero stem biomass."
    WSTNZ(WST): nonzero_stem_biomass         => (WST != 0u"g/m^2")  ~ flag

    "Nonzero root biomass."
    WRTNZ(WRT): nonzero_root_biomass         => (WRT != 0u"g/m^2")  ~ flag

    "Nonzero storage biomass."
    WSONZ(WSO): nonzero_storage_biomass      => (WSO != 0u"g/m^2")  ~ flag

    "Safe green leaf biomass divisor corresponding to source `NOTNUL(WLVG)`."
    WLVGNN(WLVG): source_notnul_green_leaf_biomass => WLVG ~ track(u"g/m^2", init = 1u"g/m^2", when = WLVGNZ)

    "Safe stem biomass divisor corresponding to source `NOTNUL(WST)`."
    WSTNN(WST): source_notnul_stem_biomass         => WST  ~ track(u"g/m^2", init = 1u"g/m^2", when = WSTNZ)

    "Safe root biomass divisor corresponding to source `NOTNUL(WRT)`."
    WRTNN(WRT): source_notnul_root_biomass         => WRT  ~ track(u"g/m^2", init = 1u"g/m^2", when = WRTNZ)

    "Safe storage biomass divisor corresponding to source `NOTNUL(WSO)`."
    WSONN(WSO): source_notnul_storage_biomass      => WSO  ~ track(u"g/m^2", init = 1u"g/m^2", when = WSONZ)

    "Nitrogen mass fraction of living green leaf biomass."
    NFLV(ANLV, WLVGNN): leaf_n_fraction       => (ANLV / WLVGNN) ~ track

    "Nitrogen mass fraction of living stem biomass."
    NFST(ANST, WSTNN): stem_n_fraction        => (ANST / WSTNN)  ~ track

    "Nitrogen mass fraction of living root biomass."
    NFRT(ANRT, WRTNN): root_n_fraction        => (ANRT / WRTNN)  ~ track

    "Nitrogen mass fraction of storage-organ biomass."
    NFSO(ANSO, WSONN): storage_n_fraction     => (ANSO / WSONN)  ~ track

    "Cumulative water incorporated as the modeled root zone expands."
    TEXPLO(EXPLOR): fst_cumulative_root_exploration_water ~ accumulate(u"mm")

    "Cumulative actual soil evaporation retained for source balance checks."
    TEVAP(EVAP): fst_cumulative_evaporation               ~ accumulate(u"mm")

    "Cumulative actual crop transpiration retained for source balance checks."
    TTRAN(TRAN): fst_cumulative_transpiration             ~ accumulate(u"mm")

    "Cumulative surface runoff retained for source balance checks."
    TRUNOF(RUNOFF): fst_cumulative_runoff                 ~ accumulate(u"mm")

    "Cumulative irrigation supplied by the active water-management relation."
    TIRRIG(IRRIG): fst_cumulative_irrigation              ~ accumulate(u"mm")

    "Cumulative drainage below the modeled root zone."
    TDRAIN(DRAIN): fst_cumulative_drainage                ~ accumulate(u"mm")

    """
    Residual of the integrated water balance written in the source form.

    A value near zero indicates closure among initial storage, cumulative
    inputs, final storage, and cumulative losses.
    """
    WATBAL(WA, WAI, TRAIN, TEXPLO, TIRRIG, TRUNOF, TTRAN, TEVAP, TDRAIN): fst_water_balance_error => begin
        WA - WAI - TRAIN - TEXPLO - TIRRIG + TRUNOF + TTRAN + TEVAP + TDRAIN
    end ~ track(u"mm")

    "Source water balance residual retained under the inherited identifier `WBAL`."
    WBAL(WATBAL): water_balance_error => WATBAL ~ track(u"mm")

    "Cumulative daily dry matter production used in the carbon balance."
    GTSUM(GTOTAL): cumulative_growth             ~ accumulate(u"g/m^2")

    "Residual between cumulative growth and living plus dead biomass pools."
    CBALAN(GTSUM, WRTLI, WLVGI, WSTI, WSOI, WTR, WDRT): carbon_balance_error => begin
        GTSUM + WRTLI + WLVGI + WSTI + WSOI - WTR - WDRT
    end ~ track(u"g/m^2")

    "Cumulative nitrogen transferred from the inorganic soil pool to the crop."
    NUPTT(NUPTR): cumulative_n_uptake            ~ accumulate(u"g/m^2")

    "Cumulative nitrogen removed with dead green leaf biomass."
    NLOSSL(RNLDLV): cumulative_leaf_n_loss       ~ accumulate(u"g/m^2")

    "Cumulative nitrogen removed with dead living-root biomass."
    NLOSSR(RNLDRT): cumulative_root_n_loss       ~ accumulate(u"g/m^2")

    "Cumulative mineralized nitrogen added to the inorganic soil pool."
    NMINERT(RTMIN): cumulative_n_mineralization  ~ accumulate(u"g/m^2")

    "Cumulative gross fertilizer nitrogen supplied at the management boundary."
    NFERTT(NFERT): cumulative_fertilizer_n_input ~ accumulate(u"g/m^2")

    "Cumulative fertilizer nitrogen entering the soil pool after recovery loss."
    NFERTRECT(FERTNS, DELT): cumulative_recovered_fertilizer_n => (FERTNS / DELT) ~ accumulate(u"g/m^2")

    "Residual of initial crop nitrogen plus uptake against crop pools and losses."
    NBALAN(ANLVI, ANSTI, ANRTI, ANSOI, NUPTT, NCROP, NLOSSL, NLOSSR): crop_n_balance_error => begin
        ANLVI + ANSTI + ANRTI + ANSOI + NUPTT - NCROP - NLOSSL - NLOSSR
    end ~ track(u"g/m^2")

    "Residual of fertilizer recovery and mineralization against soil nitrogen and uptake."
    NSBALAN(FERTNPI, NFERTRECT, NMINERT, TNSOIL, NUPTT): soil_n_balance_error => begin
        FERTNPI + NFERTRECT + NMINERT - TNSOIL - NUPTT
    end ~ track(u"g/m^2")
end

"""
LINTUL-3 water and nitrogen limited spring wheat model in Cropbox.

LINTUL (Light INTerception and UtiLization) was developed within the Wageningen
crop modeling tradition. LINTUL-3 extends the potential and water limited
members of the family with soil and crop nitrogen dynamics. This implementation
reconstructs the selected FST (Fortran Simulation Translator) test program for
spring wheat using Flevoland parameters. The published LINTUL3 article below
describes the wider model family through a rice application and is not the
exact crop parameterization reconstructed here.

# References

- Shibu, M. E., Leffelaar, P. A., van Keulen, H. and Aggarwal, P. K. (2010).
  LINTUL3, a simulation model for nitrogen-limited situations: Application to
  rice. *European Journal of Agronomy*, 32, 255-271.
  <https://doi.org/10.1016/j.eja.2010.01.003>
- van Ittersum, M. K. et al. (2003). On approaches and applications of the
  Wageningen crop models. *European Journal of Agronomy*, 18, 201-234.
  <https://doi.org/10.1016/S1161-0301(02)00106-5>
"""
@system Lintul3Model(Weather, Lintul3System, Controller)

"""
Source correspondence implementation of LINTUL-3 that retains the saved `DVS1`
state in the selected FST program's `SUBDVS` routine.

FST evaluates `SUBDVS` before advancing `TSUM`. The update condition therefore
uses the next thermal time state so that `DVS1` remains at the last value
observed before the anthesis threshold is crossed. This system is retained for
source correspondence tests, whereas `Lintul3Model` uses the continuous
development relation described above.
"""
@system Lintul3SavedStateSystem(Lintul3System) begin

    """
    Condition reproducing the selected FST `SUBDVS` update boundary.

    The next thermal time state is tested before updating `DVS1`. Consequently,
    the saved preanthesis component retains its last value below one when a
    step crosses the anthesis threshold.
    """
    DVS1UP(TSUM, RTSUMP, DELT, TSUMAN): saved_dvs1_update_condition => begin
        TSUM + RTSUMP * DELT <= TSUMAN
    end ~ flag

    """
    Development stage rate integrated only while the source `SUBDVS` boundary
    remains active.
    """
    RDVS1(RTSUMP, TSUMAN): saved_dvs1_rate => begin
        RTSUMP / TSUMAN
    end ~ track(u"d^-1", when = DVS1UP)

    "Saved preanthesis development component retained for FST correspondence."
    DVS1(RDVS1): pre_anthesis_development_stage ~ accumulate(init = DVSI)
end

"""
Source correspondence variant of the selected LINTUL-3 spring wheat model.

This runnable model retains the saved preanthesis `DVS1` behavior observed in
the FST `SUBDVS` routine. It is provided to compare the Cropbox reconstruction
with the executed source program and should not be interpreted as a fourth
LINTUL model or a separate crop parameterization.
"""
@system Lintul3SavedStateModel(Weather, Lintul3SavedStateSystem, Controller)
