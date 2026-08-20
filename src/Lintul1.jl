"""
Crop processes represented by the selected LINTUL-1 spring wheat program.

The system describes potential production under nonlimiting water and nitrogen.
It includes thermal time development, Beer-Lambert light interception, dry
matter production with constant light use efficiency, partitioning that varies
with development, leaf area growth, and leaf senescence. Weather provision and
run control are supplied when this system is composed as `Lintul1Model`.
"""
@system Lintul1System begin
    "Daily mean air temperature supplied by the shared weather system."
    TEMP ~ hold

    "Daily global solar radiation supplied by the shared weather system."
    RAD  ~ hold

    "Active clock step corresponding to the FST integration interval `DELT`."
    DELT(context.clock.step): integration_timestep ~ preserve(u"d")

    "Daily mean temperature retained under the FST identifier `DAVTMP`."
    DAVTMP(TEMP): daily_average_temperature => TEMP ~ track(u"°C")

    "Daily global radiation retained under the FST identifier `DTR`."
    DTR(RAD): daily_global_radiation        => RAD  ~ track(u"MJ/m^2/d")

    """
    Photosynthetically active radiation available above the canopy.

    The selected source assumes that PAR is one half of daily global solar
    radiation.
    """
    PAR(DTR): photosynthetically_active_radiation => (0.5DTR) ~ track(u"MJ/m^2/d")

    "Base temperature for thermal time accumulation."
    TBASE: base_temperature                    => 0            ~ preserve(parameter, u"°C")

    "Thermal time required from emergence to anthesis."
    TSUMAN: thermal_time_to_anthesis           => 1110         ~ preserve(parameter, u"K*d")

    "Thermal time required from anthesis to physiological maturity."
    TSUMMT: thermal_time_anthesis_to_maturity  => 970          ~ preserve(parameter, u"K*d")

    "Dry matter produced per unit of intercepted photosynthetically active radiation."
    LUE: light_use_efficiency                  => 3            ~ preserve(parameter, u"g/MJ")

    "Canopy light extinction coefficient in the Beer-Lambert interception relation."
    K: extinction_coefficient                  => 0.6          ~ preserve(parameter)

    "Green leaf area per unit of green leaf dry mass."
    SLA: specific_leaf_area                    => 0.022        ~ preserve(parameter, u"m^2/g")

    "Relative rate coefficient for exponential juvenile leaf area expansion."
    RGRL: relative_leaf_growth_rate            => 0.009        ~ preserve(parameter, u"K^-1*d^-1")

    "Thermal time limit of the juvenile exponential leaf area phase."
    TSUMJ: juvenile_lai_thermal_time_limit     => 330          ~ preserve(parameter, u"K*d")

    "Leaf area index limit of the juvenile exponential leaf area phase."
    LAIJ: juvenile_lai_threshold               => 0.75         ~ preserve(parameter)

    "Leaf area index above which canopy shading can cause leaf death."
    LAICR: self_shading_lai_threshold          => 4            ~ preserve(parameter)

    "Upper bound on the relative leaf death rate caused by canopy shading."
    RDRSHM: maximum_shading_death_rate         => 0.03         ~ preserve(parameter, u"d^-1")

    "Leaf area index established at the source simulation boundary."
    LAII: initial_leaf_area_index              => 0.012        ~ preserve(parameter)

    "Initial green leaf dry mass inferred from initial leaf area index and SLA."
    WLVI(LAII, SLA): initial_green_leaf_weight => (LAII / SLA) ~ preserve(u"g/m^2")

    "Initial living root dry mass."
    WRTI: initial_root_weight                  => 0            ~ preserve(parameter, u"g/m^2")

    "Initial living stem dry mass."
    WSTI: initial_stem_weight                  => 0            ~ preserve(parameter, u"g/m^2")

    "Initial storage organ dry mass."
    WSOI: initial_storage_organ_weight         => 0            ~ preserve(parameter, u"g/m^2")

    "Tabulated fraction of daily dry matter allocated to roots over thermal time."
    FRTTB => [
           0  0.5
         110  0.5
         275  0.34
         555  0.12
         780  0.07
        1055  0.03
        1160  0.02
        1305  0
        2500  0
    ] ~ interpolate(parameter, knotunit = u"K*d")

    "Tabulated fraction of daily dry matter allocated to green leaves over thermal time."
    FLVTB => [
           0  0.33
         110  0.33
         275  0.46
         555  0.44
         780  0.14
        1055  0
        2500  0
    ] ~ interpolate(parameter, knotunit = u"K*d")

    "Tabulated fraction of daily dry matter allocated to stems over thermal time."
    FSTTB => [
           0  0.17
         110  0.17
         275  0.2
         555  0.44
         780  0.79
        1055  0.97
        1160  0
        2500  0
    ] ~ interpolate(parameter, knotunit = u"K*d")

    "Tabulated fraction of daily dry matter allocated to storage organs over thermal time."
    FSOTB => [
           0  0
        1055  0
        1160  0.98
        1305  1
        2500  1
    ] ~ interpolate(parameter, knotunit = u"K*d")

    "Tabulated relative leaf death rate as a function of daily mean temperature."
    RDRT  => [
        -10  0.03
         10  0.03
         15  0.04
         30  0.09
         50  0.09
    ] ~ interpolate(parameter, u"d^-1", knotunit = u"°C")

    "Daily temperature increment above the base temperature."
    DTEFF(DAVTMP, TBASE): effective_temperature => (DAVTMP - TBASE) ~ track(u"K", min = 0)

    "Effective temperature rate supplied to the thermal time state."
    RTSUM(DTEFF): thermal_time_rate             => DTEFF            ~ track(u"K")

    "Effective thermal time accumulated from the boundary after emergence."
    TSUM(RTSUM): thermal_time                                      ~ accumulate(u"K*d")

    "Source termination flag raised after the selected program's thermal time limit."
    FINISH(TSUM): source_finish_condition => (TSUM > 2080u"K*d") ~ flag

    "Flag distinguishing the preanthesis and postanthesis development phases."
    PRE_ANTHESIS(TSUM, TSUMAN): before_anthesis         => (TSUM < TSUMAN)  ~ flag

    "Development stage rate before anthesis."
    DVR1(DTEFF, TSUMAN): pre_anthesis_development_rate  => (DTEFF / TSUMAN) ~ track(u"d^-1", when = PRE_ANTHESIS)

    "Development stage rate after anthesis."
    DVR2(DTEFF, TSUMMT): post_anthesis_development_rate => (DTEFF / TSUMMT) ~ track(u"d^-1", when = !PRE_ANTHESIS)

    "Development stage rate assembled from the mutually exclusive phase rates."
    DVR(DVR1, DVR2): development_rate                   => (DVR1 + DVR2)    ~ track(u"d^-1")

    """
    Preanthesis development stage, bounded at one at anthesis.

    Together with `DVS2`, this expresses the conventional LINTUL development
    scale of zero at emergence, one at anthesis, and two at maturity.
    """
    DVS1(TSUM, TSUMAN): pre_anthesis_development_stage          => (TSUM / TSUMAN)            ~ track(max = 1)

    "Postanthesis development stage, starting from zero at anthesis."
    DVS2(TSUM, TSUMAN, TSUMMT): post_anthesis_development_stage => ((TSUM - TSUMAN) / TSUMMT) ~ track(min = 0)

    "Combined development stage on the scale from emergence through anthesis to maturity."
    DVS(DVS1, DVS2): development_stage                          => (DVS1 + DVS2)              ~ track

    """
    Fraction of incident radiation intercepted by the canopy.

    This is the Beer-Lambert relation `1 - exp(-K * LAI)` used by the source
    radiation use efficiency model.
    """
    FINT(K, LAI): intercepted_fraction     => (1 - exp(-K * LAI)) ~ track

    "Daily intercepted PAR obtained by multiplying incident PAR by `FINT`."
    PARINT(PAR, FINT): intercepted_par     => (PAR * FINT)        ~ track(u"MJ/m^2/d")

    """
    Potential daily dry matter production.

    The LINTUL radiation use efficiency relation is
    `GTOTAL = LUE * PARINT`.
    """
    GTOTAL(PARINT, LUE): daily_growth_rate => (LUE * PARINT)      ~ track(u"g/m^2/d")

    """
    Thermal time supplied to the source partitioning tables.

    The source finish condition keeps this query within the selected table
    domain. A configured replacement table must cover the active thermal time
    range, so an incomplete table fails visibly instead of being silently
    clamped to the original endpoints.
    """
    TSUM_AFGEN(TSUM): partitioning_table_thermal_time         => TSUM                    ~ track(u"K*d")

    "Fraction of daily growth allocated to living roots."
    FRT(FRTTB, TSUM_AFGEN): root_allocation_fraction          => FRTTB(TSUM_AFGEN)       ~ track

    "Fraction of daily growth allocated to green leaves."
    FLV(FLVTB, TSUM_AFGEN): leaf_allocation_fraction          => FLVTB(TSUM_AFGEN)       ~ track

    "Fraction of daily growth allocated to stems."
    FST(FSTTB, TSUM_AFGEN): stem_allocation_fraction          => FSTTB(TSUM_AFGEN)       ~ track

    "Fraction of daily growth allocated to storage organs."
    FSO(FSOTB, TSUM_AFGEN): storage_organ_allocation_fraction => FSOTB(TSUM_AFGEN)       ~ track

    "Fraction of daily growth allocated to all aboveground organs."
    FSH(FRT): shoot_allocation_fraction                       => (1 - FRT)               ~ track

    "Sum of organ allocation fractions used to diagnose partitioning closure."
    FRACT(FRT, FLV, FST, FSO): allocation_sum                 => (FRT + FLV + FST + FSO) ~ track

    "Gross green leaf growth before senescence losses."
    GLV(GTOTAL, FLV): gross_leaf_growth_rate     => (GTOTAL * FLV) ~ track(u"g/m^2/d")

    "Daily dry matter growth allocated to roots."
    RWRT(GTOTAL, FRT): root_growth_rate          => (GTOTAL * FRT) ~ track(u"g/m^2/d")

    "Daily dry matter growth allocated to stems."
    RWST(GTOTAL, FST): stem_growth_rate          => (GTOTAL * FST) ~ track(u"g/m^2/d")

    "Daily dry matter growth allocated to storage organs."
    RWSO(GTOTAL, FSO): storage_organ_growth_rate => (GTOTAL * FSO) ~ track(u"g/m^2/d")

    """
    Flag selecting linear extrapolation below the leaf death table.

    The explicit below, inside, and above domains reproduce the selected FST
    `AFGEN` terminal segment extrapolation while keeping the table itself declarative.
    """
    RDRT_LO(TEMP): rdr_temperature_below_table => (TEMP < -10u"°C") ~ flag

    "Flag selecting linear extrapolation above the leaf death table."
    RDRT_HI(TEMP): rdr_temperature_above_table => (TEMP > 50u"°C")  ~ flag

    "Flag selecting interpolation within the leaf death table domain."
    RDRT_IN(RDRT_LO, RDRT_HI): rdr_temperature_inside_table => begin
        !RDRT_LO && !RDRT_HI
    end ~ flag

    "Relative leaf death rate interpolated within the temperature table domain."
    RDRTMP_TABLE(RDRT, TEMP): interpolated_temperature_leaf_death_rate => RDRT(TEMP)  ~ track(u"d^-1", when = RDRT_IN)

    "Relative leaf death rate below the table domain."
    RDRTMP_LOW: low_temperature_leaf_death_extrapolation        => 0.03u"d^-1" ~ track(u"d^-1", when = RDRT_LO)

    "Relative leaf death rate above the table domain."
    RDRTMP_HIGH: high_temperature_leaf_death_extrapolation      => 0.09u"d^-1" ~ track(u"d^-1", when = RDRT_HI)

    "Leaf death rate determined by temperature across the three table domains."
    RDRTMP(RDRTMP_TABLE, RDRTMP_LOW, RDRTMP_HIGH): temperature_leaf_death_rate => begin
        RDRTMP_TABLE + RDRTMP_LOW + RDRTMP_HIGH
    end ~ track(u"d^-1")

    "Aging rate determined by temperature and activated after anthesis."
    RDRDV(RDRTMP): aging_leaf_death_rate => RDRTMP ~ track(u"d^-1", when = !PRE_ANTHESIS)

    "Relative leaf death rate caused by canopy shading."
    RDRSH(LAI, LAICR, RDRSHM): shading_leaf_death_rate => begin
        RDRSHM * (LAI - LAICR) / LAICR
    end ~ track(u"d^-1", min = 0, max = RDRSHM)

    "Dominant relative leaf death rate from aging or canopy shading."
    RDR(RDRDV): relative_leaf_death_rate        => RDRDV        ~ track(u"d^-1", min = RDRSH)

    "Leaf area index lost per day through senescence."
    DLAI(LAI, RDR): leaf_area_death_rate        => (LAI * RDR)  ~ track(u"d^-1")

    "Green leaf dry mass lost per day through senescence."
    DLV(WLVG, RDR): leaf_biomass_death_rate     => (WLVG * RDR) ~ track(u"g/m^2/d")

    "Net daily change in green leaf dry mass."
    RWLVG(GLV, DLV): net_green_leaf_growth_rate => (GLV - DLV)  ~ track(u"g/m^2/d")

    "Living root biomass integrated from the net root growth rate."
    WRT(RWRT): root_biomass                           ~ accumulate(init = WRTI, u"g/m^2")

    "Living green leaf biomass after growth and senescence."
    WLVG(RWLVG): green_leaf_biomass                   ~ accumulate(init = WLVI, u"g/m^2")

    "Cumulative dead leaf biomass removed from the green leaf pool."
    WLVD(DLV): dead_leaf_biomass                      ~ accumulate(u"g/m^2")

    "Total living and dead leaf dry mass."
    WLV(WLVG, WLVD): total_leaf_biomass      => (WLVG + WLVD)     ~ track(u"g/m^2")

    "Living stem biomass integrated from allocated daily growth."
    WST(RWST): stem_biomass                           ~ accumulate(init = WSTI, u"g/m^2")

    "Storage organ biomass integrated from allocated daily growth."
    WSO(RWSO): storage_organ_biomass                  ~ accumulate(init = WSOI, u"g/m^2")

    "Total aboveground dry mass in leaves, stems, and storage organs."
    WAD(WLV, WST, WSO): above_ground_biomass => (WLV + WST + WSO) ~ track(u"g/m^2")

    "Whole crop dry mass including living roots."
    WTR(WRT, WAD): total_biomass             => (WRT + WAD)       ~ track(u"g/m^2")

    """
    Early exponential LAI growth expressed as an equivalent rate.

    The source applies exponential expansion over each integration interval.
    `DELT` makes that interval explicit so the relation remains consistent when
    the clock step changes.
    """
    GLAE(LAI, RGRL, DTEFF, DELT): exponential_lai_growth_rate => begin
        LAI * expm1(RGRL * DTEFF * DELT) / DELT
    end ~ track(u"d^-1")

    "LAI growth calculated from new leaf biomass and specific leaf area."
    GLAL(GLV, SLA): linear_lai_growth_rate => (GLV * SLA) ~ track(u"d^-1")

    "Flag selecting juvenile exponential LAI growth by thermal time and LAI."
    JUVLAI(TSUM, TSUMJ, LAI, LAIJ): juvenile_lai_condition => begin
        TSUM < TSUMJ && LAI < LAIJ
    end ~ flag

    "Flag requesting the source initial LAI when the current LAI is zero."
    SEEDLAI(LAI): source_lai_seeding_condition             => (LAI == 0) ~ flag

    "Flag allowing ordinary LAI growth after LAI has been established."
    REGLAI(SEEDLAI): regular_lai_growth_condition          => !SEEDLAI   ~ flag

    "Flag activating juvenile exponential LAI growth."
    JUVGROW(JUVLAI, REGLAI): active_juvenile_lai_growth    => begin
        JUVLAI && REGLAI
    end ~ flag

    "Flag activating mature LAI growth from new leaf biomass."
    MATGROW(JUVLAI, REGLAI): active_mature_lai_growth      => begin
        !JUVLAI && REGLAI
    end ~ flag

    "Juvenile LAI growth contribution."
    GLAIJ(GLAE): juvenile_lai_growth_rate             => GLAE          ~ track(u"d^-1", when = JUVGROW)

    "Mature LAI growth contribution derived from new leaf biomass."
    GLAIM(GLAL): mature_lai_growth_rate               => GLAL          ~ track(u"d^-1", when = MATGROW)

    "Equivalent rate that establishes the source initial LAI in one integration step."
    GLAI_SEED(LAII, DELT): emergence_lai_seeding_rate => (LAII / DELT) ~ track(u"d^-1", when = SEEDLAI)

    "Gross LAI growth assembled from seeding, juvenile, and mature contributions."
    GLAI(GLAIJ, GLAIM, GLAI_SEED): gross_lai_growth_rate => begin
        GLAIJ + GLAIM + GLAI_SEED
    end ~ track(u"d^-1")

    "Net LAI rate after subtracting leaf area death."
    RLAI(GLAI, DLAI): net_lai_change_rate => (GLAI - DLAI) ~ track(u"d^-1")

    "Green leaf area per unit ground area, integrated from growth and death rates."
    LAI(RLAI): leaf_area_index ~ accumulate
end

"""
LINTUL-1 potential production model in Cropbox.

LINTUL (Light INTerception and UtiLization) originated in the Wageningen crop
modeling tradition as a concise model of dry matter production from intercepted
radiation and a constant light use efficiency. LINTUL-1 represents the
potential production member of the model family, where water and nutrient
supply do not limit growth. This implementation reconstructs the selected FST
(Fortran Simulation Translator) spring wheat program.

# References

- Spitters, C. J. T. and Schapendonk, A. H. C. M. (1990). Evaluation of
  breeding strategies for drought tolerance in potato by means of crop growth
  simulation. *Plant and Soil*, 123, 193-203.
  <https://doi.org/10.1007/BF00011268>
- van Ittersum, M. K. et al. (2003). On approaches and applications of the
  Wageningen crop models. *European Journal of Agronomy*, 18, 201-234.
  <https://doi.org/10.1016/S1161-0301(02)00106-5>
"""
@system Lintul1Model(Weather, Lintul1System, Controller)
