using Mimi

@defcomp slrcost begin
    # Define all variables, parameters and indices used by this module

    # --- Indices ---
    ciam_country = Index()
    segments = Index()
    adaptPers = Index()

    # --- Region / segment mapping ---
    segID = Parameter(index = [segments])   # Unique segment numeric identifier
    xsc = Parameter{Dict{Int,Tuple{Int,Int,Int}}}()  # Region to segment mapping (dictionary) to keep track of which segments belong to each region
    rcp = Parameter{Int}()                   # RCP being run (metadata; not used in run)
    percentile = Parameter{Int}()            # Percentile of RCP being run (metadata; not used in run)
    ssp = Parameter{Int}()                   # SSP being used (0 for base case)

    # ---Time-related Parameters---
    tstep = Parameter()                     # Length of individual time-step (years)
    at = Parameter(index = [adaptPers])    # Array of time indices that mark starts of adaptation periods
    ntsteps = Parameter{Int}()              # Number of time-steps

    # ---Model Parameters ---
    fixed = Parameter{Bool}()               # Run model as fixed (T) or flexible (F) with respect to adaptation
    noRetreat = Parameter{Bool}()           # Default (F). If T, segments will either protect or not adapt.
    allowMaintain = Parameter{Bool}()       # Default F. If T, segments will have the option to maintain current defenses

    # ---Socioeconomic Parameters---
    popinput = Parameter{Int}()           # Input for population data source: 0 (default), 1 (Jones & O'Neill, 2016), 2 (Merkens et al, 2016)
    pop = Parameter(index = [time, ciam_country], unit = "million")           # Population of region (million people) (from MERGE or SSPs)
    refpopdens = Parameter(index = [ciam_country], unit = "persons/km2")         # Reference population density of region (people / km^2)
    rgn_ind_usa = Parameter{Int}()                     # Lookup parameter for USA region index, used in refpopdens and ypc
    #    for USA benchmark in vsl, rho and fundland calculations
    popdens = Parameter(index = [segments], unit = "persons/km2")           # Pop density of segment in time t = 1 (people/km^2)
    ypcc = Parameter(index = [time, ciam_country], unit = "US\$2010/yr/person")          # GDP per capita per region ($2010 per capita)

    popdens_seg = Variable(index = [time, segments], unit = "persons/km2")          # Population density of segment extrapolated forward in time (people / km^2)
    #popdens_seg_jones = Parameter(index=[time,segments], unit = "persons/km2")      # Holder for Jones and O'Neill population density (not currently supported)
    #popdens_seg_merkens = Parameter(index=[time,segments], unit = "persons/km2")    # Holder for Merkens et al population density (not currently supported)
    ypc_seg = Variable(index = [time, segments], unit = "US\$2010/yr/person")              # GDP per capita by segment ($2010 per capita) (multiplied by scaling factor)
    refA_R = Parameter(index = [segments])                # Reference retreat level of adaptation in 0 period
    refA_H = Parameter(index = [segments])                # Reference height for adaptation in 0 period

    # ---Land Parameters---
    landinput = Parameter{Bool}()                   # Set to T for FUND or F for GTAP

    gtapland = Parameter(index = [ciam_country], unit = "million US\$2010/km2")        # GTAP land value in 2007 (million 2010$ / km^2)
    dvbm = Parameter(unit = "million US\$2010/km2")                            # FUND value of OECD dryland per Darwin et al 1995 converted from $1995 ($2010M per sqkm) (5.376)
    kgdp = Parameter()                              # Capital output ratio (per MERGE) (3 by default)
    discountrate = Parameter()                      # Discount rate (0.04 by default)
    depr = Parameter()                              # Fraction of capital that has not been depreciated over adaptation period (retreat cases)


    landdata = Variable(index = [ciam_country], unit = "million US\$2010/km2")         # Takes on value of either fundland or gtapland
    fundland = Variable(index = [ciam_country], unit = "million US\$2010/km2")         # FUND land value in 1995 (calculated in run_timestep) (million 2010$ / km^2),

    rgn_ind_canada = Parameter{Int}()               # Region index for Canada (Used as reference for Greenland land appreciation)
    land_appr = Variable(index = [time, ciam_country])   # Land appreciation rate (calculated as regression by Yohe ref Abraham and Hendershott)
    coastland = Variable(index = [time, segments], unit = "million US\$2010/km2")  # Coastal land value (function of interior land value * scaling factor) ($2010M per sqkm)
    landvalue = Variable(index = [time, segments], unit = "million US\$2010/km2")  # Total endowment value of land ($2010M per sqkm)
    landrent = Variable(index = [time, segments], unit = "million US\$2010/km2/yr")      # Annual rental value of land ($2010M/sqkm/year)

    ρ = Variable(index = [time, ciam_country])           # Country-wide resilience parameter (logistic function related to GDP)
    capital = Variable(index = [time, segments], unit = "million US\$2010/km2")    # Total endowment value of capital stock (million $2010 / km^2)
    discountfactor = Variable(index = [time])         # Discount factor (derived from discount rate)

    # ---Coastal Parameters---
    length = Parameter(index = [segments], unit = "km")          # Segment length (km)

    # ---Protection Parameters---
    cci = Parameter(index = [ciam_country])
    pcfixed = Parameter()                   # Fraction of protection cost that is fixed (not variable in height) (0.3)
    mc = Parameter()                        # Maintenance cost (Hillen et al, 2010) (2%/yr)
    pc0 = Parameter(unit = "million US\$2010/km/m2")                       # Reference cost of protection (million 2010$ / km / vert m^2) (6.02 by default)

    # ---Retreat / No Adapt Parameters---
    mobcapfrac = Parameter()                # Fraction of capital that is mobile (0.25)
    movefactor = Parameter()                # Cost to relocate people as a factor of annual income (Tol 3x RMendelsohn 0.5x) (1)
    capmovefactor = Parameter()             # Cost to relocate mobile capital as a fraction of asset value (0.1)
    democost = Parameter()                  # Cost to demolish immobile capital as fraction of asset (0.05)

    # # ---Surge Exposure Parameters---
    # Protection case
    psig0 = Parameter(index = [segments])
    psig0coef = Parameter(index = [segments])
    psigA = Parameter(index = [segments])                       # psigA in GAMS code
    psigB = Parameter(index = [segments])                       # psigB in GAMS code

    # Retreat / No Adapt Cases
    rsig0 = Parameter(index = [segments])
    rsigA = Parameter(index = [segments])                       # rsigA in GAMS code
    rsigB = Parameter(index = [segments])                       # rsigB in GAMS code

    # ---Storm damage parameters---
    floodmortality = Parameter()                # Flood deaths as percent of exposed population; (Jonkman Vrijling 2008) (0.01)

    # if TRUE, then VSL is exogenously calculated for each country and set in vsl_ciam_country, then each segment is set by looking up its country
    # if FALSE, then VSL is endogenously calculated using vslel and vslmult as well as other endogenous socioeconomics for each segment
    vsl_exogenous = Parameter{Bool}(default=true) 

    vsl_ciam_country = Parameter(index = [time, ciam_country], unit = "million US\$2010/yr") # Value of statistical life (million 2010$) (only used for exogenous calculation of vsl)
    vslel = Parameter(default = 0.5)    # Elasticity of vsl (0.5) (only used for endogenous calculation of vsl)
    vslmult = Parameter(default = 216)  # multiplier on USA GDP (216)(only used for endogenous calculation of vsl)

    vsl = Variable(index = [time, segments], unit = "million US\$2010/yr")     # Value of statistical life (million 2010$)

    # ---Wetland Loss Parameters---
    wvbm = Parameter(default = 0.376, unit = "million US\$2010/km2/yr")                                      # Annual value of wetland services (million 2010$ / km^2 / yr); (Brander et al 2006)  (0.376)
    wetland = Parameter(index = [segments], unit = "km2")                # Initial wetland area in coastal segment (km^2)
    wmaxrate = Parameter(default = 0.01, unit = "m/yr")                                  # Maximum rate of wetland accretion (m per yr) per Kirwan et al 2010 (0.01)
    wvel = Parameter(default = 1.16)                                      # income elasticity of wetland value (1.16) (Brander et al, 2006)
    wvpdl = Parameter(default = 0.47)                                     # Population density elasticity of wetland value (0.47) (Brander et al, 2006)

    wetlandservice = Variable(index = [time, ciam_country])      # Annual value of wetland services adjusted for income and density (Brander et al 2006) ($2010M/km^2/year)
    wetlandloss = Variable(index = [time, segments])        # Fractional loss of wetland due to slr


    # ---Sea Level Rise Parameters---
    lslr = Parameter(index = [time, segments], unit = "m")                # Local sea level rise (m)
    adaptoptions = Parameter(index = [6])                     # Index of available adaptation levels for protect and retreat (0 is no adaptation)
    surgeexposure = Parameter{Float64}(index = [segments, 5])# Storm surge exposure levels (corresponding to each designated adaptation option)


    # ---Coastal Area Parameters---
    area1 = Parameter(index = [segments])
    area2 = Parameter(index = [segments])
    area3 = Parameter(index = [segments])
    area4 = Parameter(index = [segments])
    area5 = Parameter(index = [segments])
    area6 = Parameter(index = [segments])
    area7 = Parameter(index = [segments])
    area8 = Parameter(index = [segments])
    area9 = Parameter(index = [segments])
    area10 = Parameter(index = [segments])
    area11 = Parameter(index = [segments])
    area12 = Parameter(index = [segments])
    area13 = Parameter(index = [segments])
    area14 = Parameter(index = [segments])
    area15 = Parameter(index = [segments])
    areaparams = Variable(index = [segments, 15])           # Nothing is computed; this is just a convenient container for area params

    coastArea = Variable(index = [time, segments], unit = "km2")            # Coast area inundated (km^2)

    # ---Intermediate Variables---
    WetlandNoAdapt = Variable(index = [time, segments])
    FloodNoAdapt = Variable(index = [time, segments])
    StormCapitalNoAdapt = Variable(index = [time, segments])
    StormPopNoAdapt = Variable(index = [time, segments])
    RelocateNoAdapt = Variable(index = [time, segments])
    StormLossNoAdapt = Variable(index = [time, segments])
    DryLandLossNoAdapt = Variable(index = [time, segments])

    Construct = Variable(index = [time, segments, 5])
    WetlandProtect = Variable(index = [time, segments])
    StormCapitalProtect = Variable(index = [time, segments, 5])
    StormPopProtect = Variable(index = [time, segments, 5])
    StormLossProtect = Variable(index = [time, segments, 5])
    FloodProtect = Variable(index = [time, segments])

    WetlandRetreat = Variable(index = [time, segments])
    StormCapitalRetreat = Variable(index = [time, segments, 6])
    StormPopRetreat = Variable(index = [time, segments, 6])
    StormLossRetreat = Variable(index = [time, segments, 6])
    FloodRetreat = Variable(index = [time, segments, 6])
    RelocateRetreat = Variable(index = [time, segments, 6])
    DryLandLossRetreat = Variable(index = [time, segments, 6])
    coastAreaRetreat = Variable(index = [time, segments, 6])
    coastAreaNoAdapt = Variable(index = [time, segments])

    # --- Decision Variables --- (evaluated brute force)
    H = Variable(index = [time, segments, 5], unit = "m")       # Height of current sea wall, no retreat (m)
    R = Variable(index = [time, segments, 6], unit = "m")       # Retreat perimeter (m)
    SIGMA = Variable(index = [time, segments, 12])  # Expected value of effective exposure area for over-topping surge (all cases)
    # Order of sigma values: 1 no adapt case, 6 retreat cases, 5 protect cases in ascending order

    # ---Outcome Variables---
    OptimalH = Variable(index = [time, segments], unit = "m")               # m; Holder to track height built across timesteps (cumulative)
    OptimalR = Variable(index = [time, segments], unit = "m")               # m; Holder to track retreat radius across timesteps (cumulative)
    WetlandLossOptimal = Variable(index = [time, segments], unit = "km2")  # km2; Cumulative wetland loss from optimal decision
    DryLandLossOptimal = Variable(index = [time, segments], unit = "km2") # km2; Cumulative loss of dry land from optimal decision

    # DrylandLost = Variable(index=[time,segments])            # km2; container to track cumulative lost dryland
    WetlandLost = Variable(index = [time, segments], unit = "km2")            # km2; container to track cumulative lost wetland

    NoAdaptCost = Variable(index = [time, segments], unit = "billion US\$2010/yr")         # Cost of not adapting (e.g. reactive retreat) (2010$)
    ProtectCost = Variable(index = [time, segments, 5], unit = "billion US\$2010/yr")      # Total cost of protection at each level
    RetreatCost = Variable(index = [time, segments, 6], unit = "billion US\$2010/yr")      # Total cost of retreat at each level
    OptimalRetreatLevel = Variable(index = [time, segments])
    OptimalProtectLevel = Variable(index = [time, segments])
    OptimalCost = Variable(index = [time, segments], unit = "billion US\$2010/yr")          # Optimal cost based on NPV relative to start of adaptation period
    OptimalLevel = Variable(index = [time, segments])         # Fixed optimal level (1,10,100,1000,10000)
    OptimalOption = Variable(index = [time, segments])        # Fixed adaptation decision (-1 - protect, -2 - retreat, -3 - no adapt)
    NPVRetreat = Variable(index = [time, segments, 6])
    NPVProtect = Variable(index = [time, segments, 5])
    NPVNoAdapt = Variable(index = [time, segments])
    NPVOptimal = Variable(index = [segments])               # NPV of cost of optimal decisions relative to t=1
    NPVOptimalTotal = Variable()                            # Total NPV of all segments from optimal decision
    StormLossOptimal = Variable(index = [time, segments], unit = "persons")  # Cumulative expected loss of life (num people) from storm surges from optimal decision

    # ---Subcategories of Optimal Choice----
    OptimalStormCapital = Variable(index = [time, segments])
    OptimalStormPop = Variable(index = [time, segments])
    OptimalConstruct = Variable(index = [time, segments])
    OptimalWetland = Variable(index = [time, segments])
    OptimalFlood = Variable(index = [time, segments])
    OptimalRelocate = Variable(index = [time, segments])


    function run_timestep(p, v, d, t)
        # This is a workaround for a type instability that should be fixed in Mimi.jl
        d_ciam_country = d.ciam_country::Vector{Int}
        d_segments = d.segments::Vector{Int}

        ti1 = TimestepIndex(1) # used a lot
        # In first period, initialize all non-adaptation dependent intermediate variables for all timesteps
        if is_first(t)
            #  1. Initialize non-region dependent intermediate variables
            for i in collect(1:Int(p.ntsteps))
                v.discountfactor[TimestepIndex(i)] = 1 / (1 + p.discountrate)^(p.tstep * (i - 1))
            end

            # 2. Initialize region-dependent intermediate variables
            for r in d_ciam_country
                # Determine land input value (true = FUND, false = GTAP)
                if p.landinput
                    v.fundland[r] = min(p.dvbm, max(0.005, p.dvbm * p.ypcc[t, r] * p.refpopdens[r] / (p.ypcc[ti1, p.rgn_ind_usa] * p.refpopdens[p.rgn_ind_usa])))
                    v.landdata[r] = v.fundland[r]
                else
                    v.landdata[r] = p.gtapland[r]
                end

                # Calculate regional wetland service, resilience (rho), and land appreciation variables for the first period and
                #   subsequent periods
                v.wetlandservice[t, r] = p.wvbm * ((p.ypcc[t, r] / p.ypcc[ti1, p.rgn_ind_usa])^p.wvel * (p.refpopdens[r] / 27.59)^p.wvpdl)
                v.ρ[t, r] = p.ypcc[t, r] / (p.ypcc[t, r] + p.ypcc[ti1, p.rgn_ind_usa])
                v.land_appr[t, r] = 1.0

                for i in collect(2:Int(p.ntsteps))
                    ti = TimestepIndex(i)
                    tim1 = TimestepIndex(i - 1)
                    v.land_appr[ti, r] = v.land_appr[tim1, r] * exp(0.565 * growthrate(p.ypcc[tim1, r], p.ypcc[ti, r]) + 0.313 * growthrate(p.pop[tim1, r], p.pop[ti, r]))
                    v.wetlandservice[ti, r] = v.land_appr[ti, r] * v.wetlandservice[ti1, r]
                    v.ρ[ti, r] = p.ypcc[ti, r] / (p.ypcc[ti, r] + p.ypcc[ti1, p.rgn_ind_usa])
                end
            end

            # 3. Initialize segment-dependent variables
            for m in d_segments
                rgn_ind = getregion(m, p.xsc)::Int # Identify the region the segment belongs to

                # Initialize first-period population density, coast area and surge parameters
                if p.popinput == 0
                    v.popdens_seg[t, m] = p.popdens[m]
                elseif p.popinput == 1
                    error("The `popinput` argument values of 1 and 2 are not supported at this time.  In the future they will indicate use of Jones and O'Neill 2016 or Merkens et al 2016 population data, respectively.")
                    # v.popdens_seg[t,m]=p.popdens_seg_jones[ti1,m]
                elseif p.popinput == 2
                    error("The `popinput` argument values of 1 and 2 are not supported at this time.  In the future they will indicate use of Jones and O'Neill 2016 or Merkens et al 2016 population data, respectively.")
                    # v.popdens_seg[t,m]=p.popdens_seg_merkens[ti1,m]
                end
                v.areaparams[m, :] = [p.area1[m] p.area2[m] p.area3[m] p.area4[m] p.area5[m] p.area6[m] p.area7[m] p.area8[m] p.area9[m] p.area10[m] p.area11[m] p.area12[m] p.area13[m] p.area14[m] p.area15[m]]

                # Greenland segments are treated differently
                if isgreenland(m, p.xsc)::Int == 1
                    v.ypc_seg[t, m] = 22642 * 1.01^1   # FLAG: assumes t is an index (1-20)
                    v.coastland[t, m] = (v.land_appr[ti1, p.rgn_ind_canada] * v.landdata[p.rgn_ind_canada]) * max(0.5, log(1 + v.popdens_seg[t, m]) / log(25))
                    v.landvalue[t, m] = min(v.coastland[t, m], (v.land_appr[ti1, p.rgn_ind_canada] * v.landdata[p.rgn_ind_canada]))
                else
                    v.ypc_seg[t, m] = p.ypcc[t, rgn_ind] * max(0.9, (v.popdens_seg[ti1, m] / 250.0)^0.05)
                    v.coastland[t, m] = max(0.5, log(1 + v.popdens_seg[t, m]) / log(25)) * (v.land_appr[t, rgn_ind] * v.landdata[rgn_ind])  # Interior * scaling factor
                    v.landvalue[t, m] = min(v.coastland[t, m], (v.land_appr[t, rgn_ind] * v.landdata[rgn_ind]))
                end

                # Calculate VSL
                if p.vsl_exogenous
                    v.vsl[t, m] = p.vsl_ciam_country[ti1, rgn_ind] # pull VSL from vsl_ciam_country parameter for the proper region
                else # endogenous calculation of VSL
                    if isgreenland(m, p.xsc)::Int == 1 # Greenland segments are treated differently
                        v.vsl[t, m] = 1e-6 * p.vslmult * p.ypcc[ti1, p.rgn_ind_usa] * (v.ypc_seg[t, m] / p.ypcc[ti1, p.rgn_ind_usa])^p.vslel
                    else
                        v.vsl[t, m] = 1e-6 * p.vslmult * p.ypcc[ti1, p.rgn_ind_usa] * (p.ypcc[t, rgn_ind] / p.ypcc[ti1, p.rgn_ind_usa])^p.vslel
                    end
                end

                v.capital[t, m] = p.kgdp * v.ypc_seg[t, m] * v.popdens_seg[t, m] * 1e-6
                v.coastArea[t, m] = calcCoastArea(view(v.areaparams, m, :), p.lslr[t, m])

                for i = 2:Int(p.ntsteps)
                    ti = TimestepIndex(i)
                    tim1 = TimestepIndex(i - 1)
                    if p.popinput == 0
                        v.popdens_seg[ti, m] = v.popdens_seg[tim1, m] * (1 + growthrate(p.pop[tim1, rgn_ind], p.pop[ti, rgn_ind]))
                    elseif p.popinput == 1
                        error("The `popinput` argument values of 1 and 2 are not supported at this time.  In the future they will indicate use of Jones and O'Neill 2016 or Merkens et al 2016 population data, respectively.")
                        # v.popdens_seg[ti,m]=p.popdens_seg_jones[ti,m]
                    elseif p.popinput == 2
                        error("The `popinput` argument values of 1 and 2 are not supported at this time.  In the future they will indicate use of Jones and O'Neill 2016 or Merkens et al 2016 population data, respectively.")
                        # v.popdens_seg[ti,m]=p.popdens_seg_merkens[ti,m]
                    end

                    # Special treatment for Greenland segments
                    if isgreenland(m, p.xsc)::Int == 1
                        v.ypc_seg[ti, m] = 22642 * 1.01^i   # FLAG: assumes i is an index (1-20)
                        v.coastland[ti, m] = (v.land_appr[ti, p.rgn_ind_canada] * v.landdata[p.rgn_ind_canada]) * max(0.5, log(1 + v.popdens_seg[ti, m]) / log(25))
                        v.landvalue[ti, m] = min(v.coastland[ti, m], (v.land_appr[ti, p.rgn_ind_canada] * v.landdata[p.rgn_ind_canada]))

                    else
                        v.ypc_seg[ti, m] = p.ypcc[ti, rgn_ind] * max(0.9, (v.popdens_seg[ti1, m] / 250.0)^0.05) # ypcc * popdens scaling factor
                        v.coastland[ti, m] = max(0.5, log(1 + v.popdens_seg[ti, m]) / log(25)) * (v.land_appr[ti, rgn_ind] * v.landdata[rgn_ind])
                        v.landvalue[ti, m] = min(v.coastland[ti, m], (v.land_appr[ti, rgn_ind] * v.landdata[rgn_ind]))

                    end

                    # Calculate VSL
                    if p.vsl_exogenous
                        v.vsl[ti, m] = p.vsl_ciam_country[ti, rgn_ind]
                    else # endogenously calculate VSL
                        if isgreenland(m, p.xsc)::Int == 1  # Special treatment for Greenland segments
                            v.vsl[ti, m] = 1e-6 * p.vslmult * p.ypcc[ti, p.rgn_ind_usa] * (v.ypc_seg[ti, m] / p.ypcc[ti, p.rgn_ind_usa])^p.vslel   
                        else
                            v.vsl[ti, m] = 1e-6 * p.vslmult * p.ypcc[ti, p.rgn_ind_usa] * (p.ypcc[ti, rgn_ind] / p.ypcc[ti, p.rgn_ind_usa])^p.vslel   
                        end
                    end

                    v.capital[ti, m] = p.kgdp * v.ypc_seg[ti, m] * v.popdens_seg[ti, m] * 1e-6
                    v.coastArea[ti, m] = calcCoastArea(view(v.areaparams, m, :), p.lslr[ti, m])
                    v.wetlandloss[tim1, m] = min(1, (localrate(p.lslr[tim1, m], p.lslr[ti, m], p.tstep) / p.wmaxrate)^2)


                end
                v.wetlandloss[TimestepIndex(p.ntsteps), m] = min(1, (localrate(p.lslr[TimestepIndex(p.ntsteps - 1), m], p.lslr[TimestepIndex(p.ntsteps), m], p.tstep) / p.wmaxrate)^2)

            end

        end

        if (gettime(t) in p.at)
            adapt_range = collect(1:length(p.adaptoptions))

            # Determine length of adaptation period ("atstep")
            g(c) = c == gettime(t)
            at_index = findall(g, p.at)[1] # Find index corresponding to adaptation period in p.at
            at_index_next = at_index + 1    # Find index corresponding to next adaptation period
            at_index_prev = max(1, at_index - 1)

            at_prev = Int(p.at[at_index_prev])      # TODO pick one - either index into vector or use the period, it's confusing

            if at_index_next <= length(p.at)
                atstep = (p.at[at_index_next] - p.at[at_index]) * p.tstep   # In years
                at_next = Int(p.at[at_index_next])
                last_t = at_next - 1
                last = 0
            else
                # Deal with special case of last adaptation period
                atstep = p.tstep * p.ntsteps - ((p.at[at_index] - 1) * p.tstep)
                at_next = p.ntsteps # Flag this assumes timesteps are indices not years (1:20 not 2010:2100)
                last_t = at_next
                last = 1
            end
            t_range = collect(gettime(t):last_t)

            for m in d_segments
                if atstep == 0
                else
                    rgn_ind = getregion(m, p.xsc)::Int

                    # ** Calculate No Adaptation Costs **
                    for i in t_range
                        ti = TimestepIndex(i)
                        tim1 = TimestepIndex(i - 1)
                        R_NoAdapt = max(0, p.lslr[ti, m])

                        # For initial state in SLR cases, make adaptation decision relative to baseline (refA_H or R)
                        #if p.rcp>0
                        if p.rcp >= 0
                            R_NoAdapt = max(R_NoAdapt, p.refA_H[m], p.refA_R[m])
                        end

                        # Incorporate any previous period adaptation
                        if p.fixed == false && !(is_first(t))

                            R_NoAdapt = max(R_NoAdapt, v.OptimalR[TimestepIndex(gettime(t) - 1), m])
                            v.WetlandNoAdapt[ti, m] = p.tstep * v.wetlandservice[ti, rgn_ind] * max(v.WetlandLossOptimal[TimestepIndex(gettime(t) - 1), m], v.wetlandloss[ti, m] * min(v.coastArea[ti, m], p.wetland[m]))
                            if i == gettime(t)
                                # For start of new adaptation period, take into account (lack of) retreat done in previous periods (i.e. if they protected instead)
                                # This results in double-costs for this period b/c no adaptation is set up to compute relative to t+1 lslr
                                v.coastAreaNoAdapt[ti, m] = calcCoastArea(view(v.areaparams, m, :), v.OptimalR[TimestepIndex(gettime(t) - 1), m])
                            else
                                v.coastAreaNoAdapt[ti, m] = calcCoastArea(view(v.areaparams, m, :), R_NoAdapt)
                            end

                        else
                            v.coastAreaNoAdapt[ti, m] = v.coastArea[ti, m]
                            v.WetlandNoAdapt[ti, m] = p.tstep * v.wetlandservice[ti, rgn_ind] * v.wetlandloss[ti, m] * min(v.coastArea[ti, m], p.wetland[m])
                        end


                        # Storm Costs
                        v.SIGMA[ti, m, 1] = p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, R_NoAdapt - p.lslr[ti, m]))) # expected value of exposure area
                        v.StormCapitalNoAdapt[ti, m] = p.tstep * (1 - v.ρ[ti, rgn_ind]) * v.SIGMA[ti, m, 1] * v.capital[ti, m]
                        v.StormPopNoAdapt[ti, m] = p.tstep * (1 - v.ρ[ti, rgn_ind]) * v.popdens_seg[ti, m] * v.vsl[ti, m] * p.floodmortality * v.SIGMA[ti, m, 1]


                        v.StormLossNoAdapt[ti, m] = p.tstep * (1 - v.ρ[ti, rgn_ind]) * v.popdens_seg[ti, m] * p.floodmortality * v.SIGMA[ti, m, 1]
                        if i == p.ntsteps
                            v.DryLandLossNoAdapt[ti, m] = max(0, v.coastAreaNoAdapt[ti, m]) # km^2
                        else
                            # In case of negative or decreasing slr, can assume that previous inundated area is reclaimed
                            v.DryLandLossNoAdapt[ti, m] = max(0, v.coastAreaNoAdapt[ti, m], v.coastArea[TimestepIndex(i + 1), m]) # includes future period loss and previous adaptation if applicable
                        end

                        # Flood and relocation costs
                        if i == p.ntsteps
                            v.FloodNoAdapt[ti, m] = v.FloodNoAdapt[tim1, m]
                            v.RelocateNoAdapt[ti, m] = v.RelocateNoAdapt[tim1, m]
                        else
                            v.FloodNoAdapt[ti, m] = p.tstep * v.landvalue[ti, m] * 0.04 * v.DryLandLossNoAdapt[ti, m] + max(0, v.DryLandLossNoAdapt[ti, m] - v.coastAreaNoAdapt[ti, m]) * (1 - p.mobcapfrac) * v.capital[ti, m]

                            v.RelocateNoAdapt[ti, m] = max(0, v.DryLandLossNoAdapt[ti, m] - v.coastAreaNoAdapt[ti, m]) * (5 * p.movefactor * v.ypc_seg[ti, m] * 1e-6 * v.popdens_seg[ti, m] + p.capmovefactor * p.mobcapfrac * v.capital[ti, m] + p.democost * (1 - p.mobcapfrac) * v.capital[ti, m])
                        end

                        # Put all costs into $Billions and divide by 10 to account for 10-y time step
                        v.WetlandNoAdapt[ti, m] = v.WetlandNoAdapt[ti, m] * 1e-4
                        if i < p.ntsteps # already occurred in previous timestep
                            v.FloodNoAdapt[ti, m] = v.FloodNoAdapt[ti, m] * 1e-4
                            v.RelocateNoAdapt[ti, m] = v.RelocateNoAdapt[ti, m] * 1e-4
                        end
                        v.StormCapitalNoAdapt[ti, m] = v.StormCapitalNoAdapt[ti, m] * 1e-4
                        v.StormPopNoAdapt[ti, m] = v.StormPopNoAdapt[ti, m] * 1e-4

                        v.NoAdaptCost[ti, m] = v.WetlandNoAdapt[ti, m] + v.FloodNoAdapt[ti, m] + v.RelocateNoAdapt[ti, m] + v.StormCapitalNoAdapt[ti, m] + v.StormPopNoAdapt[ti, m]


                    end

                    if is_first(t)
                        v.NPVNoAdapt[t, m] = sum([v.discountfactor[TimestepIndex(j)] * v.NoAdaptCost[TimestepIndex(j), m] * 10 for j in t_range])
                    else
                        # Compute NPV Relative to planner's perspective (discounting relative to time t)
                        v.NPVNoAdapt[t, m] = sum(v.discountfactor[TimestepIndex(findind(j, t_range))] * v.NoAdaptCost[TimestepIndex(j), m] * 10 for j in t_range)
                        #v.NPVNoAdapt[gettime(t)-1,m] + sum( [ v.discountfactor[j] * v.NoAdaptCost[j,m] for j in t_range] )
                    end

                    for j in t_range
                        v.NPVNoAdapt[TimestepIndex(j), m] = v.NPVNoAdapt[t, m]
                    end



                    # ** Calculate Protection and Retreat Costs for Each Adaptation Option **
                    lslrPlan_at = p.lslr[TimestepIndex(at_next), m]
                    lslrPlan_atprev = p.lslr[t, m]

                    for i = 1:length(p.adaptoptions)
                        if is_first(t)
                            Rprev = calcHorR(-2, p.adaptoptions[i], p.lslr[ti1, m], view(p.surgeexposure, m, :), p.adaptoptions)
                            v.R[t, m, i] = calcHorR(-2, p.adaptoptions[i], lslrPlan_at, view(p.surgeexposure, m, :), p.adaptoptions)
                        else
                            if p.fixed == false
                                Rprev = v.OptimalR[TimestepIndex(gettime(t) - 1), m]
                                # Assumption: prior protection does not count because it is no longer maintained
                                v.R[t, m, i] = max(v.OptimalR[TimestepIndex(gettime(t) - 1), m], calcHorR(-2, p.adaptoptions[i], lslrPlan_at, view(p.surgeexposure, m, :), p.adaptoptions))
                            else
                                Rprev = v.R[TimestepIndex(convert(Int, p.at[at_index_prev])), m, i]
                                v.R[t, m, i] = calcHorR(-2, p.adaptoptions[i], lslrPlan_at, view(p.surgeexposure, m, :), p.adaptoptions)
                            end

                        end
                        v.SIGMA[t, m, i+1] = (p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, v.R[t, m, i] - p.lslr[t, m]))))
                        v.coastAreaRetreat[t, m, i] = calcCoastArea(view(v.areaparams, m, :), v.R[t, m, i])

                        v.FloodRetreat[t, m, i] = (p.tstep / atstep) * (atstep * v.landvalue[t, m] * 0.04 * calcCoastArea(view(v.areaparams, m, :), v.R[t, m, i]) +
                                                                        max(0, calcCoastArea(view(v.areaparams, m, :), v.R[t, m, i]) - calcCoastArea(view(v.areaparams, m, :), Rprev)) *
                                                                        (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[t, m]) * 1e-4

                        v.RelocateRetreat[t, m, i] = (p.tstep / atstep) *
                                                     max(0, calcCoastArea(view(v.areaparams, m, :), v.R[t, m, i]) - calcCoastArea(view(v.areaparams, m, :), Rprev)) * (p.movefactor * v.ypc_seg[t, m] * 1e-6 * v.popdens_seg[t, m] + p.capmovefactor * p.mobcapfrac * v.capital[t, m] + p.democost * (1 - p.mobcapfrac) * v.capital[t, m]) * 1e-4

                        v.DryLandLossRetreat[t, m, i] = max(0, v.coastAreaRetreat[t, m, i]) # Already takes into account prior adaptation

                        if p.adaptoptions[i] >= 10 || p.adaptoptions[i] == 0

                            if is_first(t)

                                # Hprev = max(p.refA_H[m],calcHorR(-1, p.adaptoptions[i], p.lslr[1,m], p.surgeexposure[m,:], p.adaptoptions))
                                Hprev = calcHorR(-1, p.adaptoptions[i], p.lslr[ti1, m], view(p.surgeexposure, m, :), p.adaptoptions)
                                v.H[t, m, i-1] = calcHorR(-1, p.adaptoptions[i], lslrPlan_at, view(p.surgeexposure, m, :), p.adaptoptions)
                                v.SIGMA[t, m, (i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0, p.lslr[t, m])) / (1.0 + p.psigA[m] * exp(p.psigB[m] * max(0, (v.H[t, m, i-1] - p.lslr[t, m]))))
                                v.FloodProtect[t, m] = 0
                            else
                                if p.fixed == false
                                    Hprev = v.OptimalH[TimestepIndex(gettime(t) - 1), m]
                                    ### Assumption: any prior retreat is credited toward required height, since not starting from original position on coast
                                    lslrPlan_Prot = lslrPlan_at - v.OptimalR[TimestepIndex(gettime(t) - 1), m]
                                    v.H[t, m, i-1] = max(v.OptimalH[TimestepIndex(gettime(t) - 1), m], calcHorR(-1, p.adaptoptions[i], lslrPlan_Prot, view(p.surgeexposure, m, :), p.adaptoptions))
                                    v.SIGMA[t, m, (i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0, p.lslr[t, m])) / (1.0 + p.psigA[m] * exp(p.psigB[m] * max(0, (v.H[t, m, i-1] + v.OptimalR[TimestepIndex(gettime(t) - 1), m] - p.lslr[t, m]))))
                                    v.FloodProtect[t, m] = p.tstep * v.landvalue[t, m] * 0.04 * v.DryLandLossOptimal[TimestepIndex(gettime(t) - 1), m]
                                else
                                    Hprev = v.H[TimestepIndex(convert(Int, p.at[at_index_prev])), m, i-1]
                                    v.H[t, m, i-1] = calcHorR(-1, p.adaptoptions[i], lslrPlan_at, view(p.surgeexposure, m, :), p.adaptoptions)
                                    v.SIGMA[t, m, (i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0, p.lslr[t, m])) / (1.0 + p.psigA[m] * exp(p.psigB[m] * max(0, (v.H[t, m, i-1] - p.lslr[t, m]))))
                                    v.FloodProtect[t, m] = 0
                                end

                            end

                            # Island protection costs are higher
                            if isisland(m, p.xsc)::Int == 1
                                pc = 2 * p.pc0 * p.cci[rgn_ind]
                            else
                                pc = p.pc0 * p.cci[rgn_ind]
                            end

                            v.Construct[t, m, i-1] = (p.tstep / atstep) *
                                                     (p.length[m] * pc * (p.pcfixed + (1 - p.pcfixed) * (v.H[t, m, i-1]^2 - Hprev^2) +
                                                                          p.mc * atstep * v.H[t, m, i-1]) + p.length[m] * 1.7 * v.H[t, m, i-1] * v.landvalue[t, m] * 0.04 / 2 * atstep) * 1e-4

                            ##
                            ## comment out this if block to match the Diaz (2016) GAMS results
                            ## This should NOT be commented out moving forward
                            ##
                            if Hprev >= v.H[t, m, i-1]
                                v.H[t, m, i-1] = Hprev
                                # Just maintenance cost + land value
                                v.Construct[t, m, i-1] = (p.tstep / atstep) * (p.length[m] * pc * p.mc * atstep * v.H[t, m, i-1] + p.length[m] * 1.7 * v.H[t, m, i-1] * v.landvalue[t, m] * 0.04 / 2 * atstep) * 1e-4
                            end

                        end

                        for j in t_range
                            tj = TimestepIndex(j)
                            v.R[tj, m, i] = v.R[t, m, i]
                            v.SIGMA[tj, m, i+1] = (p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, v.R[tj, m, i] - p.lslr[tj, m]))))
                            v.coastAreaRetreat[tj, m, i] = v.coastAreaRetreat[t, m, i]

                            if p.fixed == false && !(is_first(t))
                                v.WetlandRetreat[tj, m] = p.tstep * v.wetlandservice[tj, rgn_ind] * max(v.WetlandLossOptimal[TimestepIndex(gettime(t) - 1), m], v.wetlandloss[TimestepIndex(i), m] * min(v.coastArea[TimestepIndex(i), m], p.wetland[m]))
                            else
                                v.WetlandRetreat[tj, m] = p.tstep * v.wetlandservice[tj, rgn_ind] * v.wetlandloss[tj, m] * min(v.coastArea[tj, m], p.wetland[m])
                            end

                            v.StormCapitalRetreat[tj, m, i] = p.tstep * (1 - v.ρ[tj, rgn_ind]) * v.SIGMA[tj, m, i+1] * v.capital[tj, m]
                            v.StormPopRetreat[tj, m, i] = p.tstep * (1 - v.ρ[tj, rgn_ind]) * v.SIGMA[tj, m, i+1] * v.popdens_seg[tj, m] * v.vsl[tj, m] * p.floodmortality
                            v.StormLossRetreat[tj, m, i] = p.tstep * (1 - v.ρ[tj, rgn_ind]) * v.SIGMA[tj, m, i+1] * v.popdens_seg[tj, m] * p.floodmortality

                            v.FloodRetreat[tj, m, i] = v.FloodRetreat[t, m, i]
                            v.RelocateRetreat[tj, m, i] = v.RelocateRetreat[t, m, i]
                            v.DryLandLossRetreat[tj, m, i] = v.DryLandLossRetreat[t, m, i]

                            # Put all other costs intp $Billions from $M and divide by 10
                            v.StormCapitalRetreat[tj, m, i] = v.StormCapitalRetreat[tj, m, i] * 1e-4
                            v.StormPopRetreat[tj, m, i] = v.StormPopRetreat[tj, m, i] * 1e-4
                            v.WetlandRetreat[tj, m] = v.WetlandRetreat[tj, m] * 1e-4

                            v.RetreatCost[tj, m, i] = v.FloodRetreat[tj, m, i] + v.RelocateRetreat[tj, m, i] + v.StormCapitalRetreat[tj, m, i] + v.StormPopRetreat[tj, m, i] + v.WetlandRetreat[tj, m]

                            if p.adaptoptions[i] >= 10 || p.adaptoptions[i] == 0
                                v.H[tj, m, i-1] = v.H[t, m, i-1]

                                if p.fixed == false && !(is_first(t))
                                    v.SIGMA[tj, m, (i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0, p.lslr[tj, m])) / (1.0 + p.psigA[m] * exp(p.psigB[m] * max(0, (v.H[tj, m, i-1] + v.OptimalR[TimestepIndex(gettime(t) - 1), m] - p.lslr[tj, m]))))
                                    v.FloodProtect[tj, m] = p.tstep * v.landvalue[tj, m] * 0.04 * v.DryLandLossOptimal[TimestepIndex(gettime(t) - 1), m]

                                else
                                    v.SIGMA[tj, m, (i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0, p.lslr[tj, m])) / (1.0 + p.psigA[m] * exp(p.psigB[m] * max(0, (v.H[tj, m, i-1] - p.lslr[tj, m]))))
                                    v.FloodProtect[tj, m] = v.FloodProtect[t, m]
                                end


                                v.WetlandProtect[tj, m] = p.tstep * p.wetland[m] * v.wetlandservice[tj, rgn_ind]

                                v.StormCapitalProtect[tj, m, i-1] = p.tstep * (1 - v.ρ[tj, rgn_ind]) * v.SIGMA[tj, m, (i-1)+7] * v.capital[tj, m]
                                v.StormPopProtect[tj, m, i-1] = p.tstep * (1 - v.ρ[tj, rgn_ind]) * v.SIGMA[tj, m, (i-1)+7] * v.popdens_seg[tj, m] * v.vsl[tj, m] * p.floodmortality
                                v.StormLossProtect[tj, m, i-1] = p.tstep * (1 - v.ρ[tj, rgn_ind]) * v.SIGMA[tj, m, (i-1)+7] * v.popdens_seg[tj, m] * p.floodmortality

                                v.Construct[tj, m, i-1] = v.Construct[t, m, i-1]

                                # Put all other costs into $Billions from $M and divide by 10
                                # Note this is an annual protect cost ($B/year)
                                v.WetlandProtect[tj, m] = v.WetlandProtect[tj, m] * 1e-4
                                v.StormCapitalProtect[tj, m, i-1] = v.StormCapitalProtect[tj, m, i-1] * 1e-4
                                v.StormPopProtect[tj, m, i-1] = v.StormPopProtect[tj, m, i-1] * 1e-4
                                v.FloodProtect[tj, m] = v.FloodProtect[tj, m] * 1e-4

                                v.ProtectCost[tj, m, i-1] = v.Construct[tj, m, i-1] + v.WetlandProtect[tj, m] + v.StormCapitalProtect[tj, m, i-1] + v.StormPopProtect[tj, m, i-1] + v.FloodProtect[tj, m]

                            end

                        end

                        if is_first(t)
                            v.NPVRetreat[t, m, i] = sum(v.discountfactor[TimestepIndex(j)] * v.RetreatCost[TimestepIndex(j), m, i] * 10 for j in t_range)
                        else
                            # Compute NPV Relative to planner's perspective (discounting relative to time t)
                            v.NPVRetreat[t, m, i] = sum(v.discountfactor[TimestepIndex(findind(j, t_range))] * v.RetreatCost[TimestepIndex(j), m, i] * 10 for j in t_range)
                            #v.NPVRetreat[gettime(t)-1,m,i] + sum([v.discountfactor[j] * v.RetreatCost[j,m,i] for j in t_range])
                        end


                        for j in t_range
                            v.NPVRetreat[TimestepIndex(j), m, i] = v.NPVRetreat[t, m, i]
                        end

                        if p.adaptoptions[i] >= 10 || p.adaptoptions[i] == 0
                            if is_first(t)
                                v.NPVProtect[t, m, i-1] = sum([v.discountfactor[TimestepIndex(j)] * v.ProtectCost[TimestepIndex(j), m, i-1] * 10 for j in t_range]) # Protect
                            else
                                # Compute NPV Relative to planner's perspective (discounting relative to time t)
                                v.NPVProtect[t, m, i-1] = sum(v.discountfactor[TimestepIndex(findind(j, t_range))] * v.ProtectCost[TimestepIndex(j), m, i-1] * 10 for j in t_range)
                                #v.NPVProtect[gettime(t)-1,m,i-1] + sum( [ v.discountfactor[j] * v.ProtectCost[j,m,i-1] for j in t_range] ) # Protect
                            end


                            for j in t_range
                                v.NPVProtect[TimestepIndex(j), m, i-1] = v.NPVProtect[t, m, i-1]
                            end
                        end
                    end

                    # ** Choose Least Cost Option **
                    if gettime(t) > 1 && p.fixed
                        # if p.fixed==T and t>1, take first-period choices
                        for j in t_range
                            tj = TimestepIndex(j)
                            v.OptimalProtectLevel[tj, m] = v.OptimalProtectLevel[ti1, m]
                            v.OptimalRetreatLevel[tj, m] = v.OptimalRetreatLevel[ti1, m]
                            v.OptimalOption[tj, m] = v.OptimalOption[ti1, m]
                            v.OptimalLevel[tj, m] = v.OptimalLevel[ti1, m]
                        end
                    else
                        # If p.fixed==F or if p.fixed==T and t==1, calculate optimal level.
                        if p.allowMaintain == true

                            protectInd = findmin(view(v.NPVProtect, TimestepIndex(Int(p.at[at_index])), m, :))[2]
                            retreatInd = findmin(view(v.NPVRetreat, TimestepIndex(Int(p.at[at_index])), m, :))[2]
                        else
                            protDims = size(v.NPVProtect)[3]
                            retDims = size(v.NPVRetreat)[3]
                            protectInd = findmin(v.NPVProtect[TimestepIndex(Int(p.at[at_index])), m, 1:protDims-1])[2]
                            retreatInd = findmin(v.NPVRetreat[TimestepIndex(Int(p.at[at_index])), m, 1:retDims-1])[2]
                        end
                        for j in t_range
                            v.OptimalProtectLevel[TimestepIndex(j), m] = p.adaptoptions[protectInd+1]
                        end

                        if p.noRetreat == true
                            minLevels = Float64[p.adaptoptions[protectInd+1], 0.0]
                            choices = Union{Missing,Float64}[v.NPVProtect[TimestepIndex(Int(p.at[at_index])), m, protectInd], v.NPVNoAdapt[TimestepIndex(Int(p.at[at_index])), m]]

                            leastcost = -1 * findmin(choices)[2]
                            if leastcost == -2
                                leastcost = -3 # Account for retreat being removed from choice set
                            end
                            leastlevel = minLevels[findmin(choices)[2]]
                        else
                            for j in t_range
                                v.OptimalRetreatLevel[TimestepIndex(j), m] = p.adaptoptions[retreatInd]
                            end
                            minLevels = Float64[p.adaptoptions[protectInd+1], p.adaptoptions[retreatInd], 0.0]

                            choices = Union{Missing,Float64}[v.NPVProtect[TimestepIndex(Int(p.at[at_index])), m, protectInd], v.NPVRetreat[TimestepIndex(Int(p.at[at_index])), m, retreatInd], v.NPVNoAdapt[TimestepIndex(Int(p.at[at_index])), m]]
                            leastcost = -1 * findmin(choices)[2]
                            leastlevel = minLevels[findmin(choices)[2]]
                        end
                        for j in t_range
                            tj = TimestepIndex(j)
                            v.OptimalOption[tj, m] = leastcost
                            v.OptimalLevel[tj, m] = leastlevel
                        end
                    end

                    # Assign costs to optimal variables
                    if v.OptimalOption[t, m] == -1

                        # Protect Cost
                        protInd = findall(i -> i == v.OptimalLevel[t, m], p.adaptoptions)[1] - 1
                        for j in t_range
                            tj = TimestepIndex(j)
                            v.OptimalCost[tj, m] = v.ProtectCost[tj, m, protInd]
                            # Assign Subcosts
                            v.OptimalStormCapital[tj, m] = v.StormCapitalProtect[tj, m, protInd]
                            v.OptimalStormPop[tj, m] = v.StormPopProtect[tj, m, protInd]
                            v.OptimalConstruct[tj, m] = v.Construct[tj, m, protInd]
                            v.OptimalWetland[tj, m] = v.WetlandProtect[tj, m]
                            v.OptimalRelocate[tj, m] = 0
                            v.OptimalFlood[tj, m] = v.FloodProtect[tj, m]

                            # Assign Alternative Metrics
                            # Assume once seawall is built, wetland area is permanently destroyed
                            v.WetlandLossOptimal[tj, m] = p.wetland[m]
                        end

                        if gettime(t) == 1
                            for i in t_range
                                ti = TimestepIndex(i)
                                v.DryLandLossOptimal[ti, m] = 0
                                v.OptimalH[ti, m] = max(0, v.H[ti, m, protInd])
                                v.OptimalR[ti, m] = 0
                                if i == 1
                                    v.StormLossOptimal[ti, m] = v.StormLossProtect[ti, m, protInd]

                                else
                                    v.StormLossOptimal[ti, m] = v.StormLossOptimal[TimestepIndex(i - 1), m] + v.StormLossProtect[ti, m, protInd]

                                end
                            end
                        else
                            for j in t_range
                                v.OptimalR[TimestepIndex(j), m] = max(0, v.OptimalR[TimestepIndex(gettime(t) - 1), m])
                            end
                            for i in t_range
                                ti = TimestepIndex(i)
                                v.OptimalH[ti, m] = max(v.H[ti, m, protInd], v.OptimalH[TimestepIndex(gettime(t) - 1), m])

                                v.DryLandLossOptimal[ti, m] = max(0, v.DryLandLossOptimal[TimestepIndex(i - 1), m])
                                v.StormLossOptimal[ti, m] = v.StormLossOptimal[TimestepIndex(i - 1), m] + v.StormLossProtect[ti, m, protInd]
                            end
                        end

                    elseif v.OptimalOption[t, m] == -2
                        # Retreat Cost

                        retInd = findall(i -> i == v.OptimalLevel[t, m], p.adaptoptions)[1]
                        for j in t_range
                            tj = TimestepIndex(j)
                            v.OptimalCost[tj, m] = v.RetreatCost[tj, m, retInd]
                            # Assign Subcosts
                            v.OptimalStormCapital[tj, m] = v.StormCapitalRetreat[tj, m, retInd]
                            v.OptimalStormPop[tj, m] = v.StormPopRetreat[tj, m, retInd]
                            v.OptimalConstruct[tj, m] = 0
                            v.OptimalWetland[tj, m] = v.WetlandRetreat[tj, m]
                            v.OptimalFlood[tj, m] = v.FloodRetreat[tj, m, retInd]
                            v.OptimalRelocate[tj, m] = v.RelocateRetreat[tj, m, retInd]
                        end

                        if is_first(t)

                            for j in t_range
                                tj = TimestepIndex(j)
                                v.DryLandLossOptimal[tj, m] = v.DryLandLossRetreat[tj, m, retInd]
                                v.OptimalH[tj, m] = 0
                            end

                            for i in t_range
                                ti = TimestepIndex(i)
                                v.OptimalR[ti, m] = max(0, v.R[ti, m, retInd])
                                v.WetlandLossOptimal[ti, m] = v.wetlandloss[ti, m] * min(v.coastArea[ti, m], p.wetland[m])
                                if i == 1
                                    v.StormLossOptimal[ti, m] = v.StormLossRetreat[ti, m, findall(k -> k == v.OptimalLevel[t, m], p.adaptoptions)[1]]
                                else
                                    v.StormLossOptimal[ti, m] = v.StormLossOptimal[TimestepIndex(i - 1), m] + v.StormLossRetreat[ti, m, retInd]
                                end
                            end
                        else
                            for j in t_range
                                v.OptimalH[TimestepIndex(j), m] = v.OptimalH[TimestepIndex(gettime(t) - 1), m]
                            end

                            for i in t_range
                                ti = TimestepIndex(i)
                                v.OptimalR[ti, m] = max(v.R[ti, m, retInd], v.OptimalR[TimestepIndex(gettime(t) - 1), m])
                                # Cumulative total wetland area lost; if protected previously, all wetland is lost
                                v.WetlandLossOptimal[ti, m] = max(v.WetlandLossOptimal[TimestepIndex(gettime(t) - 1), m], v.wetlandloss[ti, m] * min(v.coastArea[ti, m], p.wetland[m]))
                                v.StormLossOptimal[ti, m] = v.StormLossOptimal[TimestepIndex(i - 1), m] + p.tstep * (1 - v.ρ[ti, rgn_ind]) * v.popdens_seg[ti, m] * p.floodmortality * v.SIGMA[ti, m, retInd]
                                v.DryLandLossOptimal[ti, m] = max(v.DryLandLossOptimal[TimestepIndex(i - 1), m], v.DryLandLossRetreat[ti, m, retInd])
                            end

                        end

                    else
                        # No Adaptation

                        for j in t_range
                            tj = TimestepIndex(j)
                            v.OptimalCost[tj, m] = v.NoAdaptCost[tj, m]
                            # Assign Subcosts
                            v.OptimalStormCapital[tj, m] = v.StormCapitalNoAdapt[tj, m]
                            v.OptimalStormPop[tj, m] = v.StormPopNoAdapt[tj, m]
                            v.OptimalConstruct[tj, m] = 0
                            v.OptimalWetland[tj, m] = v.WetlandNoAdapt[tj, m]
                            v.OptimalFlood[tj, m] = v.FloodNoAdapt[tj, m]
                            v.OptimalRelocate[tj, m] = v.RelocateNoAdapt[tj, m]
                        end

                        if is_first(t)
                            for j in t_range
                                tj = TimestepIndex(j)
                                v.OptimalH[tj, m] = 0
                                v.DryLandLossOptimal[tj, m] = v.DryLandLossNoAdapt[tj, m]
                            end
                            v.OptimalR[t, m] = max(0, p.lslr[t, m])

                            for i in t_range
                                ti = TimestepIndex(i)
                                if i > 1
                                    v.OptimalR[ti, m] = max(v.OptimalR[t, m], v.OptimalR[TimestepIndex(i - 1), m], p.lslr[ti, m])
                                else
                                    v.OptimalR[ti, m] = max(v.OptimalR[t, m], p.lslr[ti, m])
                                end
                                v.WetlandLossOptimal[ti, m] = v.wetlandloss[ti, m] * min(v.coastArea[ti, m], p.wetland[m])
                                if i == 1
                                    v.StormLossOptimal[ti, m] = v.StormLossNoAdapt[ti, m]
                                else
                                    v.StormLossOptimal[ti, m] = v.StormLossOptimal[TimestepIndex(i - 1), m] + v.StormLossNoAdapt[ti, m]
                                end
                            end
                        else
                            for j in t_range
                                v.OptimalH[TimestepIndex(j), m] = v.OptimalH[TimestepIndex(gettime(t) - 1), m]
                            end

                            for i in t_range
                                ti = TimestepIndex(i)
                                v.OptimalR[ti, m] = max(v.OptimalR[TimestepIndex(gettime(t) - 1), m], v.OptimalR[TimestepIndex(i - 1), m], p.lslr[ti, m])

                                v.WetlandLossOptimal[ti, m] = max(v.WetlandLossOptimal[TimestepIndex(gettime(t) - 1), m], v.wetlandloss[ti, m] * min(v.coastArea[ti, m], p.wetland[m]))
                                v.StormLossOptimal[ti, m] = v.StormLossOptimal[TimestepIndex(i - 1), m] + v.StormLossNoAdapt[ti, m]
                                v.DryLandLossOptimal[ti, m] = max(v.DryLandLossOptimal[TimestepIndex(i - 1), m], v.DryLandLossNoAdapt[ti, m])
                            end

                        end
                    end

                    if last == 1
                        v.NPVOptimal[m] = sum([v.discountfactor[TimestepIndex(j)] * v.OptimalCost[TimestepIndex(j), m] * 10 for j = 1:p.ntsteps])

                    end
                end
            end  # end segment loop

            if last == 1
                v.NPVOptimalTotal = sum(v.NPVOptimal)
            end
        end # end if t in adaptation period statement
    end # end function
end # end module
