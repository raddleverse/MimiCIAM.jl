# Catherine Ledna
# February 3, 2020
##------------------------------------------------------------------------
# CIAM Model
#------------------------------------------------------------------------
# Implements CIAM model adapted from Diaz, 2016. 
#------------------------------------------------------------------------

using Mimi

@defcomp slrcost begin
    # Define all variables, parameters and indices used by this module 
    # --- Indices ---
    regions = Index()
    segments = Index()                      
    adaptPers = Index()

    # --- Region / segment mapping ---
    segID = Parameter( index = [segments])   # Unique segment numeric identifier 
    #xsc = Parameter{Dict{Int32, Tuple{Int32, Bool, Bool}}}()  # Region to segment mapping (dictionary) to keep track of which segments belong to each region    # TWmod
    xsc = Parameter{Dict{Any, Any}}()  # Region to segment mapping (dictionary) to keep track of which segments belong to each region    # TWmod
    rcp = Parameter{Int}()                   # RCP being run (metadata; not used in run)
    percentile = Parameter{Int}()            # Percentile of RCP being run (metadata; not used in run)
    ssp = Parameter{Int}()                   # SSP being used (0 for base case)

    # ---Time-related Parameters---
    tstep = Parameter()                     # Length of individual time-step (years)
    at = Parameter( index = [adaptPers])    # Array of time indices that mark starts of adaptation periods 
    ntsteps = Parameter{Int}()              # Number of time-steps     

    # ---Model Parameters ---
    fixed = Parameter{Bool}()               # Run model as fixed (T) or flexible (F) with respect to adaptation
    noRetreat = Parameter{Bool}()           # Default (F). If T, segments will either protect or not adapt.
    allowMaintain = Parameter{Bool}()       # Default T. If T, segments will have the option to maintain current defenses 

    # ---Socioeconomic Parameters---
    popinput=Parameter()                               # Input for population data source: 0 (default), 1 (Jones & O'Neill, 2016), 2 (Merkens et al, 2016)
    pop = Parameter(index = [time, regions])           # Population of region (million people) (from MERGE or SSPs)
    refpopdens = Parameter( index = [regions])         # Reference population density of region (people / km^2)
    rgn_ind_usa = Parameter{Int}()                     # Lookup parameter for USA region index, used in refpopdens and ypc  
                                                       #    for USA benchmark in vsl, rho and fundland calculations
    popdens = Parameter( index = [segments])           # Pop density of segment in time t = 1 (people/km^2)
    ypcc = Parameter(index = [time, regions])          # GDP per capita per region ($2010 per capita)
     
    popdens_seg = Variable(index = [time, segments])          # Population density of segment extrapolated forward in time (people / km^2)    
    popdens_seg_jones = Parameter(index=[time,segments])      # Holder for Jones and O'Neill population density
    popdens_seg_merkens = Parameter(index=[time,segments])
    ypc_seg = Variable(index = [time, segments])              # GDP per capita by segment ($2010 per capita) (multiplied by scaling factor)
    refA_R = Parameter(index = [segments])                # Reference retreat level of adaptation in 0 period 
    refA_H = Parameter(index=[segments])                # Reference height for adaptation in 0 period

    # ---Land Parameters---  
    landinput = Parameter{Bool}()                   # Set to T for FUND or F for GTAP
         
    gtapland = Parameter( index = [regions])        # GTAP land value in 2007 (million 2010$ / km^2)
    dvbm = Parameter()                              # FUND value of OECD dryland per Darwin et al 1995 converted from $1995 ($2010M per sqkm) (5.376)
    kgdp = Parameter()                              # Capital output ratio (per MERGE) (3 by default) 
    discountrate = Parameter()                      # Discount rate (0.04 by default)
    depr = Parameter()                              # Fraction of capital that has not been depreciated over adaptation period (retreat cases)
    
    
    landdata = Variable( index = [regions])         # Takes on value of either fundland or gtapland
    fundland = Variable( index = [regions])         # FUND land value in 1995 (calculated in run_timestep) (million 2010$ / km^2), 
                                                    #   Q maybe import directly? 

    rgn_ind_canada = Parameter{Int}()               # Region index for Canada (Used as reference for Greenland land appreciation)
    land_appr = Variable(index = [time, regions])   # Land appreciation rate (calculated as regression by Yohe ref Abraham and Hendershott) 
    coastland = Variable(index = [time, segments])  # Coastal land value (function of interior land value * scaling factor) ($2010M per sqkm)
    landvalue = Variable(index = [time, segments])  # Total endowment value of land ($2010M per sqkm)
    landrent = Variable(index=[time,segments])      # Annual rental value of land ($2010M/sqkm/year) 


    ρ = Variable(index = [time, regions])           # Country-wide resilience parameter (logistic function related to GDP)
    capital = Variable(index = [time, segments])    # Total endowment value of capital stock (million $2010 / km^2)
    discountfactor = Variable(index=[time])         # Discount factor (derived from discount rate)
    
    # ---Coastal Parameters---
    length = Parameter(index = [segments])          # Segment length (km)

    # ---Protection Parameters---
    cci = Parameter(index = [regions])
    pcfixed = Parameter()                   # Fraction of protection cost that is fixed (not variable in height) (0.3)
    mc = Parameter()                        # Maintenance cost (Hillen et al, 2010) (2%/yr) 
    pc0 = Parameter()                       # Reference cost of protection (million 2010$ / km / vert m^2) (6.02 by default)
    

    # ---Retreat / No Adapt Parameters---
    mobcapfrac = Parameter()                # Fraction of capital that is mobile (0.25)
    movefactor = Parameter()                # Cost to relocate people as a factor of annual income (Tol 3x RMendelsohn 0.5x) (1)
    capmovefactor = Parameter()             # Cost to relocate mobile capital as a fraction of asset value (0.1)
    democost = Parameter()                  # Cost to demolish immobile capital as fraction of asset (0.05)

    # # ---Surge Exposure Parameters---
    # Protection case
    psig0 = Parameter( index = [segments])
    psig0coef = Parameter(index = [segments])
    psigA = Parameter(index = [segments])                       # psigA in GAMS code
    psigB = Parameter(index = [segments])                       # psigB in GAMS code

    # Retreat / No Adapt Cases
    rsig0 = Parameter( index = [segments])
    rsigA = Parameter( index = [segments])                       # rsigA in GAMS code
    rsigB = Parameter( index = [segments])                       # rsigB in GAMS code
    

    # ---Storm damage parameters---
    floodmortality = Parameter()                # Flood deaths as percent of exposed population; (Jonkman Vrijling 2008) (0.01) 
    vslel = Parameter()                         # Elasticity of vsl (0.5)
    vslmult = Parameter()                       # multiplier on USA GDP (216)
    vsl = Variable(index = [time, segments])     # Value of statistical life (million 2010$)
                                           
    # ---Wetland Loss Parameters---
    wvbm = Parameter()                                      # Annual value of wetland services (million 2010$ / km^2 / yr); (Brander et al 2006)  (0.376) 
    wetland = Parameter( index = [segments])                # Initial wetland area in coastal segment (km^2)
    wmaxrate = Parameter()                                  # Maximum rate of wetland accretion (m per yr) per Kirwan et al 2010 (0.01)
    wvel = Parameter()                                      # income elasticity of wetland value (1.16) (Brander et al, 2006)
    wvpdl = Parameter()                                     # Population density elasticity of wetland value (0.47) (Brander et al, 2006)
    
    wetlandservice = Variable(index = [time, regions])      # Annual value of wetland services adjusted for income and density (Brander et al 2006) ($2010M/km^2/year)
    wetlandloss = Variable(index = [time, segments])        # Fractional loss of wetland due to slr


    # ---Sea Level Rise Parameters---
    lslr = Parameter(index = [time, segments])                # Local sea level rise (m) 

    adaptoptions = Parameter(index = [6])                     # Index of available adaptation levels for protect and retreat (0 is no adaptation)
    surgeexposure = Parameter{Float64}( index = [segments, 5])# Storm surge exposure levels (corresponding to each designated adaptation option)
    

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

    coastArea = Variable(index=[time, segments])            # Coast area inundated (km^2)
    
    
    # ---Intermediate Variables---
    WetlandNoAdapt = Variable(index = [time, segments])
    FloodNoAdapt = Variable(index = [time, segments])
    StormCapitalNoAdapt = Variable(index = [time, segments])
    StormPopNoAdapt = Variable(index = [time, segments])
    RelocateNoAdapt = Variable(index = [time, segments])
    StormLossNoAdapt = Variable(index = [time,segments])
    DryLandLossNoAdapt = Variable(index=[time,segments])

    Construct = Variable(index = [time, segments, 5])
    WetlandProtect = Variable(index = [time, segments])    
    StormCapitalProtect = Variable(index = [time, segments,5])
    StormPopProtect = Variable(index = [time, segments,5])
    StormLossProtect = Variable(index = [time,segments,5])
    FloodProtect = Variable(index = [time,segments])
    
    WetlandRetreat = Variable(index = [time, segments])
    StormCapitalRetreat = Variable(index = [time, segments,6])
    StormPopRetreat = Variable(index = [time, segments,6])
    StormLossRetreat = Variable(index = [time,segments,6])
    FloodRetreat = Variable(index = [time, segments, 6])
    RelocateRetreat = Variable(index = [time, segments, 6])
    DryLandLossRetreat = Variable(index=[time,segments,6])
    coastAreaRetreat = Variable(index = [time,segments,6])
    coastAreaNoAdapt=Variable(index=[time,segments])
    
    # --- Decision Variables --- (evaluated brute force)
    H = Variable(index = [time, segments, 5])       # Height of current sea wall, no retreat (m)
    R = Variable(index = [time, segments, 6])       # Retreat perimeter (m)
    SIGMA = Variable(index = [time, segments, 12])  # Expected value of effective exposure area for over-topping surge (all cases)
                                                    # Order of sigma values: 1 no adapt case, 6 retreat cases, 5 protect cases in ascending order

    # ---Outcome Variables---   
    OptimalH = Variable(index=[time,segments])               # m; Holder to track height built across timesteps (cumulative)
    OptimalR = Variable(index=[time,segments])               # m; Holder to track retreat radius across timesteps (cumulative)
    WetlandLossOptimal = Variable( index = [time,segments])  # km2; Cumulative wetland loss from optimal decision
    DryLandLossOptimal = Variable( index = [time, segments]) # km2; Cumulative loss of dry land from optimal decision 

   # DrylandLost = Variable(index=[time,segments])            # km2; container to track cumulative lost dryland 
    WetlandLost = Variable(index=[time,segments])            # km2; container to track cumulative lost wetland 

    NoAdaptCost = Variable(index = [time, segments])         # Cost of not adapting (e.g. reactive retreat) (2010$)
    ProtectCost = Variable(index = [time, segments, 5])      # Total cost of protection at each level      
    RetreatCost = Variable(index = [time, segments, 6])      # Total cost of retreat at each level   
    OptimalRetreatLevel = Variable(index = [time, segments])
    OptimalProtectLevel = Variable(index = [time, segments])
    OptimalCost = Variable(index = [time, segments])          # Optimal cost based on NPV relative to start of adaptation period 
    OptimalLevel = Variable(index = [time, segments])         # Fixed optimal level (1,10,100,1000,10000)
    OptimalOption = Variable(index = [time, segments])        # Fixed adaptation decision (-1 - protect, -2 - retreat, -3 - no adapt) 
    NPVRetreat = Variable(index = [time,segments, 6])        
    NPVProtect = Variable(index = [time,segments,  5])
    NPVNoAdapt = Variable(index = [time,segments])
    NPVOptimal = Variable(index = [segments])               # NPV of cost of optimal decisions relative to t=1
    NPVOptimalTotal = Variable()                            # Total NPV of all segments from optimal decision 
    StormLossOptimal = Variable( index = [time, segments])  # Cumulative expected loss of life (num people) from storm surges from optimal decision
  
    # ---Subcategories of Optimal Choice----
    OptimalStormCapital = Variable(index = [time, segments])
    OptimalStormPop = Variable(index = [time, segments])
    OptimalConstruct = Variable(index = [time, segments])
    OptimalWetland = Variable(index = [time, segments])
    OptimalFlood = Variable(index = [time, segments])
    OptimalRelocate = Variable(index = [time, segments])
    

    function run_timestep(p, v, d, t)   
        println(gettime(t))
        # In first period, initialize all non-adaptation dependent intermediate variables for all timesteps
        if is_first(t)
          #  1. Initialize non-region dependent intermediate variables
            for i in collect(1:Int(p.ntsteps)) 
               v.discountfactor[TimestepIndex(i)] = 1/(1 + p.discountrate)^(p.tstep * (i-1)) # TWmod
            end
    
            # 2. Initialize region-dependent intermediate variables
            for r in d.regions
                # Determine land input value (true = FUND, false = GTAP)
                if p.landinput 
                    v.fundland[r] = min(p.dvbm, max(0.005, p.dvbm * p.ypcc[t,r] * p.refpopdens[r] / (p.ypcc[TimestepIndex(1),p.rgn_ind_usa] * p.refpopdens[p.rgn_ind_usa])))
                    v.landdata[r] = v.fundland[r]
                else
                    v.landdata[r] = p.gtapland[r]
                end
    
                # Calculate regional wetland service, resilience (rho), and land appreciation variables for the first period and 
                #   subsequent periods 
                v.wetlandservice[t,r] = p.wvbm * ((p.ypcc[t,r] / p.ypcc[TimestepIndex(1),p.rgn_ind_usa])^p.wvel * (p.refpopdens[r] /27.59)^p.wvpdl) 
                v.ρ[t,r] = p.ypcc[t,r] / (p.ypcc[t,r] + p.ypcc[TimestepIndex(1),p.rgn_ind_usa])
                v.land_appr[t,r] = 1.
    
                for i in collect(2:Int(p.ntsteps))
                    v.land_appr[TimestepIndex(i),r] = v.land_appr[TimestepIndex(i-1),r] * exp(0.565 * growthrate(p.ypcc[TimestepIndex(i-1),r], p.ypcc[TimestepIndex(i),r]) + 0.313 * growthrate(p.pop[TimestepIndex(i-1),r], p.pop[TimestepIndex(i),r]))
                    v.wetlandservice[TimestepIndex(i),r] = v.land_appr[TimestepIndex(i),r] * v.wetlandservice[TimestepIndex(1),r]
                    v.ρ[TimestepIndex(i),r] = p.ypcc[TimestepIndex(i),r] / (p.ypcc[TimestepIndex(i),r] + p.ypcc[TimestepIndex(1),p.rgn_ind_usa]) 
                end    
            end
    
            # 3. Initialize segment-dependent variables 
            for m in d.segments
                rgn_ind = getregion(m, p.xsc) # Identify the region the segment belongs to
                
                # Initialize first-period population density, coast area and surge parameters
                if p.popinput==0    
                    v.popdens_seg[t,m] = p.popdens[m]
                elseif p.popinput==1
                    v.popdens_seg[t,m]=p.popdens_seg_jones[TimestepIndex(1),m]
                elseif p.popinput==2
                    v.popdens_seg[t,m]=p.popdens_seg_merkens[TimestepIndex(1),m]
                end
                v.areaparams[m,:] = [p.area1[m] p.area2[m] p.area3[m] p.area4[m] p.area5[m] p.area6[m] p.area7[m] p.area8[m] p.area9[m] p.area10[m] p.area11[m] p.area12[m] p.area13[m] p.area14[m] p.area15[m]] 

                # Greenland segments are treated differently 
                if isgreenland(m,p.xsc)==1
                    v.ypc_seg[t,m] =22642*1.01^1   # FLAG: assumes t is an index (1-20)
                    v.vsl[t,m] = 1e-6 * p.vslmult * p.ypcc[TimestepIndex(1),p.rgn_ind_usa] * (v.ypc_seg[t,m]/p.ypcc[TimestepIndex(1),p.rgn_ind_usa])^p.vslel
                    v.coastland[t,m] = (v.land_appr[TimestepIndex(1),p.rgn_ind_canada] * v.landdata[p.rgn_ind_canada]) * max(0.5, log(1+v.popdens_seg[t,m])/log(25))
                    v.landvalue[t,m] = min(v.coastland[t,m], (v.land_appr[TimestepIndex(1),p.rgn_ind_canada] * v.landdata[p.rgn_ind_canada]))
                else
                    v.ypc_seg[t,m] = p.ypcc[t,rgn_ind] * max(0.9, (v.popdens_seg[TimestepIndex(1),m]/250.)^0.05)
                    v.vsl[t,m] = 1e-6 * p.vslmult * p.ypcc[TimestepIndex(1),p.rgn_ind_usa] * (p.ypcc[t,rgn_ind]/p.ypcc[TimestepIndex(1),p.rgn_ind_usa])^p.vslel
                    v.coastland[t,m] = max(0.5, log(1+v.popdens_seg[t,m])/log(25)) * (v.land_appr[t,rgn_ind] * v.landdata[rgn_ind])  # Interior * scaling factor
                    v.landvalue[t,m] = min(v.coastland[t,m], (v.land_appr[t,rgn_ind] * v.landdata[rgn_ind]))
                end
    
                v.capital[t,m] = p.kgdp * v.ypc_seg[t,m] * v.popdens_seg[t,m] * 1e-6
                v.coastArea[t,m] = calcCoastArea(v.areaparams[m,:], p.lslr[t,m])  
                
                for i in 2:Int(p.ntsteps)  
                    if p.popinput==0      
                        v.popdens_seg[TimestepIndex(i),m] = v.popdens_seg[TimestepIndex(i-1),m] * (1 + growthrate(p.pop[TimestepIndex(i-1),rgn_ind], p.pop[TimestepIndex(i),rgn_ind])) 
                    elseif p.popinput==1
                        v.popdens_seg[TimestepIndex(i),m]=p.popdens_seg_jones[TimestepIndex(i),m]
                    elseif p.popinput==2
                        v.popdens_seg[TimestepIndex(i),m]=p.popdens_seg_merkens[TimestepIndex(i),m]
                    end
     
                    # Special treatment for Greenland segments 
                    if isgreenland(m,p.xsc)==1
                        v.ypc_seg[TimestepIndex(i),m] =22642*1.01^i   # FLAG: assumes i is an index (1-20)
                        v.vsl[TimestepIndex(i),m] = 1e-6 * p.vslmult * p.ypcc[TimestepIndex(i),p.rgn_ind_usa] * (v.ypc_seg[TimestepIndex(i),m]/p.ypcc[TimestepIndex(i),p.rgn_ind_usa])^p.vslel 
                        v.coastland[TimestepIndex(i),m] = (v.land_appr[TimestepIndex(i),p.rgn_ind_canada] * v.landdata[p.rgn_ind_canada]) * max(0.5, log(1+v.popdens_seg[TimestepIndex(i),m])/log(25))
                        v.landvalue[TimestepIndex(i),m] = min(v.coastland[TimestepIndex(i),m], (v.land_appr[TimestepIndex(i),p.rgn_ind_canada] * v.landdata[p.rgn_ind_canada]))
    
                    else
                        v.ypc_seg[TimestepIndex(i),m] = p.ypcc[TimestepIndex(i),rgn_ind] * max(0.9, (v.popdens_seg[TimestepIndex(1),m]/250.)^0.05) # ypcc * popdens scaling factor
                        v.coastland[TimestepIndex(i),m] = max(0.5, log(1+v.popdens_seg[TimestepIndex(i),m])/log(25)) * (v.land_appr[TimestepIndex(i),rgn_ind] * v.landdata[rgn_ind])
                        v.vsl[TimestepIndex(i),m] = 1e-6 * p.vslmult * p.ypcc[TimestepIndex(i),p.rgn_ind_usa] * (p.ypcc[TimestepIndex(i),rgn_ind]/p.ypcc[TimestepIndex(i),p.rgn_ind_usa])^p.vslel  
                        v.landvalue[TimestepIndex(i),m] = min(v.coastland[TimestepIndex(i),m], (v.land_appr[TimestepIndex(i),rgn_ind] * v.landdata[rgn_ind]))
     
                    end
    
                    v.capital[TimestepIndex(i),m] = p.kgdp * v.ypc_seg[TimestepIndex(i),m] * v.popdens_seg[TimestepIndex(i),m] * 1e-6 
                    v.coastArea[TimestepIndex(i),m] = calcCoastArea(v.areaparams[m,:], p.lslr[TimestepIndex(i),m])
                    v.wetlandloss[TimestepIndex(i-1),m] = min(1, (localrate(p.lslr[TimestepIndex(i-1),m], p.lslr[TimestepIndex(i),m], p.tstep)/p.wmaxrate)^2)
                    
    
                end
                v.wetlandloss[TimestepIndex(p.ntsteps),m] = min(1, (localrate(p.lslr[TimestepIndex(p.ntsteps-1),m], p.lslr[TimestepIndex(p.ntsteps),m], p.tstep)/p.wmaxrate)^2)  
    
            end
    
        end
    
        if (gettime(t) in p.at)           
            adapt_range = collect(1:length(p.adaptoptions))         
    
            # Determine length of adaptation period ("atstep")
            g(c) = c == gettime(t)
            at_index = findall(g, p.at)[1] # Find index corresponding to adaptation period in p.at
            at_index_next = at_index + 1    # Find index corresponding to next adaptation period
            at_index_prev = max(1,at_index - 1)
                        
            at_prev = Int(p.at[at_index_prev])      # TODO pick one - either index into vector or use the period, it's confusing
    
            if at_index_next <= length(p.at)
                atstep = (p.at[at_index_next] - p.at[at_index])*p.tstep   # In years
                at_next = Int(p.at[at_index_next])  
                last_t = at_next-1
                last = 0
            else
                # Deal with special case of last adaptation period
                atstep = p.tstep*p.ntsteps - ((p.at[at_index]-1) * p.tstep) 
                at_next = p.ntsteps # Flag this assumes timesteps are indices not years (1:20 not 2010:2100)
                last_t = at_next
                last = 1
            end
            t_range = collect(gettime(t):last_t)
#println(typeof(t_range), " - ", t_range, " - ", TimestepIndex(t_range))
            
            for m in d.segments
                if atstep==0
                else
                    rgn_ind = getregion(m, p.xsc)
     
                    # ** Calculate No Adaptation Costs **
                    for i in t_range
                        R_NoAdapt = max(0, p.lslr[TimestepIndex(i),m])

                        # For initial state in SLR cases, make adaptation decision relative to baseline (refA_H or R)
                        if p.rcp>0
                            R_NoAdapt = max(R_NoAdapt, p.refA_H[m],p.refA_R[m])
                        end

                        # Incorporate any previous period adaptation
                        if p.fixed==false && !(is_first(t))
                            
                            R_NoAdapt = max(R_NoAdapt, v.OptimalR[TimestepIndex(gettime(t)-1),m])
                            v.WetlandNoAdapt[TimestepIndex(i),m] = p.tstep * v.wetlandservice[TimestepIndex(i),rgn_ind] * max(v.WetlandLossOptimal[TimestepIndex(gettime(t)-1),m],v.wetlandloss[TimestepIndex(i),m] * min(v.coastArea[TimestepIndex(i),m], p.wetland[m]))
                            if i==gettime(t)   
                                # For start of new adaptation period, take into account (lack of) retreat done in previous periods (i.e. if they protected instead)
                                # This results in double-costs for this period b/c no adaptation is set up to compute relative to t+1 lslr 
                                v.coastAreaNoAdapt[TimestepIndex(i),m] = calcCoastArea(v.areaparams[m,:],v.OptimalR[TimestepIndex(gettime(t)-1),m])
                            else
                                v.coastAreaNoAdapt[TimestepIndex(i),m] = calcCoastArea(v.areaparams[m,:],R_NoAdapt)
                            end
                            
                        else
                            v.coastAreaNoAdapt[TimestepIndex(i),m]= v.coastArea[TimestepIndex(i),m]
                            v.WetlandNoAdapt[TimestepIndex(i),m] = p.tstep * v.wetlandservice[TimestepIndex(i),rgn_ind] * v.wetlandloss[TimestepIndex(i),m] * min(v.coastArea[TimestepIndex(i),m], p.wetland[m]) 
                        end
                        
                        
                        # Storm Costs 
                        v.SIGMA[TimestepIndex(i),m,1] = p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, R_NoAdapt - p.lslr[TimestepIndex(i),m]))) # expected value of exposure area 
                        v.StormCapitalNoAdapt[TimestepIndex(i),m] = p.tstep * (1 - v.ρ[TimestepIndex(i),rgn_ind ]) * v.SIGMA[TimestepIndex(i),m,1] * v.capital[TimestepIndex(i),m]
                        v.StormPopNoAdapt[TimestepIndex(i),m] = p.tstep * (1 - v.ρ[TimestepIndex(i),rgn_ind ]) * v.popdens_seg[TimestepIndex(i),m] * v.vsl[TimestepIndex(i),m] * p.floodmortality * v.SIGMA[TimestepIndex(i),m,1] 
                        
                         
                        v.StormLossNoAdapt[TimestepIndex(i),m] = p.tstep * (1 - v.ρ[TimestepIndex(i),rgn_ind ]) * v.popdens_seg[TimestepIndex(i),m] * p.floodmortality * v.SIGMA[TimestepIndex(i),m,1]
                        if i==p.ntsteps
                            v.DryLandLossNoAdapt[TimestepIndex(i),m] = max(0,v.coastAreaNoAdapt[TimestepIndex(i),m]) # km^2
                        else
                            # In case of negative or decreasing slr, can assume that previous inundated area is reclaimed 
                            v.DryLandLossNoAdapt[TimestepIndex(i),m] = max(0,v.coastAreaNoAdapt[TimestepIndex(i),m],v.coastArea[TimestepIndex(i+1),m]) # includes future period loss and previous adaptation if applicable
                        end

                        # Flood and relocation costs 
                        if i==p.ntsteps
                            v.FloodNoAdapt[TimestepIndex(i),m] = v.FloodNoAdapt[TimestepIndex(i-1),m]
                            v.RelocateNoAdapt[TimestepIndex(i),m] = v.RelocateNoAdapt[TimestepIndex(i-1),m]
                        else
                            v.FloodNoAdapt[TimestepIndex(i),m]  = p.tstep * v.landvalue[TimestepIndex(i),m]*.04 * v.DryLandLossNoAdapt[TimestepIndex(i),m] + max(0,v.DryLandLossNoAdapt[TimestepIndex(i),m] - v.coastAreaNoAdapt[TimestepIndex(i),m]) * 
                                (1 - p.mobcapfrac) * v.capital[TimestepIndex(i),m]                    
                            
                            v.RelocateNoAdapt[TimestepIndex(i),m] = max(0,v.DryLandLossNoAdapt[TimestepIndex(i),m] - v.coastAreaNoAdapt[TimestepIndex(i),m]) * (5 * p.movefactor * v.ypc_seg[TimestepIndex(i),m]*1e-6*v.popdens_seg[TimestepIndex(i),m] +
                                p.capmovefactor * p.mobcapfrac * v.capital[TimestepIndex(i),m] + p.democost * (1 - p.mobcapfrac) * v.capital[TimestepIndex(i),m])
                        end
                                
                        # Put all costs into $Billions and divide by 10
                        v.WetlandNoAdapt[TimestepIndex(i),m] = v.WetlandNoAdapt[TimestepIndex(i),m] * 1e-4
                        if i<p.ntsteps # already occurred in previous timestep
                            v.FloodNoAdapt[TimestepIndex(i),m] = v.FloodNoAdapt[TimestepIndex(i),m] * 1e-4
                            v.RelocateNoAdapt[TimestepIndex(i),m] = v.RelocateNoAdapt[TimestepIndex(i),m] * 1e-4
                        end
                        v.StormCapitalNoAdapt[TimestepIndex(i),m] = v.StormCapitalNoAdapt[TimestepIndex(i),m] * 1e-4
                        v.StormPopNoAdapt[TimestepIndex(i),m] = v.StormPopNoAdapt[TimestepIndex(i),m] * 1e-4
    
                        v.NoAdaptCost[TimestepIndex(i),m] = v.WetlandNoAdapt[TimestepIndex(i),m] + v.FloodNoAdapt[TimestepIndex(i),m] +  v.RelocateNoAdapt[TimestepIndex(i),m] + v.StormCapitalNoAdapt[TimestepIndex(i),m] + v.StormPopNoAdapt[TimestepIndex(i),m]
        
    
                    end

                    if is_first(t)
                        v.NPVNoAdapt[t,m] = sum( [ v.discountfactor[TimestepIndex(j)] * v.NoAdaptCost[TimestepIndex(j),m]*10 for j in t_range] )
                    else
                        # Compute NPV Relative to planner's perspective (discounting relative to time t)
                        v.NPVNoAdapt[t,m] = sum([v.discountfactor[TimestepIndex(findind(j,t_range))] *v.NoAdaptCost[TimestepIndex(j),m]*10 for j in t_range])
                        #v.NPVNoAdapt[gettime(t)-1,m] + sum( [ v.discountfactor[j] * v.NoAdaptCost[j,m] for j in t_range] )
                    end

                    for j in t_range
                        v.NPVNoAdapt[TimestepIndex(j),m]=v.NPVNoAdapt[t,m]
                    end


    
                    # ** Calculate Protection and Retreat Costs for Each Adaptation Option **
                    lslrPlan_at = p.lslr[TimestepIndex(at_next),m]
                    lslrPlan_atprev = p.lslr[t,m]
                    
                    for i in 1:length(p.adaptoptions)
                        if is_first(t)
                            Rprev = calcHorR(-2, p.adaptoptions[i], p.lslr[TimestepIndex(1),m], p.surgeexposure[m,:], p.adaptoptions) 
                            v.R[t, m, i] = calcHorR(-2, p.adaptoptions[i], lslrPlan_at, p.surgeexposure[m,:], p.adaptoptions)
                        else
                            if p.fixed==false
                                Rprev=v.OptimalR[TimestepIndex(gettime(t)-1),m]
                                # Assumption: prior protection does not count because it is no longer maintained 
                                v.R[t,m,i] = max(v.OptimalR[TimestepIndex(gettime(t)-1),m], calcHorR(-2, p.adaptoptions[i], lslrPlan_at, p.surgeexposure[m,:], p.adaptoptions))
                            else
                                Rprev = v.R[convert(Int,p.at[at_index_prev]),m, i]
                                v.R[t, m, i] = calcHorR(-2, p.adaptoptions[i], lslrPlan_at, p.surgeexposure[m,:], p.adaptoptions)
                            end
                            
                        end
                        v.SIGMA[t,m, i+1] = (p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, v.R[t,m,i] - p.lslr[t,m]))))
                        v.coastAreaRetreat[t,m,i] = calcCoastArea(v.areaparams[m,:], v.R[t,m,i])

                        v.FloodRetreat[t,m,i] = (p.tstep/atstep) * (atstep * v.landvalue[t,m]*.04 * calcCoastArea(v.areaparams[m,:], v.R[t,m,i]) +          
                            max(0,calcCoastArea(v.areaparams[m,:], v.R[t,m,i]) - calcCoastArea(v.areaparams[m,:], Rprev))* 
                            (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[t,m]) * 1e-4
                        
                        v.RelocateRetreat[t,m,i] = (p.tstep / atstep) * 
                            max(0, calcCoastArea(v.areaparams[m,:], v.R[t,m,i]) - calcCoastArea(v.areaparams[m,:], Rprev)) * 
                            (p.movefactor * v.ypc_seg[t,m] * 1e-6 * v.popdens_seg[t,m] +
                            p.capmovefactor * p.mobcapfrac * v.capital[t,m] + p.democost * (1 - p.mobcapfrac ) * v.capital[t,m]) * 1e-4
           
                        v.DryLandLossRetreat[t,m,i] = max(0,v.coastAreaRetreat[t,m,i]) # Already takes into account prior adaptation  
                       
                        if p.adaptoptions[i] >= 10 || p.adaptoptions[i]==0
                            
                            if is_first(t)
                               
                               # Hprev = max(p.refA_H[m],calcHorR(-1, p.adaptoptions[i], p.lslr[1,m], p.surgeexposure[m,:], p.adaptoptions))
                                Hprev = calcHorR(-1, p.adaptoptions[i], p.lslr[TimestepIndex(1),m], p.surgeexposure[m,:], p.adaptoptions)
                                v.H[t,m, i-1] = calcHorR(-1, p.adaptoptions[i], lslrPlan_at, p.surgeexposure[m,:], p.adaptoptions)
                                v.SIGMA[t,m,(i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0,p.lslr[t,m])) / (1. + p.psigA[m] * exp(p.psigB[m] * max(0,(v.H[t,m, i-1] - p.lslr[t,m]))))
                                v.FloodProtect[t,m] = 0
                            else
                                if p.fixed==false
                                    Hprev = v.OptimalH[TimestepIndex(gettime(t)-1),m]
                                    ### Assumption: any prior retreat is credited toward required height, since not starting from original position on coast
                                    lslrPlan_Prot = lslrPlan_at - v.OptimalR[TimestepIndex(gettime(t)-1),m]
                                    v.H[t,m,i-1] = max(v.OptimalH[TimestepIndex(gettime(t)-1),m],calcHorR(-1,p.adaptoptions[i], lslrPlan_Prot, p.surgeexposure[m,:], p.adaptoptions))
                                    v.SIGMA[t,m,(i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0,p.lslr[t,m])) / (1. + p.psigA[m] * exp(p.psigB[m] * max(0,(v.H[t,m, i-1]+v.OptimalR[TimestepIndex(gettime(t)-1),m] - p.lslr[t,m]))))
                                    v.FloodProtect[t,m] = p.tstep * v.landvalue[t,m]*.04 * v.DryLandLossOptimal[TimestepIndex(gettime(t)-1),m]
                                else
                                    Hprev = v.H[convert(Int,p.at[at_index_prev]),m,i-1]
                                    v.H[t,m, i-1] = calcHorR(-1, p.adaptoptions[i], lslrPlan_at, p.surgeexposure[m,:], p.adaptoptions)
                                    v.SIGMA[t,m,(i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0,p.lslr[t,m])) / (1. + p.psigA[m] * exp(p.psigB[m] * max(0,(v.H[t,m, i-1] - p.lslr[t,m]))))
                                    v.FloodProtect[t,m] = 0
                                end
                                
                            end
                          
                            # Island protection costs are higher
                            if isisland(m,p.xsc)==1
                                pc = 2*p.pc0*p.cci[rgn_ind]
                            else
                                pc = p.pc0 * p.cci[rgn_ind]
                            end
    
                            v.Construct[t,m,i-1] = (p.tstep/atstep) * 
                                (p.length[m] * pc * (p.pcfixed + (1- p.pcfixed)*(v.H[t,m, i-1]^2 - Hprev^2) + 
                                p.mc*atstep*v.H[t,m, i-1]) + p.length[m] * 1.7 * v.H[t,m, i-1] * v.landvalue[t,m]*.04/2*atstep) * 1e-4
                            
                            if Hprev >= v.H[t,m,i-1]
                                v.H[t,m,i-1] = Hprev
                                # Just maintenance cost + land value
                                v.Construct[t,m,i-1] = (p.tstep/atstep) * (p.length[m] * pc *p.mc*atstep*v.H[t,m, i-1]+ p.length[m] * 1.7 * v.H[t,m, i-1] * v.landvalue[t,m]*.04/2*atstep) * 1e-4 
                            end
      
                        end
    
                        for j in t_range
                            v.R[TimestepIndex(j),m,i] = v.R[t,m,i]
                            v.SIGMA[TimestepIndex(j),m,i+1] = (p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, v.R[TimestepIndex(j),m,i] - p.lslr[TimestepIndex(j),m]))))
                            v.coastAreaRetreat[TimestepIndex(j),m,i] = v.coastAreaRetreat[t,m,i]
                            
                            if p.fixed==false && !(is_first(t))
                                v.WetlandRetreat[TimestepIndex(j),m] = p.tstep * v.wetlandservice[TimestepIndex(j),rgn_ind]* max(v.WetlandLossOptimal[TimestepIndex(gettime(t)-1),m],v.wetlandloss[TimestepIndex(i),m] * min(v.coastArea[TimestepIndex(i),m], p.wetland[m]))
                            else
                                v.WetlandRetreat[TimestepIndex(j),m] = p.tstep * v.wetlandservice[TimestepIndex(j),rgn_ind] * v.wetlandloss[TimestepIndex(j),m] * min(v.coastArea[TimestepIndex(j),m], p.wetland[m])
                            end
    
                            v.StormCapitalRetreat[TimestepIndex(j),m,i] = p.tstep * (1 - v.ρ[TimestepIndex(j),rgn_ind]) * v.SIGMA[TimestepIndex(j),m,i+1]* v.capital[TimestepIndex(j),m] 
                            v.StormPopRetreat[TimestepIndex(j),m,i] =  p.tstep * (1 - v.ρ[TimestepIndex(j),rgn_ind]) * v.SIGMA[TimestepIndex(j),m,i+1]* v.popdens_seg[TimestepIndex(j),m] * v.vsl[TimestepIndex(j),m] * p.floodmortality
                            v.StormLossRetreat[TimestepIndex(j),m,i] = p.tstep * (1 - v.ρ[TimestepIndex(j),rgn_ind]) * v.SIGMA[TimestepIndex(j),m,i+1]* v.popdens_seg[TimestepIndex(j),m] * p.floodmortality
    
                            v.FloodRetreat[TimestepIndex(j),m, i] = v.FloodRetreat[t,m,i]
                            v.RelocateRetreat[TimestepIndex(j),m,i] = v.RelocateRetreat[t,m,i]
                            v.DryLandLossRetreat[TimestepIndex(j),m,i] = v.DryLandLossRetreat[t,m,i]
    
                            # Put all other costs intp $Billions from $M and divide by 10
                            v.StormCapitalRetreat[TimestepIndex(j),m,i]  = v.StormCapitalRetreat[TimestepIndex(j),m,i]  * 1e-4
                            v.StormPopRetreat[TimestepIndex(j),m,i]  = v.StormPopRetreat[TimestepIndex(j),m,i]  * 1e-4
                            v.WetlandRetreat[TimestepIndex(j),m] = v.WetlandRetreat[TimestepIndex(j),m] * 1e-4
                                        
                            v.RetreatCost[TimestepIndex(j),m, i] = v.FloodRetreat[TimestepIndex(j),m,i] + v.RelocateRetreat[TimestepIndex(j),m,i] + v.StormCapitalRetreat[TimestepIndex(j),m,i] + v.StormPopRetreat[TimestepIndex(j),m,i] + v.WetlandRetreat[TimestepIndex(j),m]
                                
                            if p.adaptoptions[i] >= 10 || p.adaptoptions[i]==0
                                v.H[TimestepIndex(j),m, i-1] = v.H[t,m, i-1]
                                
                                if p.fixed==false && !(is_first(t))
                                    v.SIGMA[TimestepIndex(j),m,(i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0,p.lslr[TimestepIndex(j),m])) / 
                                    (1. + p.psigA[m] * exp(p.psigB[m] * max(0,(v.H[TimestepIndex(j),m, i-1]+v.OptimalR[TimestepIndex(gettime(t)-1),m] - p.lslr[TimestepIndex(j),m]))))
                                    v.FloodProtect[TimestepIndex(j),m]=p.tstep * v.landvalue[TimestepIndex(j),m]*.04 * v.DryLandLossOptimal[TimestepIndex(gettime(t)-1),m]
                                   
                                else
                                    v.SIGMA[TimestepIndex(j),m,(i-1)+7] = (p.psig0[m] + p.psig0coef[m] * max(0,p.lslr[TimestepIndex(j),m])) / 
                                    (1. + p.psigA[m] * exp(p.psigB[m] * max(0,(v.H[TimestepIndex(j),m, i-1] - p.lslr[TimestepIndex(j),m]))))
                                    v.FloodProtect[TimestepIndex(j),m]=v.FloodProtect[t,m]
                                end
                               
    
                                v.WetlandProtect[TimestepIndex(j),m] = p.tstep * p.wetland[m] .* v.wetlandservice[TimestepIndex(j),rgn_ind]
                                        
                                v.StormCapitalProtect[TimestepIndex(j),m,i-1] = p.tstep * (1 - v.ρ[TimestepIndex(j),rgn_ind]) * v.SIGMA[TimestepIndex(j),m,(i-1)+7] * v.capital[TimestepIndex(j),m]                                              
                                v.StormPopProtect[TimestepIndex(j),m,i-1] =  p.tstep * (1 - v.ρ[TimestepIndex(j),rgn_ind]) * v.SIGMA[TimestepIndex(j),m,(i-1)+7] * v.popdens_seg[TimestepIndex(j),m] * v.vsl[TimestepIndex(j),m] * p.floodmortality
                                v.StormLossProtect[TimestepIndex(j),m,i-1] = p.tstep * (1 - v.ρ[TimestepIndex(j),rgn_ind]) * v.SIGMA[TimestepIndex(j),m,(i-1)+7] * v.popdens_seg[TimestepIndex(j),m] * p.floodmortality
                                    
                                v.Construct[TimestepIndex(j),m,i-1] = v.Construct[t,m, i-1]
    
                                # Put all other costs intp $Billions from $M and divide by 10
                                # Note this is an annual protect cost ($B/year)
                                v.WetlandProtect[TimestepIndex(j),m] = v.WetlandProtect[TimestepIndex(j),m] * 1e-4
                                v.StormCapitalProtect[TimestepIndex(j),m,i-1] = v.StormCapitalProtect[TimestepIndex(j),m,i-1] * 1e-4
                                v.StormPopProtect[TimestepIndex(j),m,i-1] = v.StormPopProtect[TimestepIndex(j),m,i-1] * 1e-4
                                v.FloodProtect[TimestepIndex(j),m] = v.FloodProtect[TimestepIndex(j),m] * 1e-4
                                                    
                                v.ProtectCost[TimestepIndex(j),m,i-1] = v.Construct[TimestepIndex(j),m,i-1] + v.WetlandProtect[TimestepIndex(j),m] + v.StormCapitalProtect[TimestepIndex(j),m,i-1] + v.StormPopProtect[TimestepIndex(j),m,i-1] + v.FloodProtect[TimestepIndex(j),m]
    
                            end
    
                        end
    
                        if is_first(t)
                            v.NPVRetreat[t,m,i] = sum([v.discountfactor[TimestepIndex(j)] * v.RetreatCost[TimestepIndex(j),m,i] * 10 for j in t_range])
                        else
                            # Compute NPV Relative to planner's perspective (discounting relative to time t)
                            v.NPVRetreat[t,m,i] = sum([v.discountfactor[TimestepIndex(findind(j,t_range))]*v.RetreatCost[TimestepIndex(j),m,i]*10 for j in t_range])
                            #v.NPVRetreat[gettime(t)-1,m,i] + sum([v.discountfactor[j] * v.RetreatCost[j,m,i] for j in t_range])
                        end

                        
                        for j in t_range
                            v.NPVRetreat[TimestepIndex(j),m,i] = v.NPVRetreat[t,m,i]
                        end
    
                        if p.adaptoptions[i] >=10 || p.adaptoptions[i]==0
                            if is_first(t)
                                v.NPVProtect[t,m,i-1] = sum( [ v.discountfactor[TimestepIndex(j)] * v.ProtectCost[TimestepIndex(j),m,i-1] *10 for j in t_range] ) # Protect
                            else
                                # Compute NPV Relative to planner's perspective (discounting relative to time t)
                                v.NPVProtect[t,m,i-1] = sum([v.discountfactor[TimestepIndex(findind(j,t_range))]*v.ProtectCost[TimestepIndex(j),m,i-1]*10 for j in t_range])
                                #v.NPVProtect[gettime(t)-1,m,i-1] + sum( [ v.discountfactor[j] * v.ProtectCost[j,m,i-1] for j in t_range] ) # Protect
                            end

                            
                            for j in t_range
                                v.NPVProtect[TimestepIndex(j),m,i-1] = v.NPVProtect[t,m,i-1]
                            end
                        end
                    end

                    # ** Choose Least Cost Option **
                    if gettime(t)>1 && p.fixed
                        # if p.fixed==T and t>1, take first-period choices 
                        for j in t_range
                            v.OptimalProtectLevel[TimestepIndex(j), m] = v.OptimalProtectLevel[TimestepIndex(1), m]
                            v.OptimalRetreatLevel[TimestepIndex(j),m] = v.OptimalRetreatLevel[TimestepIndex(1),m] 
                            v.OptimalOption[TimestepIndex(j),m] = v.OptimalOption[TimestepIndex(1),m]
                            v.OptimalLevel[TimestepIndex(j),m] = v.OptimalLevel[TimestepIndex(1),m]
                        end
                    else
                        # If p.fixed==F or if p.fixed==T and t==1, calculate optimal level.
                        if p.allowMaintain==true
                            
                            protectInd = findmin(v.NPVProtect[TimestepIndex(Int(p.at[at_index])),m,:])[2]
                            retreatInd = findmin(v.NPVRetreat[TimestepIndex(Int(p.at[at_index])),m,:])[2]
                        else
                            protDims=size(v.NPVProtect)[3]
                            retDims=size(v.NPVRetreat)[3]
                            protectInd = findmin(v.NPVProtect[TimestepIndex(Int(p.at[at_index])),m,1:protDims-1])[2]
                            retreatInd = findmin(v.NPVRetreat[TimestepIndex(Int(p.at[at_index])),m,1:retDims-1])[2]
                        end
                        for j in t_range
                            v.OptimalProtectLevel[TimestepIndex(j), m] = p.adaptoptions[protectInd+1]
                        end

                        if p.noRetreat==true 
                            minLevels = [p.adaptoptions[protectInd+1], 0]
                            choices = [v.NPVProtect[TimestepIndex(Int(p.at[at_index])),m,protectInd], v.NPVNoAdapt[TimestepIndex(Int(p.at[at_index])),m]]

                            leastcost = -1 * findmin(choices)[2]
                            if leastcost==-2
                                leastcost=-3 # Account for retreat being removed from choice set 
                            end
                            leastlevel = minLevels[findmin(choices)[2]]
                        else
                            for j in t_range
                                v.OptimalRetreatLevel[TimestepIndex(j),m] = p.adaptoptions[retreatInd]
                            end
                            minLevels = [p.adaptoptions[protectInd+1], p.adaptoptions[retreatInd], 0]
                            
                            choices = [v.NPVProtect[TimestepIndex(Int(p.at[at_index])),m,protectInd], v.NPVRetreat[TimestepIndex(Int(p.at[at_index])),m,retreatInd], v.NPVNoAdapt[TimestepIndex(Int(p.at[at_index])),m]]
                            leastcost = -1 * findmin(choices)[2]
                            leastlevel = minLevels[findmin(choices)[2]]
                        end
                        for j in t_range
                            v.OptimalOption[TimestepIndex(j),m] = leastcost
                            v.OptimalLevel[TimestepIndex(j),m] = leastlevel
                        end
                    end
                    
                    # Assign costs to optimal variables
                    if v.OptimalOption[t,m]==-1
              
                        # Protect Cost
                        protInd = findall(i->i==v.OptimalLevel[t,m], p.adaptoptions)[1]-1
                        for j in t_range
                            v.OptimalCost[TimestepIndex(j),m] = v.ProtectCost[TimestepIndex(j),m,protInd] 
                            # Assign Subcosts 
                            v.OptimalStormCapital[TimestepIndex(j),m] = v.StormCapitalProtect[TimestepIndex(j),m,protInd]
                            v.OptimalStormPop[TimestepIndex(j),m] = v.StormPopProtect[TimestepIndex(j),m,protInd]
                            v.OptimalConstruct[TimestepIndex(j),m] = v.Construct[TimestepIndex(j),m,protInd]
                            v.OptimalWetland[TimestepIndex(j),m] = v.WetlandProtect[TimestepIndex(j),m]
                            v.OptimalRelocate[TimestepIndex(j),m] = 0
                            v.OptimalFlood[TimestepIndex(j),m] = v.FloodProtect[TimestepIndex(j),m] 

                            # Assign Alternative Metrics 
                            # Assume once seawall is built, wetland area is permanently destroyed 
                            v.WetlandLossOptimal[TimestepIndex(j),m] = p.wetland[m]
                        end                        

                        if gettime(t)==1
                            for i in t_range
                                v.DryLandLossOptimal[TimestepIndex(i),m] = 0
                                v.OptimalH[TimestepIndex(i),m]=max(0,v.H[TimestepIndex(i),m,protInd])
                                v.OptimalR[TimestepIndex(i),m]=0
                                if i==1
                                    v.StormLossOptimal[TimestepIndex(i),m] = v.StormLossProtect[TimestepIndex(i),m,protInd]
                                    
                                else 
                                    v.StormLossOptimal[TimestepIndex(i),m] = v.StormLossOptimal[TimestepIndex(i-1),m] +  v.StormLossProtect[TimestepIndex(i),m,protInd]
                                    
                                end
                            end
                        else
                            for j in t_range
                                v.OptimalR[TimestepIndex(j),m]=max(0,v.OptimalR[TimestepIndex(gettime(t)-1),m])
                            end
                            for i in t_range
                                v.OptimalH[TimestepIndex(i),m] = max(v.H[TimestepIndex(i),m,protInd],v.OptimalH[TimestepIndex(gettime(t)-1),m])
                                
                                v.DryLandLossOptimal[TimestepIndex(i),m] = max(0, v.DryLandLossOptimal[TimestepIndex(i-1),m])
                                v.StormLossOptimal[TimestepIndex(i),m] = v.StormLossOptimal[TimestepIndex(i-1),m] + v.StormLossProtect[TimestepIndex(i),m,protInd]
                            end
                        end

                    elseif v.OptimalOption[t,m]==-2
                        # Retreat Cost
                      
                        retInd=findall(i->i==v.OptimalLevel[t,m], p.adaptoptions)[1]
                        for j in t_range
                            v.OptimalCost[TimestepIndex(j),m] = v.RetreatCost[TimestepIndex(j),m, retInd]
                            # Assign Subcosts  
                            v.OptimalStormCapital[TimestepIndex(j),m] = v.StormCapitalRetreat[TimestepIndex(j),m, retInd]
                            v.OptimalStormPop[TimestepIndex(j),m] = v.StormPopRetreat[TimestepIndex(j),m, retInd]
                            v.OptimalConstruct[TimestepIndex(j),m] = 0
                            v.OptimalWetland[TimestepIndex(j),m] = v.WetlandRetreat[TimestepIndex(j),m]
                            v.OptimalFlood[TimestepIndex(j),m] = v.FloodRetreat[TimestepIndex(j),m,retInd]
                            v.OptimalRelocate[TimestepIndex(j),m] = v.RelocateRetreat[TimestepIndex(j),m,retInd]
                        end
                      
                        if is_first(t)
                            
                            for j in t_range
                                v.DryLandLossOptimal[TimestepIndex(j),m] = v.DryLandLossRetreat[TimestepIndex(j),m,retInd]
                                v.OptimalH[TimestepIndex(j),m]=0
                            end
                           
                            for i in t_range
                                v.OptimalR[TimestepIndex(i),m]=max(0,v.R[TimestepIndex(i),m,retInd])
                                v.WetlandLossOptimal[TimestepIndex(i),m] = v.wetlandloss[TimestepIndex(i),m] * min(v.coastArea[TimestepIndex(i),m], p.wetland[m])
                                if i==1
                                    v.StormLossOptimal[TimestepIndex(i),m] = v.StormLossRetreat[TimestepIndex(i),m,findall(k->k==v.OptimalLevel[t,m], p.adaptoptions)[1]] 
                                else 
                                    v.StormLossOptimal[TimestepIndex(i),m] = v.StormLossOptimal[TimestepIndex(i-1),m] + v.StormLossRetreat[TimestepIndex(i),m,retInd] 
                                end
                            end
                        else
                            for j in t_range
                                v.OptimalH[TimestepIndex(j),m]=v.OptimalH[TimestepIndex(gettime(t)-1),m]
                            end
                            
                            for i in t_range
                                v.OptimalR[TimestepIndex(i),m]=max(v.R[TimestepIndex(i),m,retInd],v.OptimalR[TimestepIndex(gettime(t)-1),m])
                                # Cumulative total wetland area lost; if protected previously, all wetland is lost 
                                v.WetlandLossOptimal[TimestepIndex(i),m] = max(v.WetlandLossOptimal[TimestepIndex(gettime(t)-1),m], v.wetlandloss[TimestepIndex(i),m]*min(v.coastArea[TimestepIndex(i),m], p.wetland[m]))
                                v.StormLossOptimal[TimestepIndex(i),m] = v.StormLossOptimal[TimestepIndex(i-1),m] + p.tstep * (1 - v.ρ[TimestepIndex(i),rgn_ind ]) * v.popdens_seg[TimestepIndex(i),m] * p.floodmortality * v.SIGMA[TimestepIndex(i),m,retInd] 
                                v.DryLandLossOptimal[TimestepIndex(i),m] = max(v.DryLandLossOptimal[TimestepIndex(i-1),m],v.DryLandLossRetreat[TimestepIndex(i),m,retInd])
                            end
                           
                        end

                    else
                        # No Adaptation 
                    
                        for j in t_range
                            v.OptimalCost[TimestepIndex(j),m] = v.NoAdaptCost[TimestepIndex(j),m]
                            # Assign Subcosts 
                            v.OptimalStormCapital[TimestepIndex(j),m] = v.StormCapitalNoAdapt[TimestepIndex(j),m]
                            v.OptimalStormPop[TimestepIndex(j),m] = v.StormPopNoAdapt[TimestepIndex(j),m]
                            v.OptimalConstruct[TimestepIndex(j),m] = 0
                            v.OptimalWetland[TimestepIndex(j),m] = v.WetlandNoAdapt[TimestepIndex(j),m]
                            v.OptimalFlood[TimestepIndex(j),m] = v.FloodNoAdapt[TimestepIndex(j),m]
                            v.OptimalRelocate[TimestepIndex(j),m] = v.RelocateNoAdapt[TimestepIndex(j),m]
                        end

                        if is_first(t)
                            for j in t_range
                                v.OptimalH[TimestepIndex(j),m]=0
                                v.DryLandLossOptimal[TimestepIndex(j),m] = v.DryLandLossNoAdapt[TimestepIndex(j),m]
                            end
                            v.OptimalR[t,m]=max(0,p.lslr[t,m])

                            for i in t_range
                                if i>1
                                    v.OptimalR[TimestepIndex(i),m]=max(v.OptimalR[t,m],v.OptimalR[TimestepIndex(i-1),m],p.lslr[TimestepIndex(i),m])
                                else
                                    v.OptimalR[TimestepIndex(i),m]=max(v.OptimalR[t,m],p.lslr[TimestepIndex(i),m])
                                end
                                v.WetlandLossOptimal[TimestepIndex(i),m] = v.wetlandloss[TimestepIndex(i),m] * min(v.coastArea[TimestepIndex(i),m], p.wetland[m])
                                if i==1
                                    v.StormLossOptimal[TimestepIndex(i),m] = v.StormLossNoAdapt[TimestepIndex(i),m]
                                else 
                                    v.StormLossOptimal[TimestepIndex(i),m] = v.StormLossOptimal[TimestepIndex(i-1),m] + v.StormLossNoAdapt[TimestepIndex(i),m]
                                end
                            end
                        else
                            for j in t_range
                                v.OptimalH[TimestepIndex(j),m]=v.OptimalH[TimestepIndex(gettime(t)-1),m]
                            end
                            
                            for i in t_range
                                
                                v.OptimalR[TimestepIndex(i),m]=max(v.OptimalR[TimestepIndex(gettime(t)-1),m],v.OptimalR[TimestepIndex(i-1),m],p.lslr[TimestepIndex(i),m])
                                
                                v.WetlandLossOptimal[TimestepIndex(i),m] = max(v.WetlandLossOptimal[TimestepIndex(gettime(t)-1),m], v.wetlandloss[TimestepIndex(i),m]*min(v.coastArea[TimestepIndex(i),m], p.wetland[m]))
                                v.StormLossOptimal[TimestepIndex(i),m] = v.StormLossOptimal[TimestepIndex(i-1),m] + v.StormLossNoAdapt[TimestepIndex(i),m]
                                v.DryLandLossOptimal[TimestepIndex(i),m] = max(v.DryLandLossOptimal[TimestepIndex(i-1),m], v.DryLandLossNoAdapt[TimestepIndex(i),m])
                            end
                           
                        end
                    end

                    if last==1
                        v.NPVOptimal[m] = sum( [ v.discountfactor[TimestepIndex(j)] * v.OptimalCost[TimestepIndex(j),m] *10 for j in 1:p.ntsteps] )
                        
                    end
                end
            end  # end segment loop 
            
            if last==1
                v.NPVOptimalTotal = sum(v.NPVOptimal)
            end
        end # end if t in adaptation period statement  
    end # end function 
end # end module 


# Helper functions
function growthrate(x1, x2)
    epsilon = 1.e-9
    return (x2 / (x1 + epsilon) - 1.)
end

function calcCoastArea(areaparams, var)

    area = (areaparams[1]*max(0,min(0.5,var-0))
    +(areaparams[1]+areaparams[2])/2*max(0,min(1,var-0.5))
    +areaparams[2]*max(0,min(0.5,var-1.5))
    +areaparams[3]*max(0,min(1,var-2))
    +areaparams[4]*max(0,min(1,var-3))
    +areaparams[5]*max(0,min(1,var-4))
    +areaparams[6]*max(0,min(1,var-5))
    +areaparams[7]*max(0,min(1,var-6))
    +areaparams[8]*max(0,min(1,var-7))
    +areaparams[9]*max(0,min(1,var-8))
    +areaparams[10]*max(0,min(1,var-9))
    +areaparams[11]*max(0,min(1,var-10))
    +areaparams[12]*max(0,min(1,var-11))
    +areaparams[13]*max(0,min(1,var-12))
    +areaparams[14]*max(0,min(1,var-13))
    +areaparams[15]*max(0,var-14))

    return area

end

function localrate(lslr1, lslr2, tstep)
    return max(0, (lslr2 - lslr1))/tstep
end

function getsegments(rgn_name, xsc)

    segs = collect(keys(filter( (k,v) -> v[1]==rgn_name,xsc)))
    return segs

end

function getregion(seg_ind, xsc)
    rgn = xsc[seg_ind][1]
    
    return rgn
end

function isgreenland(seg_ind, xsc)
    greenland = xsc[seg_ind][2]
    return greenland
end

function isisland(seg_ind, xsc)
    island = xsc[seg_ind][3]
    return island
end

function calcHorR(option, level, lslrPlan, surgeExpLevels, adaptOptions)
    
    ind = findind(level, adaptOptions)

    if option==-1 && level ==10
        # Protect height differs from retreat radius only in case of 10 yr surge exposure
        H = max(0, lslrPlan + surgeExpLevels[ind] / 2)
        return H
    elseif level==0 # Maintain existing defenses 
        H_R=0
        return H_R
    else
        H_R = max(0, lslrPlan + surgeExpLevels[ind])
        return H_R
    end
end

function pos(x)
    return max(0,x)
end