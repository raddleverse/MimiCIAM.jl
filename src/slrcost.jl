# Catherine Ledna
# Jan 24, 2019
##------------------------------------------------------------------------
# CIAM Model
#------------------------------------------------------------------------
# Implements CIAM model adapted from Diaz, 2016. 
#------------------------------------------------------------------------
# TODO update to work with Julia 1.1 and Mimi 

using Mimi

@defcomp slrcost begin
    # --- Indices ---
    regions = Index()
    segments = Index()                      
    adaptPers = Index()

   
    # --- Region / segment mapping ---
    xsc::Dict{Any, Any} = Parameter()        # Region to segment mapping    

    # ---Time-related Parameters---
    tstep = Parameter()                     # Length of individual time-step (years)
    at = Parameter( index = [adaptPers])    # Array of time indices that mark starts of adaptation periods 
    ntsteps::Int = Parameter()              # Number of time-steps     

    # ---Model Parameters ---
    fixed::Bool = Parameter()               # Run model as fixed (T) or flexible (F) with respect to adaptation

    # ---Socioeconomic Parameters---
    pop = Parameter(index = [time, regions])           # Population of region (million people) (from MERGE)
    refpopdens = Parameter( index = [regions])         # Reference population density of region (people / km^2)
    refpopdens_usa = Parameter()                       # Reference population density of USA (people/km^2) 
    popdens = Parameter( index = [segments])           # Pop density of segment in time t = 1 (people/km^2)
    ypcc = Parameter(index = [time, regions])          # GDP per capita per region ($2010 per capita)
    ypc_usa = Parameter(index = [time])                # GDP per capita in USA; used as benchmark ($2010 per capita)
     
    popdens_seg = Variable(index = [time, segments])          # Population density of segment extrapolated forward in time (people / km^2)    
    ypc_seg = Variable(index = [time, segments])              # GDP per capita by segment ($2010 per capita) (multiplied by scaling factor)

    # ---Land Parameters---  
    landinput::Bool = Parameter()                   # Set to T for FUND or F for GTAP
         
    gtapland = Parameter( index = [regions])        # GTAP land value in 2007 (million 2010$ / km^2)
    gtapland_canada = Parameter()                   # GTAP land value in 2007 for Canada (million 2010$ / km^2)
    fundland_canada = Parameter()                   # FUND land value in 1995 for Canada (million 2010$ / km^2)
    dvbm = Parameter()                              # FUND value of OECD dryland per Darwin et al 1995 converted from $1995 ($2010M per sqkm) (5.376)
    kgdp = Parameter()                              # Capital output ratio (per MERGE) (3 by default) 
    discountrate = Parameter()                      # Discount rate (0.04 by default)
    depr = Parameter()                              # Fraction of capital that has not been depreciated over adaptation period (retreat cases)
    
    
    landdata = Variable( index = [regions])         # Takes on value of either fundland or gtapland
    landdata_canada = Variable()
    fundland = Variable( index = [regions])         # FUND land value in 1995 (calculated in run_timestep) (million 2010$ / km^2), 
                                                    #   Q maybe import directly? 

    land_appr_canada = Parameter(index = [time])    # Canada land appreciation rate (used for Greenland)
    land_appr = Variable(index = [time, regions])   # Land appreciation rate (calculated as regression by Yohe ref Abraham and Hendershott) 
    coastland = Variable(index = [time, regions])  # Coastal land value (function of interior land value * scaling factor) ($2010M per sqkm)
    landvalue = Variable(index = [time, regions])  # Total endowment value of land ($2010M per sqkm)


    ρ = Variable(index = [time, regions])           # Country-wide resilience parameter (logistic function related to GDP)
    capital = Variable(index = [time, regions])    # Total endowment value of capital stock (million $2010 / km^2)
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
    vsl = Variable(index = [time, segments])     # Value of statistica life (million 2010$)
                                           
    # ---Wetland Loss Parameters---
    wbvm = Parameter()                                      # Annual value of wetland services (million 2010$ / km^2 / yr); (Brander et al 2006)  (0.376) 
    wetland = Parameter( index = [segments])            # Initial wetland area in coastal segment (km^2)
    wmaxrate = Parameter()                                  # Maximum rate of wetland accretion (m per yr) per Kirwan et al 2010 (0.01)
    
    wetlandservice = Variable(index = [time, regions])      # Annual value of wetland services TODO change doc
    wetlandloss = Variable(index = [time, segments])        # Fractional loss of wetland due to slr


    # ---Sea Level Rise Parameters---
    lslr = Parameter(index = [time, segments])                # Local sea level rise (m) 

    adaptoptions = Parameter(index = [5])                     # Index of available adaptation levels for protect and retreat (0 is no adaptation)
    s10 = Parameter(index = [segments])
    s100 = Parameter(index = [segments])
    s1000 = Parameter(index = [segments])
    smax = Parameter(index = [segments])
    surgeExposure::Float64 = Variable( index = [segments, 4])     # Storm surge exposure levels (corresponding to each designated adaptation option)
    

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

    coastArea = Variable(index=[time, segments])            # Calculate remaining coastal area after slr (m^2)
    
    # ---Intermediate Variables---
    WetlandNoAdapt = Variable(index = [time, segments])
    FloodNoAdapt = Variable(index = [time, segments])
    StormNoAdapt = Variable(index = [time, segments])
    RelocateNoAdapt = Variable(index = [time, segments])

    Construct = Variable(index = [time, segments, 4])
    WetlandProtect = Variable(index = [time, segments])
    StormProtect = Variable(index = [time, segments,4])
    H = Variable(index = [time, segments, 4])
    
    WetlandRetreat = Variable(index = [time, segments])
    StormRetreat = Variable(index = [time, segments,5])
    FloodRetreat = Variable(index = [time, segments, 5])
    RelocateRetreat = Variable(index = [time, segments, 5])
    R = Variable(index = [time, segments, 5])

    # ---Decision Variables---   
    NoAdaptCost = Variable(index = [time, segments])         # Cost of not adapting (e.g. reactive retreat) (2010$)
    ProtectCost = Variable(index = [time, segments, 4])      # Total cost of protection at each level      
    RetreatCost = Variable(index = [time, segments, 5])      # Total cost of retreat at each level   
    OptimalFixedCost = Variable(index = [time, segments])        # Fixed optimal cost based on NPV in period 1   
    OptimalFixedLevel = Variable(index = [segments])        # Fixed optimal level (1,10,100,1000,10000)
    OptimalFixedOption = Variable(index = [segments])       # Fixed adaptation decision (-1 - protect, -2 - retreat, -3 - no adapt) 
    NPVRetreat = Variable(index = [segments, adaptPers, 5])
    NPVProtect = Variable(index = [segments, adaptPers, 4])
    NPVNoAdapt = Variable(index = [segments, adaptPers])                     
  
    # ---Outcome Variables---
    AdaptationDecision = Variable(index = [time, segments])   # Option chosen for adaptation period
    AdaptationCost = Variable(index = [time, segments])       # Cost of option chosen for adaptation period (B 2010$ / yr) 
    AdaptationLevel = Variable(index = [time, segments])      # Level of protect or retreat (if chosen)
    RegionalCost = Variable(index = [time, regions])          # Cost of adaptation at level of region

    function run_timestep(p, v, d, t)    
        # In first period, initialize all non-adaptation dependent intermediate variables for all timesteps
        if is_first(t)
          #  1. Initialize non-region dependent intermediate variables
    
            for i in collect(1:Int(p.ntsteps)) 
               v.discountfactor[i] = 1/(1 + p.discountrate)^(p.tstep * (i-1))
            end
    
            # 2. Initialize region-dependent intermediate variables
            for r in d.regions
                # Determine land input value (true = FUND, false = GTAP)
                # TODO switch to importing fund land values as param
                if p.landinput 
                    v.fundland[r] = min(p.dvbm, max(0.005, p.dvbm * p.ypcc[t,r] * p.refpopdens[r] / (p.ypc_usa[1] * p.refpopdens_usa)))
                    v.landdata[r] = v.fundland[r]
                    v.landdata_canada = p.fundland_canada
                else
                    v.landdata[r] = p.gtapland[r]
                    v.landdata_canada = p.gtapland_canada
                end
    
                v.wetlandservice[t,r] = p.wbvm * ((p.ypcc[t,r] / p.ypc_usa[1])^1.16 * (p.refpopdens[r] /27.59)^0.47) 
                v.ρ[t,r] = p.ypcc[t,r] / (p.ypcc[t,r] + p.ypc_usa[1])
                v.land_appr[t,r] = 1.
    
                for i in collect(2:Int(p.ntsteps))
                    v.land_appr[i,r] = v.land_appr[i-1,r] * exp(0.565 * growthrate(p.ypcc[i-1,r], p.ypcc[i,r]) + 0.313 * growthrate(p.pop[i-1,r], p.pop[i,r]))
                    v.wetlandservice[i,r] = v.land_appr[i,r] * v.wetlandservice[1,r]
                    v.ρ[i,r] = p.ypcc[i,r] / (p.ypcc[i,r] + p.ypc_usa[1]) 
                end    
            end
    
            # 3. Initialize segment-dependent variables 
            for m in d.segments
                rgn_ind = getregion(m, p.xsc)
     
                v.popdens_seg[t,m] = p.popdens[m]
                v.areaparams[m,:] = [p.area1[m] p.area2[m] p.area3[m] p.area4[m] p.area5[m] p.area6[m] p.area7[m] p.area8[m] p.area9[m] p.area10[m] p.area11[m] p.area12[m] p.area13[m] p.area14[m] p.area15[m]]
                v.surgeExposure[m,:] = [p.s10[m] p.s100[m] p.s1000[m] p.smax[m]] 

                # Greenland segments treated differently 
                if isgreenland(m,p.xsc)==1
                    v.ypc_seg[t,m] =22642*1.01^t   # FLAG: assumes t is an index (1-20)
                    v.vsl[t,m] = 1e-6 * 216 * p.ypc_usa[1] * (v.ypc_seg[t,m]/p.ypc_usa[1])^0.5
                    v.coastland[t,m] = (p.land_appr_canada[1] * v.landdata_canada) * max(0.5, log(1+v.popdens_seg[t,m])/log(25))
                    v.landvalue[t,m] = min(v.coastland[t,m], (p.land_appr_canada[1] * v.landdata_canada))
                else
                    v.ypc_seg[t,m] = p.ypcc[t,rgn_ind] * max(0.9, (p.popdens[m]/250.)^0.05)
                    v.vsl[t,m] = 1e-6 * 216 * p.ypc_usa[1] * (p.ypcc[t,rgn_ind]/p.ypc_usa[1])^0.5
                    v.coastland[t,m] = max(0.5, log(1+v.popdens_seg[t,m])/log(25)) * (v.land_appr[t,rgn_ind] * v.landdata[rgn_ind])  # Interior * scaling factor
                    v.landvalue[t,m] = min(v.coastland[t,m], (v.land_appr[t,rgn_ind] * v.landdata[rgn_ind]))
                end
    
                v.capital[t,m] = p.kgdp * v.ypc_seg[t,m] * v.popdens_seg[t,m] * 1e-6
                v.coastArea[t,m] = calcCoastArea(v.areaparams[m,:], p.lslr[t,m])  
                
                
                for i in 2:Int(p.ntsteps)      
                    v.popdens_seg[i,m] = v.popdens_seg[i-1,m] * (1 + growthrate(p.pop[i-1,rgn_ind], p.pop[i,rgn_ind])) 
     
                    # Special treatment for Greenland segments 
                    if isgreenland(m,p.xsc)==1
                        v.ypc_seg[i,m] =22642*1.01^i   # FLAG: assumes i is an index (1-20)
                        v.vsl[i,m] = 1e-6 * 216 * p.ypc_usa[i] * (v.ypc_seg[i,m]/p.ypc_usa[i])^0.5  
                        v.coastland[i,m] = (p.land_appr_canada[i] * v.landdata_canada) * max(0.5, log(1+v.popdens_seg[i,m])/log(25))
                        v.landvalue[i,m] = min(v.coastland[i,m], (p.land_appr_canada[i] * v.landdata_canada))
    
                    else
                        v.ypc_seg[i,m] = p.ypcc[i,rgn_ind] * max(0.9, (p.popdens[m]/250.)^0.05) # ypcc * popdens scaling factor
                        v.coastland[i,m] = max(0.5, log(1+v.popdens_seg[i,m])/log(25)) * (v.land_appr[i,rgn_ind] * v.landdata[rgn_ind])
                        v.vsl[i,m] = 1e-6 * 216 * p.ypc_usa[i] * (p.ypcc[i,rgn_ind]/p.ypc_usa[i])^0.5     
                        v.landvalue[i,m] = min(v.coastland[i,m], (v.land_appr[i,rgn_ind] * v.landdata[rgn_ind]))
     
                    end
    
                    v.capital[i,m] = p.kgdp * v.ypc_seg[m, i] * v.popdens_seg[i,m] * 1e-6 
                    v.coastArea[i,m] = calcCoastArea(v.areaparams[m,:], p.lslr[i,m])
                    v.wetlandloss[i-1,m] = min(1, (localrate(p.lslr[i-1,m], p.lslr[i,m], p.tstep)/p.wmaxrate)^2)
                   
    
                end
    
                v.wetlandloss[p.ntsteps,m] = min(1, (localrate(p.lslr[p.ntsteps-1,m], p.lslr[p.ntsteps,m], p.tstep)/p.wmaxrate)^2)  
    
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
            println(typeof(t))
            for m in d.segments
                if atstep==0
                    v.AdaptationCost[t,m] = v.AdaptationCost[t-1,m]
                    v.AdaptationDecision[t,m] = v.AdaptationDecision[t-1,m]
                    v.AdaptationLevel[t,m] = v.AdaptationLevel[t-1,m]
                else
                    rgn_ind = getregion(m, p.xsc)
     
                    # ** Calculate No Adaptation Costs **
                    for i in t_range
                        R_NoAdapt = max(0, p.lslr[m,i])
    
                        v.StormNoAdapt[i,m] = p.tstep * (1 - v.ρ[i,rgn_ind ]) * (p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, R_NoAdapt - p.lslr[i,m])))) * 
                            (v.capital[i,m] + v.popdens_seg[i,m] * v.vsl[i,m] * p.floodmortality)
                                
                        v.WetlandNoAdapt[i,m] = p.tstep * v.wetlandservice[i,rgn_ind] * v.wetlandloss[i,m] * min(v.coastArea[i,m], p.wetland[m])
    
                        if i==p.ntsteps
                            v.FloodNoAdapt[i,m] = p.tstep * v.landvalue[i-1,m]*.04 * max(0, v.coastArea[i,m]) + (max(0, v.coastArea[i,m]) - max(0, v.coastArea[i-1,m])) * 
                                (1 - p.mobcapfrac) * v.capital[i-1,m]
            
                            v.RelocateNoAdapt[i,m] = (max(0, v.coastArea[i,m]) - max(0,v.coastArea[i-1,m])) * (5 * p.movefactor * v.ypc_seg[i-1,m]*1e-6*v.popdens_seg[i-1,m] +
                                p.capmovefactor * p.mobcapfrac * v.capital[i-1,m] + p.democost * (1 - p.mobcapfrac) * v.capital[i-1,m])
                        else
                            v.FloodNoAdapt[i,m]  = p.tstep * v.landvalue[i,m]*.04 * max(0, v.coastArea[i+1,m]) + (max(0, v.coastArea[i+1,m]) - max(0, v.coastArea[i,m])) * 
                                (1 - p.mobcapfrac) * v.capital[i,m]                    
                            
                            v.RelocateNoAdapt[i,m] = (max(0, v.coastArea[i+1,m]) - max(0,v.coastArea[i,m])) * (5 * p.movefactor * v.ypc_seg[i,m]*1e-6*v.popdens_seg[i,m] +
                                p.capmovefactor * p.mobcapfrac * v.capital[i,m] + p.democost * (1 - p.mobcapfrac) * v.capital[i,m])
                        end
                                
                        # Put all costs into $Billions and divide by 10
                        v.WetlandNoAdapt[i,m] = v.WetlandNoAdapt[i,m] * 1e-4
                        v.FloodNoAdapt[i,m] = v.FloodNoAdapt[i,m] * 1e-4
                        v.RelocateNoAdapt[i,m] = v.RelocateNoAdapt[i,m] * 1e-4
                        v.StormNoAdapt[i,m] = v.StormNoAdapt[i,m] * 1e-4
    
                        v.NoAdaptCost[i,m] = v.WetlandNoAdapt[i,m] + v.FloodNoAdapt[i,m] +  v.RelocateNoAdapt[i,m] + v.StormNoAdapt[i,m]
        
    
                    end
                    v.NPVNoAdapt[at_index,m] = sum( [ v.discountfactor[j] * v.NoAdaptCost[j,m] for j in t_range] )
    
                    # ** Calculate Protectio and Retreat Costs for Each Adaptation Option **
                    lslrPlan_at = p.lslr[at_next,m]
                    println(lslrPlan_at)
                    lslrPlan_atprev = p.lslr[at_prev,m]
                    println(lslrPlan_at," ", lslrPlan_atprev)

                    for i in 1:length(p.adaptoptions)
                        println(p.adaptoptions[i],v.surgeExposure[m,:],p.adaptoptions)
                        v.R[t, m, i] = calcHorR(-2, p.adaptoptions[i], lslrPlan_at, v.surgeExposure[m,:], p.adaptoptions)
                            
                        if is_first(t)
                            Rprev = calcHorR(-2, p.adaptoptions[i], p.lslr[1,m], v.surgeExposure[m,:], p.adaptoptions) 
                        else
                            Rprev = v.R[convert(Int,p.at[at_index_prev]),m, i]
                        end
    
    
                        v.FloodRetreat[t,m,i] = (p.tstep/atstep) * (atstep * v.landvalue[t,m]*.04 * calcCoastArea(v.areaparams[m,:], v.R[t,m,i]) +          
                            max(0,calcCoastArea(v.areaparams[m,:], v.R[t,m,i]) - calcCoastArea(v.areaparams[m,:], Rprev))* 
                            (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[t,m]) * 1e-4
    
                        v.RelocateRetreat[t,m,i] = (p.tstep / atstep) * 
                            max(0, calcCoastArea(v.areaparams[m,:], v.R[t,m,i]) - calcCoastArea(v.areaparams[m,:], Rprev)) * 
                            (p.movefactor * v.ypc_seg[t,m] * 1e-6 * v.popdens_seg[t,m] +
                            p.capmovefactor * p.mobcapfrac * v.capital[t,m] + p.democost * (1 - p.mobcapfrac ) * v.capital[t,m]) * 1e-4
           
    
                        if p.adaptoptions[i] >= 10
                            v.H[t,m, i-1] = calcHorR(-1, p.adaptoptions[i], lslrPlan_at, v.surgeExposure[m,:], p.adaptoptions)
    
                            if is_first(t)
                                Hprev = calcHorR(-1, p.adaptoptions[i], p.lslr[1,m], v.surgeExposure[m,:], p.adaptoptions)
                            else
                                Hprev = v.H[convert(Int,p.at[at_index_prev]),m,i-1]
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
                                    
                        end
    
                        for j in t_range
                            v.R[j,m,i] = v.R[t,m,i]
                            v.WetlandRetreat[j,m] = p.tstep * v.wetlandservice[j,rgn_ind] * v.wetlandloss[j,m] * min(v.coastArea[j,m], p.wetland[m])
    
                            v.StormRetreat[j,m,i] = p.tstep * (1 - v.ρ[j,rgn_ind]) * 
                                (p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, v.R[j,m,i] - p.lslr[j,m])))) * 
                                (v.capital[j,m] + v.popdens_seg[j,m] * v.vsl[j,m] * p.floodmortality)
    
                            v.FloodRetreat[j,m, i] = v.FloodRetreat[t,m i]
                            v.RelocateRetreat[j,m,i] = v.RelocateRetreat[t,m,i]
    
                            # Put all other costs intp $Billions from $M and divide by 10
                            v.StormRetreat[j,m,i]  = v.StormRetreat[j,m,i]  * 1e-4
                            v.WetlandRetreat[j,m] = v.WetlandRetreat[j,m] * 1e-4
                                        
                            v.RetreatCost[j,m, i] = v.FloodRetreat[j,m,i] + v.RelocateRetreat[j,m,i] + v.StormRetreat[j,m,i] + v.WetlandRetreat[j,m]
                                
                            if p.adaptoptions[i] >= 10
                                v.H[j,m, i-1] = v.H[t,m, i-1]
    
                                v.WetlandProtect[j,m] = p.tstep * p.wetland[m] .* v.wetlandservice[j,rgn_ind]
                                        
                                v.StormProtect[j,m,i-1] = p.tstep * (1 - v.ρ[j,rgn_ind]) * (p.psig0[m] + p.psig0coef[m] * max(0,p.lslr[j,m])) / 
                                                        (1. + p.psigA[m] * exp(p.psigB[m] * max(0,(v.H[j,m, i-1] - p.lslr[j,m])))) *
                                                        (v.capital[j,m] + v.popdens_seg[j,m] * v.vsl[j,m] * p.floodmortality)
                                    
                                v.Construct[j,m,i-1] = v.Construct[t,m, i-1]
    
                                # Put all other costs intp $Billions from $M and divide by 10
                                v.WetlandProtect[j,m] = v.WetlandProtect[j,m] * 1e-4
                                v.StormProtect[j,m,i-1] = v.StormProtect[j,m,i-1] * 1e-4
                                                    
                                v.ProtectCost[j,m,i-1] = v.Construct[j,m,i-1] + v.WetlandProtect[j,m] + v.StormProtect[j,m,i-1]
    
                            end
    
                        end
    
                        v.NPVRetreat[at_index,m,i] = sum([v.discountfactor[j] * v.RetreatCost[findind(j,t_range),m,i] for j in t_range])
    
                        if p.adaptoptions[i] >=10
                            v.NPVProtect[at_index,m,i-1] = sum( [ v.discountfactor[j] * v.ProtectCost[findind(j,t_range),m,i-1] for j in t_range] ) # Protect
                        end
                    end
    
                    # ** Choose Least Cost Option **
                    if is_first(t) && p.fixed
                        protectInd = indmin(v.NPVProtect[at_index,m,:])
                        retreatInd = indmin(v.NPVRetreat[at_index,m,:])
        
                        minLevels = [p.adaptoptions[protectInd+1], p.adaptoptions[retreatInd], 0]
                        choices = [v.NPVProtect[at_index,m,protectInd], v.NPVRetreat[at_index,m,retreatInd], v.NPVNoAdapt[at_index,m]]
                        leastcost = -1 * indmin(choices)
                        leastlevel = minLevels[indmin(choices)]
            
                        v.OptimalFixedOption[m] = leastcost
                        v.OptimalFixedLevel[m] = leastlevel
                            
                    end
                        
                    if v.OptimalFixedOption[m]==-1
                        v.OptimalFixedCost[t_range,m] = v.ProtectCost[t_range,m,(findall(i->i==v.OptimalFixedLevel[m], p.adaptoptions)-1)] 
                    elseif v.OptimalFixedOption[m]==-2
                        v.OptimalFixedCost[t_range,m] = v.RetreatCost[t_range,m, findall(i->i==v.OptimalFixedLevel[m], p.adaptoptions)] 
                    else
                        v.OptimalFixedCost[t_range,m] = v.NoAdaptCost[t_range,m]
                    end
    
                end
            end     
        end
    end
end


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
    return max(0, (lslr2 - lslr1)/tstep)
end

function getsegments(rgn_name, xsc)

    segs = collect(keys(filter( (k,v) -> v[1]==rgn_name,xsc)))
    return segs

end

function findind(val, vec)
    # Look up index corresponding to name
    # vec - a vector of region or segment names (strings)
    # val - a string corresponding to value in 'vec'
    h(u) = u == val
    name_ind = findall(h, vec)[1]
    return name_ind

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

function calcHorR(option, level, lslrPlan, surgeExpLevels, adaptoptions)

    ind = findind(level, adaptoptions)

    if option==-1 && level ==10
        # Protect height differs from retreat radius only in case of 10 yr surge exposure
        H = max(0, lslrPlan + surgeExpLevels[ind] / 2)

        return H
    else
        H_R = max(0, lslrPlan + surgeExpLevels[ind])
        return H_R
    end
end

function pos(x)
    return max(0,x)
end