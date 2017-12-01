# Catherine Ledna
# November 8, 2017
#
# CIAM (Coastal Impact and Adaptation Model)
# Adapted from Delavane Diaz (2016) 
#-------------------------------------------------------------------------------
# CIAM computes adaptation costs from sea level rise (slr) in a geographically 
# disaggregated cost-minimization framework. For details and documentation, 
# see Diaz (2016) and https://github.com/delavane/CIAM/.  
#
#-------------------------------------------------------------------------------

#------------
# TO DO: Fix pc / pc0 / cci issue in protect case -- Done
# TO DO: Update to time x region, not region x time matrices 
#       pop, ypcc done
#       to do: variables too for consistency (this will take awhile)
# TO DO: change names to fit params from ciam - mostly done?
# TO DO: implement full data importing paradigm (dictionary and all) on subset of segments
#       Deal with xsc and general region-segment issues 
# TO DO: Alter test code to reflect main code changes
# TO DO: Get results from GAMS 
# TO DO: Time everything 
# TO DO: Debug everything
#   Current status 11/17: 
#   Think I've changed everything. Changed areaparams back to 1-15 area; using variable container (may rethink this) 
#   Changed to time x region for at least 2 params in defcomp but haven't changed in code yet, to do
#   In process of changing surgeExposure to variable container  / removing level index (unnecessary) - Done 
#------------

using Mimi

@defcomp ciam begin
    # --- Indices ---
    regions = Index()
    segments = Index()                      # Q: Have all the segments upfront?
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
    
    popdens_seg = Variable(index = [segments, time])          # Population density of segment extrapolated forward in time (people / km^2)    
    ypc_seg = Variable(index = [segments, time])              # GDP per capita by segment ($2010 per capita) (multiplied by scaling factor)


    # ---Land Parameters---  
    landinput::Bool = Parameter()                   # Set to T for FUND or F for GTAP
         
    gtapland = Parameter( index = [regions])        # GTAP land value in 2007 (million 2010$ / km^2)
    dvbm = Parameter()                              # FUND value of OECD dryland per Darwin et al 1995 converted from $1995 ($2010M per sqkm) (5.376)
    kgdp = Parameter()                              # Capital output ratio (per MERGE) (3 by default) 
    discountrate = Parameter()                      # Discount rate (0.04 by default)
    depr = Parameter()                              # Fraction of capital that has not been depreciated over adaptation period (retreat cases)
    
    
    landdata = Variable( index = [regions])         # Takes on value of either fundland or gtapland
    fundland = Variable( index = [regions])         # FUND land value in 1995 (calculated in run_timestep) (million 2010$ / km^2), 
                                                    #   Q maybe import directly? 
    land_appr = Variable(index = [regions, time])   # Land appreciation rate (calculated as regression by Yohe ref Abraham and Hendershott) 
    coastland = Variable(index = [segments, time])  # Coastal land value (function of interior land value * scaling factor) ($2010M per sqkm)
    landvalue = Variable(index = [segments, time])  # Total endowment value of land ($2010M per sqkm)


    ρ = Variable(index = [regions, time])           # Country-wide resilience parameter (logistic function related to GDP)
    capital = Variable(index = [segments, time])    # Total endowment value of capital stock (million $2010 / km^2)
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
    vsl = Variable(index = [regions, time])     # Value of statistica life (million 2010$)
                                           
    # ---Wetland Loss Parameters---
    wbvm = Parameter()                                      # Annual value of wetland services (million 2010$ / km^2 / yr); (Brander et al 2006)  (0.376) 
    wetland = Parameter( index = [segments])            # Initial wetland area in coastal segment (km^2)
    wmaxrate = Parameter()                                  # Maximum rate of wetland accretion (m per yr) per Kirwan et al 2010 (0.01)
    
    wetlandservice = Variable(index = [regions, time])      # Annual value of wetland services TODO change doc
    wetlandloss = Variable(index = [segments, time])        # Fractional loss of wetland due to slr


    # ---Sea Level Rise Parameters---
    lslr = Parameter(index = [segments, time])                # Local sea level rise (m) 

    adaptOptions = Parameter(index = [5])                     # Index of available adaptation levels for protect and retreat (0 is no adaptation)
    s10 = Parameter(index = [segments])
    s100 = Parameter(index = [segments])
    s1000 = Parameter(index = [segments])
    smax = Parameter(index = [segments])
    surgeExposure::Float64 = Variable( index = [segments, 5])     # Storm surge exposure levels (corresponding to each designated adaptation option)
    

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

    coastArea = Variable(index=[segments, time])            # Calculate remaining coastal area after slr (m^2)
    
    # ---Intermediate Variables---
    R = Variable(index = [segments, adaptPers] )             # Matrix of retreat radii for each adaptation option/exposure level
    H = Variable(index = [segments, adaptPers] )             # Matrix of construction heights for each adaptation option

    # ---Decision Variables---   
    NoAdaptCost = Variable(index = [segments, time])         # Cost of not adapting (e.g. reactive retreat) (2010$)                                                
    ProtectCost = Variable(index = [segments, time])         # Cost of protection at each timestep
    RetreatCost = Variable(index = [segments, time])         # Cost of retreat at each timestep  

    # ---Outcome Variables---
    AdaptationDecision = Variable(index = [segments, time])   # Option chosen for adaptation period
    AdaptationCost = Variable(index = [segments, time])       # Cost of option chosen for adaptation period (B 2010$ / yr) 
    AdaptationLevel = Variable(index = [segments, time])      # Level of protect or retreat (if chosen)
    RegionalCost = Variable(index = [regions, time])          # Cost of adaptation at level of region

end

function run_timestep(s::ciam, t::Int)
    p = s.Parameters
    v = s.Variables
    d = s.Dimensions

    # In first period, initialize all non-adaptation dependent intermediate variables for all timesteps
    if t==1
      #  1. Initialize non-region dependent intermediate variables

        for i in collect(t:Int(p.ntsteps)) 
           v.discountfactor[i] = 1/(1 + p.discountrate)^(p.tstep * (i-1))
        end

        # Put area parameters into an array (using Variable structure to preserve across calls)
        v.areaparams = [p.area1 p.area2 p.area3 p.area4 p.area5 p.area6 p.area7 p.area8 p.area9 p.area10 p.area11 p.area12 p.area13 p.area14 p.area15]
        v.surgeExposure = [zeros(length(p.s10)) p.s10 p.s100 p.s1000 p.smax] 

        # 2. Initialize region-dependent intermediate variables
        for r in d.regions
            # Determine land input value (true = FUND, false = GTAP)
            # TODO switch to importing fund land values as param
            if p.landinput 
                v.fundland[r] = min(p.dvbm, max(0.005, p.dvbm * p.ypcc[1,r] * p.refpopdens[r] / (p.ypc_usa[1] * p.refpopdens[r])))
                v.landdata[r] = v.fundland[r]
            else
                v.landdata[r] = p.gtapland[r]
            end

            v.land_appr[r, t] = 1.
            v.wetlandservice[r, t] = p.wbvm * ((p.ypcc[t,r] / p.ypc_usa[1])^1.16 * (p.refpopdens[r] /27.59)^0.47) 
            v.ρ[r, t] = p.ypcc[t,r] / (p.ypcc[t,r] + p.ypc_usa[1])
            v.vsl[r, t] = 1e-6 * 216 * p.ypc_usa[t] * (p.ypcc[t,r]/p.ypc_usa[t])^0.05

            
            for i in collect(2:Int(p.ntsteps))
                v.land_appr[r, i] = v.land_appr[r, i-1] * exp(0.565 * growthrate(p.ypcc[i-1,r], p.ypcc[i,r]) + 0.313 * growthrate(p.pop[i-1,r], p.pop[i,r]))
                v.wetlandservice[r,i] = v.land_appr[r,i] * v.wetlandservice[r,1]
                v.ρ[r, i] = p.ypcc[i,r] / (p.ypcc[i,r] + p.ypc_usa[1]) 
                v.vsl[r, i] = 1e-6 * 216 * p.ypc_usa[i] * (p.ypcc[i,r]/p.ypc_usa[i])^0.05  
            end    
        end

        # 3. Initialize segment-dependent variables 
        for m in d.segments
            rgn_ind = getregion(m, p.xsc)
            # At beginning, initialize H and R at period 1 as 0. May get changed later when optimization is performed. 
            v.H[m, 1] = 0
            v.R[m, 1] = 0

            v.popdens_seg[m, t] = p.popdens[m]
            v.ypc_seg[m, t] = p.ypcc[t, rgn_ind] * max(0.9, (p.popdens[m]/250.)^0.05)
            v.capital[m, t] = p.kgdp * v.ypc_seg[m, t] * v.popdens_seg[m, t] * 1e-6
            v.coastland[m, t] = max(0.5, log(1+v.popdens_seg[m, t])/log(25)) * (v.land_appr[rgn_ind, t] * v.landdata[rgn_ind])  # Interior * scaling factor
            v.landvalue[m, t] = min(v.coastland[m, t], (v.land_appr[rgn_ind, t] * v.landdata[rgn_ind]))
            v.coastArea[m, t] = calcCoastArea(v.areaparams[m,:], p.lslr[m, t])  

            
            for i in collect(2:Int(p.ntsteps))
                v.popdens_seg[m,i] = v.popdens_seg[m, i-1] * (1 + growthrate(p.pop[ i-1, rgn_ind], p.pop[i, rgn_ind])) 
                v.ypc_seg[m, i] = p.ypcc[i,rgn_ind] * max(0.9, (p.popdens[m]/250.)^0.05) # ypcc * popdens scaling factor
                v.capital[m, i] = p.kgdp * v.ypc_seg[m, i] * v.popdens_seg[m, i] * 1e-6 
                v.coastland[m, i] = max(0.5, log(1+v.popdens_seg[m, i])/log(25)) * (v.land_appr[rgn_ind, i] * v.landdata[rgn_ind])
                v.landvalue[m, i] = min(v.coastland[m, i], (v.land_appr[rgn_ind, i] * v.landdata[rgn_ind]))
                v.coastArea[m, i] = calcCoastArea(v.areaparams[m,:], p.lslr[m, i])
                v.wetlandloss[m, i-1] = min(1, (localrate(p.lslr[m, i-1], p.lslr[m, i], p.tstep)/p.wmaxrate)^2)
            end

            v.wetlandloss[m, p.ntsteps] = min(1, (localrate(p.lslr[m, p.ntsteps-1], p.lslr[m, p.ntsteps], p.tstep)/p.wmaxrate)^2)  

        end

    end

    # # ------- Calculate adaptation decisions in each segment -------
    # # Only calculate at start of adaptation period, for all t in adaptation period
    if (t in p.at)           
        adapt_range = collect(1:length(p.adaptOptions))         

        # Determine length of adaptation period ("atstep")
        g(c) = c == t
        at_index = find(g, p.at)[1] # WEIRD - changing to g solved issue
        at_index_next = at_index + 1
        at_index_prev = max(1,at_index - 1)
        
        at_prev = Int(p.at[at_index_prev])

        if at_index_next <= length(p.at)
            atstep = (p.at[at_index_next] - p.at[at_index])*p.tstep   # years
            at_next = Int(p.at[at_index_next])  
            last_t = at_next-1
            last = 0
        else
            # Deal with special case of last adaptation period
            atstep = p.tstep*p.ntsteps - (p.at[at_index] * p.tstep) 
            at_next = p.ntsteps # Flag this assumes timesteps are indices not years (1:20 not 2010:2100)
            last_t = at_next
            last = 1
        end
        t_range = collect(t:last_t)
         

        # Calculate Adaptation Cost and Decision for each segment
        for m in d.segments
            if atstep==0
                v.AdaptationCost[m,t] = v.AdaptationCost[m, t-1]
                v.AdaptationDecision[m,t] = v.AdaptationDecision[m, t-1]
                v.AdaptationLevel[m, t] = v.AdaptationLevel[m, t-1]
            else
                rgn_ind = getregion(m, p.xsc)
                
                # ** Calculate No Adaptation Costs **

                for i in t_range
                    StormNoAdapt = p.tstep * (1 - v.ρ[rgn_ind , i]) * (p.rsig0[m] / (1 + p.rsigA[m] )) * 
                            (v.capital[m, i] + v.popdens_seg[m, i] * v.vsl[rgn_ind, i] * p.floodmortality)
                        
                    WetlandNoAdapt = p.tstep * v.wetlandservice[rgn_ind, i] * v.wetlandloss[m, i] * 
                            min(v.coastArea[m, i], p.wetland[m])

                    if i==p.ntsteps
                        FloodNoAdapt = p.tstep * v.landvalue[m,i-1]*.04 * max(0, v.coastArea[m, i]) + (max(0, v.coastArea[m, i]) - max(0, v.coastArea[m, i-1])) * 
                            (1 - p.mobcapfrac) * v.capital[m, i-1]
    
                        RelocateNoAdapt = (max(0, v.coastArea[m, i]) - max(0,v.coastArea[m, i-1])) * (5 * p.movefactor * p.ypcc[i-1, rgn_ind]*1e-6*v.popdens_seg[m, i-1] +
                            p.capmovefactor * p.mobcapfrac * v.capital[m, i-1] + p.democost * (1 - p.mobcapfrac) * v.capital[m, i-1])
                    else
                        FloodNoAdapt = p.tstep * v.landvalue[m,i]*.04 * max(0, v.coastArea[m, i+1]) + (max(0, v.coastArea[m, i+1]) - max(0, v.coastArea[m,i])) * 
                            (1 - p.mobcapfrac) * v.capital[m,i]                    
                       
                        RelocateNoAdapt = (max(0, v.coastArea[m,i+1]) - max(0,v.coastArea[m,i])) * (5 * p.movefactor * p.ypcc[i, rgn_ind]*1e-6*v.popdens_seg[m,i] +
                            p.capmovefactor * p.mobcapfrac * v.capital[m,i] + p.democost * (1 - p.mobcapfrac) * v.capital[m,i])
                    end
                                                    
                    v.NoAdaptCost[m,i] = WetlandNoAdapt + FloodNoAdapt + RelocateNoAdapt + StormNoAdapt

                end
                
                NPVNoAdapt = sum( [ v.discountfactor[j] * v.NoAdaptCost[m,j] for j in t_range] )

               
                # ** Calculate Protect and Retreat Costs ** 
                lslrPlan_at = p.lslr[m, at_next]

                if (t > 1 && p.fixed)

                    option = v.AdaptationDecision[m, at_prev]
                    level = v.AdaptationLevel[m, at_prev]

                                                                    # TODO 1/10 index for protect
                    if option==-1 # Protect
                        v.H[m,at_index] = calcHorR(option, level, lslrPlan_at, v.surgeExposure[m,:], p.adaptOptions)
    
                        Construct = (p.tstep/atstep) * p.length[m] * p.pc0 * p.cci[rgn_ind] * (p.pcfixed + (1- p.pcfixed)*(v.H[m,at_index]^2 - v.H[m, at_index_prev]^2) + 
                                p.mc*atstep*v.H[m,at_index]) + p.length[m] * 1.7 * v.H[m,at_index] * v.landvalue[m,t]*.04/2*atstep

                        for i in t_range
                            WetlandProtect = p.tstep * p.wetland[m] * v.wetlandservice[rgn_ind, i]
                            StormProtect = p.tstep * (1 - v.ρ[rgn_ind, i]) * (p.psig0[m] + p.psig0coef[m] * p.lslr[m, i]) / (1. + p.psigA[m] * exp(p.psigB[m] * max(0,(v.H[m,at_index] - p.lslr[m,i])))) *
                                                        (v.capital[m,i] + v.popdens_seg[m,i] * v.vsl[rgn_ind, i] * p.floodmortality)
    
                            v.AdaptationDecision[m, i] = v.AdaptationDecision[m, at_prev]
                            v.AdaptationLevel[m, i] =  v.AdaptationLevel[m, at_prev]
                            v.AdaptationCost[m, i] = (Construct + WetlandProtect + StormProtect) .* 1e-6 ./ p.tstep
                        end
                        
                    elseif option==-2 # Retreat
                        v.R[m,at_index] = calcHorR(option, level, lslrPlan_at, v.surgeExposure[m,:], p.adaptOptions)
           
                        RelocateRetreat = (p.tstep / atstep) * max(0, calcCoastArea(v.areaparams[m,:], v.R[m,at_index]) - calcCoastArea(v.areaparams[m,:], v.R[m,at_index_prev])) * 
                                                p.movefactor * v.ypc_seg[m,t] * 1e-6 * v.popdens_seg[m,t] +
                                                p.capmovefactor * p.mobcapfrac * v.capital[m,t] + p.democost * (1 - p.mobcapfrac ) * v.capital[m,t]
    
                                                
                        FloodRetreat = (p.tstep/atstep) * atstep * v.landvalue[m,t]*.04 * calcCoastArea(v.areaparams[m,:], v.R[m,at_index]) + 
                                             max(0,calcCoastArea(v.areaparams[m,:], v.R[m,at_index]) - calcCoastArea(v.areaparams[m,:], v.R[m,at_index_prev]))* 
                                             (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[m,t] # todo check this summation
                        
                        for i in t_range
                            WetlandRetreat = p.tstep * v.wetlandservice[rgn_ind, i] * v.wetlandloss[m, i] * 
                                min(v.coastArea[m, i], p.wetland[m])

                            StormRetreat = p.tstep * (1 - v.ρ[rgn_ind, i]) * 
                                (p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, v.R[m, at_index] - p.lslr[m, i])))) * 
                                (v.capital[m, i] + v.popdens_seg[m, i] * v.vsl[rgn_ind, i] * p.floodmortality) 
                                
                            v.AdaptationDecision[m, i] = v.AdaptationDecision[m, at_prev]
                            v.AdaptationLevel[m, i] =  v.AdaptationLevel[m, at_prev]
                            v.AdaptationCost[m, i] = (FloodRetreat + RelocateRetreat + StormRetreat + WetlandRetreat) .* 1e-6 ./ p.tstep
                        end
    
                    else # No Adapt
                        for i in t_range
                            v.AdaptationDecision[m, i] = v.AdaptationDecision[m, at_prev]
                            v.AdaptationLevel[m, i] =  v.AdaptationLevel[m, at_prev]
                            v.AdaptationCost[m, i] = v.NoAdaptCost[m,i] .* 1e-6 ./ p.tstep
                        end
                    end
                
                # ** Non-Fixed Option (Choose Least Cost) **
                # ** Calculate both protect and retreat and choose option based on least cost NPV
                #      if in first period or not using fixed model **
                else

                    # ** Calculate Protect and Retreat Cost and NPV by level **
                    # ** Initialize Temp Variables ** 
                    NPVAdapt = NaN .* ones(2, length(p.adaptOptions))
                    RetreatTot = NaN .* ones(length(p.adaptOptions), length(t_range))
                    ProtectTot = NaN .* ones(length(p.adaptOptions), length(t_range))
    

                    for i in collect(1:length(p.adaptOptions))
                        # ** Calculate Retreat Costs by Adaptation Level **                   
                        R = calcHorR(-2, p.adaptOptions[i], lslrPlan_at, v.surgeExposure[m,:], p.adaptOptions)
    
                        FloodRetreat = (p.tstep/atstep) * atstep * v.landvalue[m,t]*.04 * calcCoastArea(v.areaparams[m,:], R) + 
                            max(0,calcCoastArea(v.areaparams[m,:], R) - calcCoastArea(v.areaparams[m,:], v.R[m,at_index_prev]))* 
                            (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[m,t] # todo check this summation
                
                            # check this as well
                        RelocateRetreat = (p.tstep / atstep) * max(0, calcCoastArea(v.areaparams[m,:], R) - calcCoastArea(v.areaparams[m,:], v.R[m,at_index_prev])) * 
                            p.movefactor * v.ypc_seg[m,t] * 1e-6 * v.popdens_seg[m,t] +
                            p.capmovefactor * p.mobcapfrac * v.capital[m,t] + p.democost * (1 - p.mobcapfrac ) * v.capital[m,t]
                        
                        for j in t_range
                            WetlandRetreat = p.tstep * v.wetlandservice[rgn_ind, i] * v.wetlandloss[m, i] * 
                                min(v.coastArea[m, i], p.wetland[m])

                            StormRetreat = p.tstep * (1 - v.ρ[rgn_ind, j]) * 
                                (p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, R - p.lslr[m, j])))) * 
                                (v.capital[m, j] + v.popdens_seg[m, j] * v.vsl[rgn_ind, j] * p.floodmortality) 
    
                            RetreatTot[i,j] = FloodRetreat + RelocateRetreat + StormRetreat + WetlandRetreat
                        end
                           
                        NPVAdapt[2, i] = sum([v.discountfactor[j] * RetreatTot[i,j] for j in t_range])
                        
                        # ** Calculate protect costs for a subset of adaptation levels ** 
                        if p.adaptOptions[i] >= 10
                            H = calcHorR(-1, p.adaptOptions[i], lslrPlan_at, v.surgeExposure[m,:], p.adaptOptions)
                            
                            # Flexible mode question: v.H[m, at_prev] for t > 1. 
                            Construct = (p.tstep/atstep) * 
                                p.length[m] * p.pc0 * p.cci[rgn_ind] * (p.pcfixed + (1- p.pcfixed)*(H^2 - v.H[m, at_index_prev]^2) + 
                                p.mc*atstep*H) + p.length[m] * 1.7 * H * v.landvalue[m,t]*.04/2*atstep
    
                            for j in t_range
                                WetlandProtect = p.tstep * p.wetland[m] .* v.wetlandservice[rgn_ind, j]
                                
                                StormProtect = p.tstep * (1 - v.ρ[rgn_ind, j]) * (p.psig0[m] + p.psig0coef[m] * p.lslr[m, j]) / 
                                                                (1. + p.psigA[m] * exp(p.psigB[m] * max(0,(H - p.lslr[m,j])))) *
                                                                (v.capital[m,j] + v.popdens_seg[m,j] * v.vsl[rgn_ind, j] * p.floodmortality)
                                        
                                ProtectTot[i,j] = Construct + WetlandProtect + StormProtect
                
                            end
    
                            NPVAdapt[1,i] = sum( [ v.discountfactor[j] * ProtectTot[i,j] for j in t_range] ) # Protect
                            
                        end           
                    end
    
                    # ** Choose least cost adaptation option **

                    protectInd = indmin(NPVAdapt[1,:])
                    retreatInd = indmin(NPVAdapt[2,:])
    
                    minLevels = [p.adaptOptions[protectInd], p.adaptOptions[retreatInd], 0]
                    choices = [NPVAdapt[1,protectInd], NPVAdapt[2,retreatInd], NPVNoAdapt]
                    leastcost = -1 * indmin(choices)
                    leastlevel = minLevels[indmin(choices)]
    
                    v.AdaptationDecision[m, t_range] = leastcost
                    v.AdaptationLevel[m, t_range] = leastlevel
                        
                    # Set H or R
                    v.H[m, at_index] = 0
                    v.R[m, at_index] = 0

                    if leastcost==-1
                        v.AdaptationCost[m, t_range] = ProtectTot[protectInd, t_range] .* 1e-6 ./ p.tstep                        
                        v.H[m, at_index] =  calcHorR(-1, p.adaptOptions[protectInd], lslrPlan_at, v.surgeExposure[m,:], p.adaptOptions)
                    elseif leastcost==-2
                        v.AdaptationCost[m, t_range] = RetreatTot[retreatInd, t_range] .* 1e-6 ./ p.tstep 
                        v.R[m, at_index] = calcHorR(-2, p.adaptOptions[retreatInd], lslrPlan_at, v.surgeExposure[m,:], p.adaptOptions)
                    end
    
                end
            end

            
        end
        # ** Sum Costs to Regional Level **
        for r in d.regions
            segs = getsegments(r, p.xsc)
            v.RegionalCost[r, t_range] = sum( [v.AdaptationCost[m, t_range] for m in segs ])
        end

    end
end


    
# Helper functions
function growthrate(x1, x2)
    epsilon = 1.e-9
    return (x2 / (x1 + epsilon) - 1.)
end

function calcCoastArea(areaparams, var)

    area = (areaparams[1] * max(0, min(0.5, var)) + areaparams[1] + areaparams[2] )/ 2 * max(0, min(1, var-0.5)) + 
        areaparams[2] * max(0, min(0.5, var - 1.5)) + areaparams[3] * max(0, min(1, var-2)) + areaparams[4] * max(0, min(1, var-3)) + 
        areaparams[5] * max(0, min(1, var - 4)) + areaparams[6] * max(0, min(1, var - 5)) + areaparams[7] * max(0, min(1, var-6)) + 
        areaparams[8] * max(0, min(1, var - 7)) + areaparams[9] * max(0, min(1, var - 8)) + areaparams[10] * max(0, min(1, var - 9)) +
        areaparams[11] * max(0, min(1, var - 10)) + areaparams[12] * max(0, min(1, var-11)) + areaparams[13] * max(0, min(1, var-12)) + 
        areaparams[14] * max(0, min(1, var - 13)) + areaparams[15] * max(0, var - 14) 

    return area

end

function localrate(lslr1, lslr2, tstep)
    return max(0, (lslr2 - lslr1)/tstep)
end

function getsegments(rgn_name, xsc)

    segs = collect(keys(filter( (k,v) -> v==rgn_name,xsc)))
    return segs

end

function findind(val, vec)
    # Look up index corresponding to name
    # vec - a vector of region or segment names (strings)
    # val - a string corresponding to value in 'vec'
    h(u) = u == val
    name_ind = find(h, vec)[1]
    return name_ind

end

function getregion(seg_ind, xsc, index=true)
    rgn = xsc[seg_ind]
    
  #  if index
  #      rgn_ind = index_lookup(rgn, rgnms)
  #      return rgn_ind
  #  end

    return rgn
end

function calcHorR(option, level, lslrPlan, surgeExpLevels, adaptOptions)

    ind = findind(level, adaptOptions)

    if option==-1 && level ==10
        # Protect height differs from retreat radius only in case of 10 yr surge exposure
        H = max(0, lslrPlan + surgeExpLevels[ind] / 2)

        return H
    else
        H_R = max(0, lslrPlan + surgeExpLevels[ind])
        return H_R
    end
end
