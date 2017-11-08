# Catherine Ledna
# October 6, 2017
#
# CIAM (Coastal Impact and Adaptation Model)
# Adapted from Delavane Diaz (2016) 
#-------------------------------------------------------------------------------
# CIAM computes adaptation costs from sea level rise (slr) in a geographically 
# disaggregated cost-minimization framework. For details and documentation, 
# see Diaz (2016) and https://github.com/delavane/CIAM/.  
#
#-------------------------------------------------------------------------------

using Mimi

@defcomp ciam begin
    # --- Indices ---
    regions = Index()
    segments = Index()                      # Q: Have all the segments upfront?
    level = Index()
    adaptPers = Index()
   
    # --- Region / segment mapping ---
   # rgn_names::Array{String,1} = Parameter()
   # seg_names::Array{String,1} = Parameter()
    xsc::Dict{Int64, Int64} = Parameter()        # Region to segment mapping    

    # ---Time-related Parameters---
    tstep = Parameter()                     # Length of individual time-step (years)
    at = Parameter( index = [adaptPers])    # Array of time indices that mark starts of adaptation periods 
    ntsteps::Int = Parameter()              # Number of time-steps     

    # ---Model Parameters ---
    fixed::Bool = Parameter()               # Run model as fixed (T) or flexible (F) with respect to adaptation

    # ---Socioeconomic Parameters---
    pop_country = Parameter(index = [regions, time])          # Population of country (million people) (from MERGE)
    refpopdens_country = Parameter( index = [regions])        # Reference population density of country (people / km^2)
    refpopdens_usa = Parameter()                              # Reference population density of USA (people/km^2) 
    popdens1_seg = Parameter( index = [segments])             # Pop density of segment in time t = 1 (people/km^2)
    popdens = Variable(index = [segments, time])              # Population density of segment extrapolated forward in time (people / km^2)

    ypc_country = Parameter(index = [regions, time]) # GDP per capita per country ($2010 per capita)
    ypc_usa = Parameter(index = [time])     # GDP per capita in USA; used as benchmark ($2010 per capita)
    ypc_seg = Variable(index = [segments, time])      # GDP per capita by segment ($2010 per capita) (multiplied by scaling factor)


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
    interior = Variable(index = [regions, time])    # Value of interior land (function of land value and appreciation rate) ($2010M per sqkm)
    coastland = Variable(index = [segments, time])  # Coastal land value (function of interior land * scaling factor) ($2010M per sqkm)
    landvalue = Variable(index = [segments, time])  # Total endowment value of land ($2010M per sqkm)


    ρ = Variable(index = [regions, time])           # Country-wide resilience parameter (logistic function related to GDP)
    capital = Variable(index = [segments, time])    # Total endowment value of capital stock (million $2010 / km^2)
    discountfactor = Variable(index=[time])         # Discount factor (derived from discount rate)
    
    # ---Coastal Parameters---
    length = Parameter(index = [segments])          # Segment length (km)

    # ---Protection Parameters---
    pc = Parameter(index = [segments])      # Cost of protection per segment (million 2010$ / km / vertical m^2 )
    pcfixed = Parameter()                   # Fraction of protection cost that is fixed (not variable in height) (0.3)
    mc = Parameter()                        # Maintenance cost (Hillen et al, 2010) (2%/yr) 
    pc0 = Parameter()                       # Reference cost of protection (million 2010$ / km / vert m^2) (6.02 by default)


    # ---Retreat /No Adapt Parameters---
    mobcapfrac = Parameter()                # Fraction of capital that is mobile (0.25)
    movefactor = Parameter()                # Cost to relocate people as a factor of annual income (Tol 3x RMendelsohn 0.5x) (1)
    capmovefactor = Parameter()             # Cost to relocate mobile capital as a fraction of asset value (0.1)
    democost = Parameter()                  # Cost to demolish immobile capital as fraction of asset (0.05)

    # # ---Surge Exposure Parameters---
    # Protection case
    pσ₀ = Parameter( index = [segments])
    pσ₀coef = Parameter(index = [segments])
    pσ₁ = Parameter(index = [segments])                       # psigA in GAMS code
    pσ₂ = Parameter(index = [segments])                       # psigB in GAMS code

    # Retreat / No Adapt Cases
    rσ₀ = Parameter( index = [segments])
    rσ₁ = Parameter( index = [segments])                       # rsigA in GAMS code
    rσ₂ = Parameter( index = [segments])                       # rsigB in GAMS code
    

    # ---Storm damage parameters---
    floodmortality = Parameter()                # Flood deaths as percent of exposed population; (Jonkman Vrijling 2008) (0.01) 
    vsl = Variable(index = [regions, time])     # Value of statistica life (million 2010$)
                                           
    # ---Wetland Loss Parameters---
    wbvm = Parameter()                                      # Annual value of wetland services (million 2010$ / km^2 / yr); (Brander et al 2006)  (0.376) 
    wetlandarea = Parameter( index = [segments])            # Initial wetland area in coastal segment (km^2)
    wmaxrate = Parameter()                                  # Maximum rate of wetland accretion (m per yr) per Kirwan et al 2010 (0.01)
    
    wetlandservice = Variable(index = [regions, time])      # Annual value of wetland services TODO change doc
    wetlandloss = Variable(index = [segments, time])        # Fractional loss of wetland due to slr
    exposedwetlandarea = Variable(index = [segments, time]) # Total exposed wetland area in timestep

    # ---Sea Level Rise Parameters---
    lslr = Parameter(index = [segments, time])                # Local sea level rise (m) 

    adaptOptions = Parameter( index = [level])      # Index of available adaptation levels for protect and retreat (0 is no adaptation)
    surgeExposure::Float64 = Parameter( index = [segments, level])     # Storm surge exposure levels (for each designated adaptation option)
    
    R = Variable( index = [segments, adaptPers] )   # Matrix of retreat radii for each adaptation option/exposure level
    H = Variable( index = [segments, adaptPers] )   # Matrix of construction heights for each adaptation option

    # ---Coastal Area Parameters---
    areaparams = Parameter( index = [segments, 15])         # Parameters used in calculating area of coast that is inundated
                                                            #   under retreat or no-adaptation scenarios

    coastArea = Variable(index=[segments, time])                      # Calculate remaining coastal area after slr (m^2)
    
    # ---Intermediate Variables---
    # No Adaptation Case
    FloodNoAdapt = Variable( index = [segments, time] )
    RelocateNoAdapt = Variable( index = [segments, time] )
    WetlandRetreatNoAdapt = Variable( index = [segments, time] )
    StormNoAdapt = Variable( index = [segments, time] )

    # Protect Case
    StormProtect = Variable( index = [segments, time] )
    WetlandProtect = Variable( index = [segments, time] )          # Wetland costs not contingent on adaptation level
    
    # Retreat Case
    StormRetreat = Variable( index = [segments, time] )
    FloodRetreat = Variable( index = [segments, time] )
    RelocateRetreat = Variable( index = [segments, time] )

    # ---Decision Variables---   
    NoAdaptCost = Variable( index = [segments, time])       # Cost of not adapting (e.g. reactive retreat) (2010$)                                                
    ProtectCost = Variable( index = [segments, time])          # Cost of protection at each timestep
    RetreatCost = Variable( index = [segments, time])          # Cost of retreat at each timestep  


    # ---Outcome Variables---
    AdaptationDecision = Variable(index = [segments, time])   # Option chosen for adaptation period
    AdaptationCost = Variable(index = [segments, time])       # Cost of option chosen for adaptation period (discretized according to tstep) 
    AdaptationLevel = Variable(index = [segments, time])      # Level of protect or retreat (if chosen)

end

function run_timestep(s::ciam, t::Int)
    p = s.Parameters
    v = s.Variables
    d = s.Dimensions

    print("\nstarting\n")
    # Define some internal params that are used for all cases
    protectOptions = p.adaptOptions[ p.adaptOptions[:] .>= 10 ]

    # In first period, initialize all non-adaptation dependent intermediate variables for all timesteps
    if t==1
      #  1. Initialize non-region dependent intermediate variables
        for i in collect(t:Int(p.ntsteps)) 
           v.discountfactor[i] = 1/(1 + p.discountrate)^(p.tstep * (i-1))
        end

        # 2. Initialize region-dependent intermediate variables
        for r in d.regions
            # Determine land input value (true = FUND, false = GTAP)
            if p.landinput 
                v.fundland[r] = min(p.dvbm, max(0.005, p.dvbm * p.ypc_country[r,1] * p.refpopdens_country[r] / (p.ypc_usa[1] * p.refpopdens_country[r])))
                v.landdata[r] = v.fundland[r]
            else
                v.landdata[r] = p.gtapland[r]
            end

            v.land_appr[r, t] = 1.
            v.wetlandservice[r, t] = p.wbvm * ((p.ypc_country[r, t] / p.ypc_usa[1])^1.16 * (p.refpopdens_country[r] /27.59)^0.47) 
            v.ρ[r, t] = p.ypc_country[r, t] / (p.ypc_country[r, t] + p.ypc_usa[1])
            v.interior[r, t] = v.land_appr[r, t] * v.landdata[r] 
            v.vsl[r, t] = 1e-6 * 216 * p.ypc_usa[t] * (p.ypc_country[r, t]/p.ypc_usa[t])^0.05
            

            for i in collect(2:Int(p.ntsteps))
                v.land_appr[r, i] = v.land_appr[r, i-1] * exp(0.565 * growthrate(p.ypc_country[r,i-1], p.ypc_country[r,i]) + 0.313 * growthrate(p.pop_country[r,i-1], p.pop_country[r,i]))
                v.wetlandservice[r,i] = v.land_appr[r,i] * v.wetlandservice[r,1]
                v.ρ[r, i] = p.ypc_country[r, i] / (p.ypc_country[r, i] + p.ypc_usa[1]) 
                v.interior[r, i] = v.land_appr[r, i] * v.landdata[r] 
                v.vsl[r, i] = 1e-6 * 216 * p.ypc_usa[i] * (p.ypc_country[r, i]/p.ypc_usa[i])^0.05  
                
            end    
        end

        # 3. Initialize segment-dependent variables 
        print("\nsegs: ", d.segments)
        for m in d.segments
            rgn_ind = getregion(m, p.xsc)
            # At beginning, initialize H and R at period 1 as 0. May get changed later when optimization is performed. 
            v.H[m, 1] = 0
            v.R[m, 1] = 0
     
            v.popdens[m, t] = p.popdens1_seg[m]
            v.ypc_seg[m, t] = p.ypc_country[rgn_ind, t] * max(0.9, (p.popdens1_seg[m]/250.)^0.05)
            v.capital[m, t] = p.kgdp * v.ypc_seg[m, t] * v.popdens[m, t] * 1e-6
            v.coastland[m, t] = max(0.5, log(1+v.popdens[m, t])/log(25)) * v.interior[rgn_ind, t] # Interior * scaling factor
            v.landvalue[m, t] = min(v.coastland[m, t], v.interior[rgn_ind, t])
            v.coastArea[m, t] = calcCoastArea(p.areaparams[m,:], p.lslr[m, t])  
            
            for i in collect(2:Int(p.ntsteps))
                v.popdens[m,i] = v.popdens[m, i-1] * (1 + growthrate(p.pop_country[rgn_ind, i-1], p.pop_country[rgn_ind, i])) 
                v.ypc_seg[m, i] = p.ypc_country[rgn_ind, i] * max(0.9, (p.popdens1_seg[m]/250.)^0.05) # ypc_country * popdens scaling factor
                v.capital[m, i] = p.kgdp * v.ypc_seg[m, i] * v.popdens[m, i] * 1e-6 
                v.coastland[m, i] = max(0.5, log(1+v.popdens[m, i])/log(25)) * v.interior[rgn_ind, i]
                v.landvalue[m, i] = min(v.coastland[m, i], v.interior[rgn_ind, i])
                v.coastArea[m, i] = calcCoastArea(p.areaparams[m,:], p.lslr[m, i])
                v.wetlandloss[m, i-1] = min(1, (localrate(p.lslr[m, i-1], p.lslr[m, i], p.tstep)/p.wmaxrate)^2)
            end

            v.wetlandloss[m, p.ntsteps] = min(1, (localrate(p.lslr[m, p.ntsteps-1], p.lslr[m, p.ntsteps], p.tstep)/p.wmaxrate)^2)  

        end
    end

    # ------- Calculate adaptation decisions in each segment -------
    # Only calculate at start of adaptation period, for all t in adaptation period
    print("\nhere")
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
            print("\nseg: ",m)
            rgn_ind = getregion(m, p.xsc)
            
            # ** Calculate No Adaptation Costs **
            print("\ntrange: ", t_range)
            for i in t_range
                v.StormNoAdapt[m, i] = p.tstep * (1 - v.ρ[rgn_ind , i]) * (p.rσ₀[m] / (1 + p.rσ₁[m] )) * 
                                         (v.capital[m, i] + v.popdens[m, i] * v.vsl[rgn_ind, i] * p.floodmortality)
                
                v.WetlandRetreatNoAdapt[m, i] = p.tstep * v.wetlandservice[rgn_ind, i] * v.wetlandloss[m, i] * min(v.coastArea[m, i], p.wetlandarea[m])


            end
            if last==1
                for i in collect(t:last_t - 1)
                    v.FloodNoAdapt[m, i] = p.tstep * v.landvalue[m,i]*.04 * max(0, v.coastArea[m, i+1]) + (max(0, v.coastArea[m, i+1]) - max(0, v.coastArea[m, i])) * 
                                         (1 - p.mobcapfrac) * v.capital[m, i]                    
                        
                    v.RelocateNoAdapt[m, i] = (max(0, v.coastArea[m, i+1]) - max(0,v.coastArea[m, i])) * (5 * p.movefactor * p.ypc_country[rgn_ind, i]*1e-6*v.popdens[m, i] +
                                         p.capmovefactor * p.mobcapfrac * v.capital[m, i] + p.democost * (1 - p.mobcapfrac) * v.capital[m, i])
                        
                    v.NoAdaptCost[m, i] = v.WetlandRetreatNoAdapt[m, i] + v.FloodNoAdapt[m, i] + v.RelocateNoAdapt[m, i] + v.StormNoAdapt[m, i]
                end
                
                v.FloodNoAdapt[m, last_t] = v.FloodNoAdapt[m, last_t - 1]
                v.RelocateNoAdapt[m, last_t] = v.RelocateNoAdapt[m, last_t - 1]
                v.NoAdaptCost[m, last_t] = v.NoAdaptCost[m, last_t - 1]
                NPVNoAdapt = sum( [ v.discountfactor[j] * v.NoAdaptCost[m,j] for j in t_range] )
            else
                for i in t_range
                    v.FloodNoAdapt[m, i] = p.tstep * v.landvalue[m,i]*.04 * max(0, v.coastArea[m, i+1]) + (max(0, v.coastArea[m, i+1]) - max(0, v.coastArea[m,i])) * 
                                         (1 - p.mobcapfrac) * v.capital[m,i]                    
                                    
                    v.RelocateNoAdapt[m,i] = (max(0, v.coastArea[m,i+1]) - max(0,v.coastArea[m,i])) * (5 * p.movefactor * p.ypc_country[rgn_ind,i]*1e-6*v.popdens[m,i] +
                                         p.capmovefactor * p.mobcapfrac * v.capital[m,i] + p.democost * (1 - p.mobcapfrac) * v.capital[m,i])
                        
                    v.NoAdaptCost[m,i] = v.WetlandRetreatNoAdapt[m,i] + v.FloodNoAdapt[m,i] + v.RelocateNoAdapt[m,i] + v.StormNoAdapt[m,i]
                
                end
                NPVNoAdapt = sum( [ v.discountfactor[j] * v.NoAdaptCost[m,j] for j in t_range] )

            end

            # ** Calculate Protect and Retreat Costs ** 
            lslrPlan_at = p.lslr[m, at_next]
            print("\nlslrPlan: ", lslrPlan_at)

            if (t > 1 && p.fixed)
                print("\nStarting fixed: ", t)
                option = v.AdaptationDecision[m, at_prev]
                level = v.AdaptationLevel[m, at_prev]
                print("\n", option, "\nlevel ", level)
                                                                # TODO 1/10 index for protect
                if option==-1 # Protect
                    v.H[m,at_index] = calcHorR(option, level, lslrPlan_at, p.surgeExposure[m,:], p.adaptOptions)

                    Construct = (p.tstep/atstep) * 
                            sum([p.length[m] * p.pc[m] * (p.pcfixed + (1- p.pcfixed)*(v.H[m,t]^2 - v.H[m, at_index_prev]^2) + 
                            p.mc*atstep*v.H[m,at_index]) + p.length[m] * 1.7 * v.H[m,at_index] * v.landvalue[m,j]*.04/2*atstep for j in t_range])
                    print("Construct", Construct)
                    for i in t_range
                        v.WetlandProtect[m, i] = p.tstep * p.wetlandarea[m] * v.wetlandservice[rgn_ind, i]
                        v.StormProtect[m, i] = p.tstep * (1 - v.ρ[rgn_ind, i]) * (p.pσ₀[m] + p.pσ₀coef[m] * p.lslr[m, i]) / (1. + p.pσ₁[m] * exp(p.pσ₂[m] * max(0,(v.H[m,at_index] - p.lslr[m,i])))) *
                                                    (v.capital[m,i] + v.popdens[m,i] * v.vsl[rgn_ind, i] * p.floodmortality)

                        v.AdaptationDecision[m, i] = v.AdaptationDecision[m, at_prev]
                        v.AdaptationLevel[m, i] =  v.AdaptationLevel[m, at_prev]
                        v.AdaptationCost[m, i] = Construct + v.WetlandProtect[m, i] + v.StormProtect[m, i]
                    end
                    
                elseif option==-2 # Retreat
                    v.R[m,at_index] = calcHorR(option, level, lslrPlan_at, p.surgeExposure[m,:], p.adaptOptions)
                    print("\nR ", v.R)                  
                    v.RelocateRetreat[m,t] = (p.tstep / atstep) * max(0, calcCoastArea(p.areaparams[m,:], v.R[m,at_index]) - calcCoastArea(p.areaparams[m,:], v.R[m,at_index_prev])) * 
                                            p.movefactor * v.ypc_seg[m,t] * 1e-6 * v.popdens[m,t] +
                                            p.capmovefactor * p.mobcapfrac * v.capital[m,t] + p.democost * (1 - p.mobcapfrac ) * v.capital[m,t]

                                            
                    v.FloodRetreat[m,t] = (p.tstep/atstep) * atstep * v.landvalue[m,t]*.04 * calcCoastArea(p.areaparams[m,:], v.R[m,at_index]) + 
                                         max(0,calcCoastArea(p.areaparams[m,:], v.R[m,at_index]) - calcCoastArea(p.areaparams[m,:], v.R[m,at_index_prev]))* 
                                         (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[m,t] # todo check this summation
                    print("Rel ", v.RelocateRetreat[m,t], "\nFlood ", v.FloodRetreat[m,t])

                    for i in t_range
                        v.StormRetreat[m,i] = p.tstep * (1 - v.ρ[rgn_ind, i]) * 
                            (p.rσ₀[m] / (1 + p.rσ₁[m] * exp(p.rσ₂[m] * max(0, v.R[m, at_index] - p.lslr[m, i])))) * 
                            (v.capital[m, i] + v.popdens[m, i] * v.vsl[rgn_ind, i] * p.floodmortality) 
                            
                        v.AdaptationDecision[m, i] = v.AdaptationDecision[m, at_prev]
                        v.AdaptationLevel[m, i] =  v.AdaptationLevel[m, at_prev]
                        v.AdaptationCost[m, i] = v.FloodRetreat[m,t] + v.RelocateRetreat[m,t] + v.StormRetreat[m, i] + v.WetlandRetreatNoAdapt[m,i]
                    end

                else # No Adapt
                    for i in t_range
                        v.AdaptationDecision[m, i] = v.AdaptationDecision[m, at_prev]
                        v.AdaptationLevel[m, i] =  v.AdaptationLevel[m, at_prev]
                        v.AdaptationCost[m, i] = v.NoAdaptCost[m,i]
                    end
                end
            
            # ** Non-Fixed Option (Choose Least Cost) **
            # ** Calculate both protect and retreat and choose option based on least cost NPV
            #      if in first period or not using fixed model **
            else
                print("\nStartingt=1\n")
                # ** Calculate Protect and Retreat Cost and NPV by level **
                v.WetlandProtect[m, t_range] = p.tstep * p.wetlandarea[m] .* v.wetlandservice[rgn_ind, t_range]
    
                # ** Initialize Temp Variables ** 
                StormRet = NaN .* ones(length(p.adaptOptions),length(t_range))   
                RelocateRet = NaN .* ones(length(p.adaptOptions))
                FloodRet = NaN .* ones(length(p.adaptOptions))         

                NPVAdapt = NaN .* ones(2, length(p.adaptOptions))
     
                RetreatTot = NaN .* ones(length(p.adaptOptions), length(t_range))
                ProtectTot = NaN .* ones(length(p.adaptOptions), length(t_range))
                H = NaN .* ones(length(p.adaptOptions))
                R = NaN .* ones(length(p.adaptOptions))
                Construct = NaN .* ones(length(p.adaptOptions))
                StormProt = -9 .* ones(length(p.adaptOptions), length(t_range))
                

                print("\nStartingForLoop\n")
                for i in collect(1:length(p.adaptOptions))
                    # ** Calculate Retreat Costs by Level **                   
                    R[i] = calcHorR(-2, p.adaptOptions[i], lslrPlan_at, p.surgeExposure[m,:], p.adaptOptions)

                    FloodRet[i] = (p.tstep/atstep) * atstep * v.landvalue[m,t]*.04 * calcCoastArea(p.areaparams[m,:], R[i]) + 
                        max(0,calcCoastArea(p.areaparams[m,:], R[i]) - calcCoastArea(p.areaparams[m,:], v.R[m,at_index_prev]))* 
                        (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[m,t] # todo check this summation
            
                    print("\nFloodRet ", FloodRet[i])
                        # check this as well
                    RelocateRet[i] = (p.tstep / atstep) * max(0, calcCoastArea(p.areaparams[m,:], R[i]) - calcCoastArea(p.areaparams[m,:], v.R[m,at_index_prev])) * 
                        p.movefactor * v.ypc_seg[m,t] * 1e-6 * v.popdens[m,t] +
                        p.capmovefactor * p.mobcapfrac * v.capital[m,t] + p.democost * (1 - p.mobcapfrac ) * v.capital[m,t]
                    
                    print("\nRelocate ", RelocateRet[i])
                    for j in t_range
                        StormRet[i, j] = p.tstep * (1 - v.ρ[rgn_ind, j]) * 
                            (p.rσ₀[m] / (1 + p.rσ₁[m] * exp(p.rσ₂[m] * max(0, R[i] - p.lslr[m, j])))) * 
                            (v.capital[m, j] + v.popdens[m, j] * v.vsl[rgn_ind, j] * p.floodmortality) 

                        RetreatTot[i,j] = FloodRet[i] + RelocateRet[i] + StormRet[i,j] + v.WetlandRetreatNoAdapt[m,j]
                    end
                    
                       
                    NPVAdapt[2, i] = sum([v.discountfactor[j] * RetreatTot[i,j] for j in t_range])
                    
                    # ** Calculate protect costs for a subset of adaptation options ** 
                    print("Prot\n")
                    if p.adaptOptions[i] >= 10
                        H[i] = calcHorR(-1, p.adaptOptions[i], lslrPlan_at, p.surgeExposure[m,:], p.adaptOptions)
                        
                        # Flexible mode question: v.H[m, at_prev] for t > 1. 
                        Construct[i] = (p.tstep/atstep) * 
                            sum([p.length[m] * p.pc[m] * (p.pcfixed + (1- p.pcfixed)*(H[i]^2 - v.H[m, at_index_prev]^2) + 
                            p.mc*atstep*H[i]) + p.length[m] * 1.7 * H[i] * v.landvalue[m,j]*.04/2*atstep for j in t_range])

                        for j in t_range
                            StormProt[i,j] = p.tstep * (1 - v.ρ[rgn_ind, j]) * (p.pσ₀[m] + p.pσ₀coef[m] * p.lslr[m, j]) / 
                                                            (1. + p.pσ₁[m] * exp(p.pσ₂[m] * max(0,(H[i] - p.lslr[m,j])))) *
                                                            (v.capital[m,j] + v.popdens[m,j] * v.vsl[rgn_ind, j] * p.floodmortality)
                                    
                            ProtectTot[i,j] = Construct[i] + v.WetlandProtect[m, j] + StormProt[i,j]
            
                        end

                        NPVAdapt[1,i] = sum( [ v.discountfactor[j] * ProtectTot[i,j] for j in t_range] ) # Protect
                        
                    end           
                end

                # ** Choose least cost adaptation option **
                print("\nChooseNPV\n")
                protectInd = indmin(NPVAdapt[1,:])
                retreatInd = indmin(NPVAdapt[2,:])

                minLevels = [p.adaptOptions[protectInd], p.adaptOptions[retreatInd], 0]
                choices = [NPVAdapt[1,protectInd], NPVAdapt[2,retreatInd], NPVNoAdapt]
                leastcost = -1 * indmin(choices)
                leastlevel = minLevels[indmin(choices)]

                v.AdaptationDecision[m, t_range] = leastcost
                v.AdaptationLevel[m, t_range] = leastlevel
                v.AdaptationCost[m, t_range] = ProtectTot[protectInd, t_range]

                # Set H or R
                v.H[m, at_index] = 0
                v.R[m, at_index] = 0

                if leastcost==-1
                    v.H[m, at_index] =  calcHorR(-1, p.adaptOptions[protectInd], lslrPlan_at, p.surgeExposure[m,:], p.adaptOptions)
                elseif leastcost==-2
                    v.R[m, at_index] = calcHorR(-2, p.adaptOptions[retreatInd], lslrPlan_at, p.surgeExposure[m,:], p.adaptOptions)
                end

            end
        end
    end
end



#         if atstep==0
#          # Cost in last period same as cost in previous period if last period starts new adaptation period
#             v.ProtectCost[adapt_range, t] = v.ProtectCost[adapt_range, t-1]
#             v.RetreatCost[adapt_range, t] = v.RetreatCost[adapt_range, t-1]
#             v.NoAdaptCost[t] = v.NoAdaptCost[t-1]
            
#             v.AdaptationCost[t] = v.AdaptationCost[t-1]
#             v.AdaptationDecision[t] = v.AdaptationDecision[t-1]
#             v.AdaptationLevel[t] = v.AdaptationLevel[t-1] 

#         elseif atstep < 0 # FLAG: can take this out for speed/performance reasons if needed
#             error("Adaptation period exceeds total time horizon")               
        
#         else     
#             t_range = collect(t:last_t)
#             print(t_range,"\n")
                        
#             # ---Decision Variables---
#             # Planned for sea level rise (lslrPlan); equivalent to level of slr reached by start of next 
#             #   adaptation period. SLR is assumed to be cumulative in m relative to a t=0 starting position.
#             lslrPlan_at = p.lslr[at_next]

#             # ** Find H and R for construct and retreat cases **
#             # Look up index for '10' adaptation level, as it's treated differently in construct case
#             f(x) = x == 10
#             ten_index = find(f, p.adaptOptions)[1]

#             for i in adapt_range, j in t_range 
#                 v.R[i, j] = max(0, lslrPlan_at + p.surgeExposure[i])
                
#                 if i==ten_index # Construction heights for options below 10 are not calculated
#                     v.H[i,j] = max(0, lslrPlan_at + p.surgeExposure[i] / 2)            
#                 elseif i>ten_index
#                     v.H[i,j] = max(0, lslrPlan_at + p.surgeExposure[i])
#                 end

#             end            

#             # ** Calculate Storm and Wetland Costs for all cases over all t in adaptation period **                    
#             for j in t_range, i in adapt_range                   
                    
#                     v.StormRetreat[i,j] = p.tstep * (1 - v.ρ[j]) * (p.rσ₀ / (1 + p.rσ₁ * exp(p.rσ₂ * max(0, v.R[i,j] - p.lslr[j])))) * 
#                         (v.capital[j] + v.popdens[j] * v.vsl[j] * p.floodmortality)
                    
#                     v.StormNoAdapt[j] = p.tstep * (1 - v.ρ[j]) * (p.rσ₀ / (1 + p.rσ₁ )) * 
#                         (v.capital[j] + v.popdens[j] * v.vsl[j] * p.floodmortality)
                
#                     v.WetlandRetreatNoAdapt[j] = p.tstep * v.wetlandservice[j] * v.wetlandloss[j] * min(v.coastArea[j], p.wetlandarea)
#                     v.WetlandProtect[j] = p.tstep * p.wetlandarea * v.wetlandservice[j]


#                     if i >= ten_index
#                         v.StormProtect[i,j] = p.tstep * (1 - v.ρ[j]) * (p.pσ₀ + p.pσ₀coef * p.lslr[j]) / (1. + p.pσ₁ * exp(p.pσ₂ * max(0,(v.H[i,j] - p.lslr[j])))) *
#                             (v.capital[j] + v.popdens[j] * v.vsl[j] * p.floodmortality)
#                     end

#             end

#             # ** Calculate Total No Adaptation Costs **           
#             if last==1
#                 for i in collect(t:last_t - 1)
#                     # TODO check these functions
#                     # Check time interval - is it correct?
#                     v.FloodNoAdapt[i] = p.tstep * v.landrent[i] * max(0, v.coastArea[i+1]) + (max(0, v.coastArea[i+1]) - max(0, v.coastArea[i])) * 
#                         (1 - p.mobcapfrac) * v.capital[i]                    
        
#                     v.RelocateNoAdapt[i] = (max(0, v.coastArea[i+1]) - max(0,v.coastArea[i])) * (5 * p.movefactor * p.ypc_country[i]*1e-6*v.popdens[i] +
#                         p.capmovefactor * p.mobcapfrac * v.capital[i] + p.democost * (1 - p.mobcapfrac) * v.capital[i])
        
#                     v.NoAdaptCost[i] = v.WetlandRetreatNoAdapt[i] + v.FloodNoAdapt[i] + v.RelocateNoAdapt[i] + v.StormNoAdapt[i]
        
#                 end
#                 v.FloodNoAdapt[last_t] = v.FloodNoAdapt[last_t -1]
#                 v.RelocateNoAdapt[last_t] = v.RelocateNoAdapt[last_t -1]
#                 v.NoAdaptCost[last_t] = v.NoAdaptCost[last_t -1]
#             else
#                 for i in t_range

#                     v.FloodNoAdapt[i] = p.tstep * v.landrent[i] * max(0, v.coastArea[i+1]) + (max(0, v.coastArea[i+1]) - max(0, v.coastArea[i])) * 
#                         (1 - p.mobcapfrac) * v.capital[i]                    
                    
#                     v.RelocateNoAdapt[i] = (max(0, v.coastArea[i+1]) - max(0,v.coastArea[i])) * (5 * p.movefactor * p.ypc_country[i]*1e-6*v.popdens[i] +
#                         p.capmovefactor * p.mobcapfrac * v.capital[i] + p.democost * (1 - p.mobcapfrac) * v.capital[i])
        
#                     v.NoAdaptCost[i] = v.WetlandRetreatNoAdapt[i] + v.FloodNoAdapt[i] + v.RelocateNoAdapt[i] + v.StormNoAdapt[i]
#                 end
#             end

#             print("flag1")
#             # ** Calculate Additional Protect and Retreat Costs **
#             if t==1
#                 for i in adapt_range, j in t_range

#                     if i >= ten_index
#                         v.Construct[i] = (p.tstep/atstep) * sum([p.length * p.pc * (p.pcfixed + (1- p.pcfixed)*(v.H[i,t]^2) + p.mc*atstep*v.H[i,t]) + p.length * 1.7 * v.H[i,t] * v.landrent[j]/2*atstep for j in t_range])                        
#                         v.ProtectCost[i, j] = v.Construct[i] + v.WetlandProtect[j] + v.StormProtect[i,j] 
#                     end 
#                     ### Q: How is time / summation being performed here? Same q for Relocate
#                     v.FloodRetreat[i] = (p.tstep/atstep) * atstep * v.landrent[t] * calcCoastArea(p.areaparams, v.R[i,t]) + 
#                         max(0,calcCoastArea(p.areaparams, v.R[i,t]))* (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[t]

#                     v.RelocateRetreat[i] = (p.tstep / atstep) * max(0, calcCoastArea(p.areaparams, v.R[i,t])) * p.movefactor * v.ypc_seg[t] * 1e-6 * v.popdens[t] +
#                         p.capmovefactor * p.mobcapfrac * v.capital[t] + p.democost * (1 - p.mobcapfrac ) * v.capital[t]

#                     v.RetreatCost[i, j] = v.FloodRetreat[i] + v.RelocateRetreat[i] + v.WetlandRetreatNoAdapt[j] + v.StormRetreat[i,j]

#                 end

#             else
#                 at_prev = Int(p.at[at_index - 1])              
#                 for i in adapt_range, j in t_range

#                     if i >= ten_index
#                         v.Construct[i] = (p.tstep/atstep) * sum([p.length * p.pc * (p.pcfixed + (1- p.pcfixed)*(v.H[i,t]^2 - v.H[i, at_prev]^2) + p.mc*atstep*v.H[i,t]) + p.length * 1.7 * v.H[i,t] * v.landrent[j]/2*atstep for j in t_range])
#                         v.ProtectCost[i, j] = v.Construct[i] + v.WetlandProtect[j] + v.StormProtect[i,j]                        
#                     end

#                     v.FloodRetreat[i] = (p.tstep/atstep) * atstep * v.landrent[t] * calcCoastArea(p.areaparams, v.R[i,t]) + 
#                         max(0,calcCoastArea(p.areaparams, v.R[i,t]) - calcCoastArea(p.areaparams, v.R[i,at_prev]))* (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[t]

#                     v.RelocateRetreat[i] = (p.tstep / atstep) * max(0, calcCoastArea(p.areaparams, v.R[i,t]) - calcCoastArea(p.areaparams, v.R[i,at_prev])) * p.movefactor * v.ypc_seg[t] * 1e-6 * v.popdens[t] +
#                         p.capmovefactor * p.mobcapfrac * v.capital[t] + p.democost * (1 - p.mobcapfrac ) * v.capital[t]

#                     v.RetreatCost[i, j] = v.FloodRetreat[i] + v.RelocateRetreat[i] + v.WetlandRetreatNoAdapt[j] + v.StormRetreat[i,j]
    
#                 end
#             end
#             print("flag2")
#             # ** Calculate Final Adaptation Decision and costs for adaptation period **
#             # TODO: Make sure costs are in correct units          
#             if (t > 1 && p.fixed)
#                 # Fixed version of model: choose adaptation option chosen in previous period
#                 v.AdaptationDecision[t_range] = v.AdaptationDecision[t-1]
#                 v.AdaptationLevel[t_range] = v.AdaptationLevel[t-1]

#                 if v.AdaptationDecision[t] == -1
#                     adaptInd = find(r -> r==v.AdaptationLevel[t], p.adaptOptions)
#                     v.AdaptationCost[t_range] = v.ProtectCost[adaptInd, t_range]

#                 elseif v.AdaptationDecision[t] == -2
#                     adaptInd = find(r -> r==v.AdaptationLevel[t], p.adaptOptions)
#                     v.AdaptationCost[t_range] = v.RetreatCost[adaptInd, t_range]

#                 else
#                     v.AdaptationCost[t_range] = v.NoAdaptCost[t_range]
                    
#                 end

#             else
#                 # Flexible version of model: Make decision based on NPV of adaptation options in period t
#                 print("flexible")
#                 NPVNoAdapt = sum( [ v.discountfactor[j] * v.NoAdaptCost[j] for j in t_range ])
#                 NPVProtect = -9. .* ones(1:length(p.adaptOptions))
#                 NPVRetreat = -9. .* ones(1:length(p.adaptOptions))

#                 for i in adapt_range
#                     if i >= ten_index
#                         NPVProtect[i] = sum( [ v.discountfactor[j] * v.ProtectCost[i,j] for j in t_range] )
#                     end
#                     NPVRetreat[i] = sum( [ v.discountfactor[j] * v.RetreatCost[i,j] for j in t_range])
#                 end
#                 NPVProtect = hcat(fill(-1, length(p.adaptOptions)), p.adaptOptions, NPVProtect)
#                 NPVProtect = NPVProtect[find(z -> z >=0, NPVProtect[:,3]),:] # Remove N/A options (-9)
#                 NPVRetreat = hcat(fill(-2, length(p.adaptOptions)), p.adaptOptions, NPVRetreat)
#                 NPVNoAdapt = hcat(-3, 0., NPVNoAdapt)
    
           
#                 # Store options and indicators in matrices
#                 # NPV Cost Matrix: [adaptationLevel adaptationOption NPV]
#                 NPVCostMatrix = [ NPVProtect;
#                                 NPVRetreat;
#                                 NPVNoAdapt
#                                 ]
    
#                 # # Choose NPV-minimizing adaptation options 
#                 minIndexAll = indmin(NPVCostMatrix[:,3])
#                 minChoiceAll = NPVCostMatrix[minIndexAll, 1]
#                 minLevelAll = NPVCostMatrix[minIndexAll, 2]
    
#                 v.AdaptationDecision[t_range] = minChoiceAll
#                 v.AdaptationLevel[t_range] = minLevelAll
    
#                 # Look up Cost of NPV minimizing adaptation for all t in adaptation period
#                 if minChoiceAll==-1
#                     # Protection Cost
#                     adaptInd = find(p -> p==minLevelAll, p.adaptOptions)
#                     v.AdaptationCost[t_range] = v.ProtectCost[adaptInd, t_range]
       
#                 elseif minChoiceAll==-2
#                     # Retreat Cost
#                     adaptInd = find(r -> r==minLevelAll, p.adaptOptions)
#                     v.AdaptationCost[t_range] = v.RetreatCost[adaptInd, t_range]
         
#                 else
#                     v.AdaptationCost[t_range] = v.NoAdaptCost[t_range] 
      
#                 end
                
#             end
            
#        end



    
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
    # Look up name corresponding to index, or index corresponding to name
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
