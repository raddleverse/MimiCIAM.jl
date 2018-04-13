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
    pop = Parameter(index = [regions, time])           # Population of region (million people) (from MERGE)
    refpopdens = Parameter( index = [regions])         # Reference population density of region (people / km^2)
    refpopdens_usa = Parameter()                       # Reference population density of USA (people/km^2) 
    popdens = Parameter( index = [segments])           # Pop density of segment in time t = 1 (people/km^2)
    ypcc = Parameter(index = [regions, time])          # GDP per capita per region ($2010 per capita)
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
    WetlandNoAdapt = Variable(index = [segments, time])
    FloodNoAdapt = Variable(index = [segments, time])
    StormNoAdapt = Variable(index = [segments, time])
    RelocateNoAdapt = Variable(index = [segments, time])

    Construct = Variable(index = [segments, time, 4])
    WetlandProtect = Variable(index = [segments, time])
    StormProtect = Variable(index = [segments, time])
    
    WetlandRetreat = Variable(index = [segments, time])
    StormRetreat = Variable(index = [segments, time])
    FloodRetreat = Variable(index = [segments, time, 5])
    RelocateRetreat = Variable(index = [segments, time, 5])

    # ---Decision Variables---   
    NoAdaptCost = Variable(index = [segments, time])         # Cost of not adapting (e.g. reactive retreat) (2010$)
    ProtectCost = Variable(index = [segments, time, 4])      # Total cost of protection at each level      
    RetreatCost = Variable(index = [segments, time, 5])      # Total cost of retreat at each level   
    OptimalFixedCost = Variable(index = [segments, time])        # Fixed optimal cost based on NPV in period 1   
    OptimalFixedLevel = Variable(index = [segments])        # Fixed optimal level (1,10,100,1000,10000)
    OptimalFixedOption = Variable(index = [segments])       # Fixed adaptation decision (-1 - protect, -2 - retreat, -3 - no adapt) 
    NPVRetreat = Variable(index = [segments, adaptPers, 5])
    NPVProtect = Variable(index = [segments, adaptPers, 4])
    NPVNoAdapt = Variable(index = [segments, adaptPers])                     
  
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
                v.fundland[r] = min(p.dvbm, max(0.005, p.dvbm * p.ypcc[r,1] * p.refpopdens[r] / (p.ypc_usa[1] * p.refpopdens_usa)))
                v.landdata[r] = v.fundland[r]
            else
                v.landdata[r] = p.gtapland[r]
            end

            v.land_appr[r, t] = 1.
            v.wetlandservice[r, t] = p.wbvm * ((p.ypcc[r,t] / p.ypc_usa[1])^1.16 * (p.refpopdens[r] /27.59)^0.47) 
            v.ρ[r, t] = p.ypcc[r,t] / (p.ypcc[r,t] + p.ypc_usa[1])
            v.vsl[r, t] = 1e-6 * 216 * p.ypc_usa[t] * (p.ypcc[r,t]/p.ypc_usa[t])^0.5

            
            for i in collect(2:Int(p.ntsteps))
                v.land_appr[r, i] = v.land_appr[r, i-1] * exp(0.565 * growthrate(p.ypcc[r,i-1], p.ypcc[r,i]) + 0.313 * growthrate(p.pop[r,i-1], p.pop[r,i]))
                v.wetlandservice[r,i] = v.land_appr[r,i] * v.wetlandservice[r,1]
                v.ρ[r, i] = p.ypcc[r,i] / (p.ypcc[r,i] + p.ypc_usa[1]) 
                v.vsl[r, i] = 1e-6 * 216 * p.ypc_usa[i] * (p.ypcc[r,i]/p.ypc_usa[i])^0.5  
            end    
        end

        # 3. Initialize segment-dependent variables 
        for m in d.segments
            rgn_ind = getregion(m, p.xsc)
 
            v.popdens_seg[m, t] = p.popdens[m]
            v.ypc_seg[m, t] = p.ypcc[rgn_ind,t] * max(0.9, (p.popdens[m]/250.)^0.05)
            v.capital[m, t] = p.kgdp * v.ypc_seg[m, t] * v.popdens_seg[m, t] * 1e-6
            v.coastland[m, t] = max(0.5, log(1+v.popdens_seg[m, t])/log(25)) * (v.land_appr[rgn_ind, t] * v.landdata[rgn_ind])  # Interior * scaling factor
            v.landvalue[m, t] = min(v.coastland[m, t], (v.land_appr[rgn_ind, t] * v.landdata[rgn_ind]))
            v.coastArea[m, t] = calcCoastArea(v.areaparams[m,:], p.lslr[m, t])  

            
            for i in collect(2:Int(p.ntsteps))
                v.popdens_seg[m,i] = v.popdens_seg[m, i-1] * (1 + growthrate(p.pop[rgn_ind, i-1], p.pop[rgn_ind,i])) 
                v.ypc_seg[m, i] = p.ypcc[rgn_ind,i] * max(0.9, (p.popdens[m]/250.)^0.05) # ypcc * popdens scaling factor
                v.capital[m, i] = p.kgdp * v.ypc_seg[m, i] * v.popdens_seg[m, i] * 1e-6 
                v.coastland[m, i] = max(0.5, log(1+v.popdens_seg[m, i])/log(25)) * (v.land_appr[rgn_ind, i] * v.landdata[rgn_ind])
                v.landvalue[m, i] = min(v.coastland[m, i], (v.land_appr[rgn_ind, i] * v.landdata[rgn_ind]))
                v.coastArea[m, i] = calcCoastArea(v.areaparams[m,:], p.lslr[m, i])
                v.wetlandloss[m, i-1] = min(1, (localrate(p.lslr[m, i-1], p.lslr[m, i], p.tstep)/p.wmaxrate)^2)
            end

            v.wetlandloss[m, p.ntsteps] = min(1, (localrate(p.lslr[m, p.ntsteps-1], p.lslr[m, p.ntsteps], p.tstep)/p.wmaxrate)^2)  

        end

    end

    if (t in p.at)           
        adapt_range = collect(1:length(p.adaptOptions))         

        # Determine length of adaptation period ("atstep")
        g(c) = c == t
        at_index = find(g, p.at)[1] # WEIRD - changing to g solved issue
        at_index_next = at_index + 1
        at_index_prev = max(1,at_index - 1)
        
        at_prev = Int(p.at[at_index_prev])      # TODO pick one - either index into vector or use the period, it's confusing

        if at_index_next <= length(p.at)
            atstep = (p.at[at_index_next] - p.at[at_index])*p.tstep   # years
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
        t_range = collect(t:last_t)

        for m in d.segments
            if atstep==0
                v.AdaptationCost[m,t] = v.AdaptationCost[m, t-1]
                v.AdaptationDecision[m,t] = v.AdaptationDecision[m, t-1]
                v.AdaptationLevel[m, t] = v.AdaptationLevel[m, t-1]
            else
                rgn_ind = getregion(m, p.xsc)
                
                # ** Calculate No Adaptation Costs **
                    for i in t_range
                        v.StormNoAdapt[m, i] = p.tstep * (1 - v.ρ[rgn_ind , i]) * (p.rsig0[m] / (1 + p.rsigA[m] )) * 
                                (v.capital[m, i] + v.popdens_seg[m, i] * v.vsl[rgn_ind, i] * p.floodmortality)
                            
                        v.WetlandNoAdapt[m,i] = p.tstep * v.wetlandservice[rgn_ind, i] * v.wetlandloss[m, i] * 
                                min(v.coastArea[m, i], p.wetland[m])

                        if i==p.ntsteps
                            v.FloodNoAdapt[m,i] = p.tstep * v.landvalue[m,i-1]*.04 * max(0, v.coastArea[m, i]) + (max(0, v.coastArea[m, i]) - max(0, v.coastArea[m, i-1])) * 
                                (1 - p.mobcapfrac) * v.capital[m, i-1]
        
                            v.RelocateNoAdapt[m,i] = (max(0, v.coastArea[m, i]) - max(0,v.coastArea[m, i-1])) * (5 * p.movefactor * v.ypc_seg[m,i-1]*1e-6*v.popdens_seg[m, i-1] +
                                p.capmovefactor * p.mobcapfrac * v.capital[m, i-1] + p.democost * (1 - p.mobcapfrac) * v.capital[m, i-1])
                        else
                            v.FloodNoAdapt[m,i]  = p.tstep * v.landvalue[m,i]*.04 * max(0, v.coastArea[m, i+1]) + (max(0, v.coastArea[m, i+1]) - max(0, v.coastArea[m,i])) * 
                                (1 - p.mobcapfrac) * v.capital[m,i]                    
                        
                            v.RelocateNoAdapt[m,i] = (max(0, v.coastArea[m,i+1]) - max(0,v.coastArea[m,i])) * (5 * p.movefactor * v.ypc_seg[m,i]*1e-6*v.popdens_seg[m,i] +
                                p.capmovefactor * p.mobcapfrac * v.capital[m,i] + p.democost * (1 - p.mobcapfrac) * v.capital[m,i])
                        end
                            
                        v.NoAdaptCost[m,i] = v.WetlandNoAdapt[m,i] + v.FloodNoAdapt[m,i] +  v.RelocateNoAdapt[m,i] + v.StormNoAdapt[m, i]
                      
                    end
                    v.NPVNoAdapt[m, at_index] = sum( [ v.discountfactor[j] * v.NoAdaptCost[m,j] for j in t_range] )

                    # ** Calculate Protectio and Retreat Costs for Each Adaptation Option **
                    lslrPlan_at = p.lslr[m, at_next]
                    lslrPlan_atprev = p.lslr[m, t]

                    for i in 1:length(p.adaptOptions)

                        R = calcHorR(-2, p.adaptOptions[i], lslrPlan_at, v.surgeExposure[m,:], p.adaptOptions)
                        H = calcHorR(-1, p.adaptOptions[i], lslrPlan_at, v.surgeExposure[m,:], p.adaptOptions)

                        if t==1
                            Rprev = calcHorR(-2, p.adaptOptions[i], p.lslr[m,1], v.surgeExposure[m,:], p.adaptOptions) 
                            Hprev = calcHorR(-1, p.adaptOptions[i], p.lslr[m,1], v.surgeExposure[m,:], p.adaptOptions)
                        else
                            Rprev = calcHorR(-2, p.adaptOptions[i], lslrPlan_atprev, v.surgeExposure[m,:], p.adaptOptions)
                            Hprev = calcHorR(-1, p.adaptOptions[i], lslrPlan_atprev, v.surgeExposure[m,:], p.adaptOptions)
                        end


                        v.FloodRetreat[m,at_index,i] = (p.tstep/atstep) * (atstep * v.landvalue[m,t]*.04 * calcCoastArea(v.areaparams[m,:], R) +          
                            max(0,calcCoastArea(v.areaparams[m,:], R) - calcCoastArea(v.areaparams[m,:], Rprev))* 
                            (1 - p.depr) * (1 - p.mobcapfrac) * v.capital[m,t])

                        v.RelocateRetreat[m,at_index,i] = (p.tstep / atstep) * 
                            max(0, calcCoastArea(v.areaparams[m,:], R) - calcCoastArea(v.areaparams[m,:], Rprev)) * 
                            (p.movefactor * v.ypc_seg[m,t] * 1e-6 * v.popdens_seg[m,t] +
                            p.capmovefactor * p.mobcapfrac * v.capital[m,t] + p.democost * (1 - p.mobcapfrac ) * v.capital[m,t])
       

                        if p.adaptOptions[i] >= 10
                            v.Construct[m,at_index,i-1] = (p.tstep/atstep) * 
                                (p.length[m] * p.pc0 * p.cci[rgn_ind] * (p.pcfixed + (1- p.pcfixed)*(H^2 - Hprev^2) + 
                                p.mc*atstep*H) + p.length[m] * 1.7 * H * v.landvalue[m,t]*.04/2*atstep)
                                
                        end

                        for j in t_range
                            v.WetlandRetreat[m,j] = p.tstep * v.wetlandservice[rgn_ind, j] * v.wetlandloss[m, j] * 
                                            min(v.coastArea[m, j], p.wetland[m])

                            v.StormRetreat[m,j] = p.tstep * (1 - v.ρ[rgn_ind, j]) * 
                                    (p.rsig0[m] / (1 + p.rsigA[m] * exp(p.rsigB[m] * max(0, R - p.lslr[m, j])))) * 
                                    (v.capital[m, j] + v.popdens_seg[m, j] * v.vsl[rgn_ind, j] * p.floodmortality)

                            v.FloodRetreat[m, j, i] = v.FloodRetreat[m, at_index, i]
                            v.RelocateRetreat[m,j,i] = v.RelocateRetreat[m,at_index,i]
                                    
                            v.RetreatCost[m, j, i] = v.FloodRetreat[m,j,i] + v.RelocateRetreat[m,j,i] + v.StormRetreat[m,j] + v.WetlandRetreat[m,j]
                            
                            if p.adaptOptions[i] >= 10
                                v.WetlandProtect[m,j] = p.tstep * p.wetland[m] .* v.wetlandservice[rgn_ind, j]
                                        
                                v.StormProtect[m,j] = p.tstep * (1 - v.ρ[rgn_ind, j]) * (p.psig0[m] + p.psig0coef[m] * p.lslr[m, j]) / 
                                                        (1. + p.psigA[m] * exp(p.psigB[m] * max(0,(H - p.lslr[m,j])))) *
                                                        (v.capital[m,j] + v.popdens_seg[m,j] * v.vsl[rgn_ind, j] * p.floodmortality)
                                
                                v.Construct[m,j,i-1] = v.Construct[m, at_index, i-1]
                                                
                                v.ProtectCost[m,j,i-1] = v.Construct[m,j,i-1] + v.WetlandProtect[m,j] + v.StormProtect[m,j]

                             end

                        end

                        v.NPVRetreat[m, at_index,i] = sum([v.discountfactor[j] * v.RetreatCost[m,findind(j,t_range),i] for j in t_range])

                        if p.adaptOptions[i] >=10
                            v.NPVProtect[m,at_index,i-1] = sum( [ v.discountfactor[j] * v.ProtectCost[m,findind(j,t_range),i-1] for j in t_range] ) # Protect
                        end
                     end

                    # ** Choose Least Cost Option **
                    if t==1 && p.fixed
                        protectInd = indmin(v.NPVProtect[m,at_index,:])
                        retreatInd = indmin(v.NPVRetreat[m,at_index,:])
        
                        minLevels = [p.adaptOptions[protectInd+1], p.adaptOptions[retreatInd], 0]
                        choices = [v.NPVProtect[m,at_index,protectInd], v.NPVRetreat[m,at_index,retreatInd], v.NPVNoAdapt[m,at_index]]
                        leastcost = -1 * indmin(choices)
                        leastlevel = minLevels[indmin(choices)]
        
                        v.OptimalFixedOption[m] = leastcost
                        v.OptimalFixedLevel[m] = leastlevel
                        
                    end
                    
                    if v.OptimalFixedOption[m]==-1
                        v.OptimalFixedCost[m, t_range] = v.ProtectCost[m, t_range, (find(i->i==v.OptimalFixedLevel[m], p.adaptOptions)-1)] 
                    elseif v.OptimalFixedOption[m]==-2
                        v.OptimalFixedCost[m, t_range] = v.RetreatCost[m, t_range, find(i->i==v.OptimalFixedLevel[m], p.adaptOptions)] 
                    else
                        v.OptimalFixedCost[m, t_range] = v.NoAdaptCost[m, t_range]
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

function pos(x)
    return max(0,x)
end