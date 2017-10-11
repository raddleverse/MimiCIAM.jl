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
# This is a truncated version for testing purposes that only evaluates the coastal
# protection option. 
#-------------------------------------------------------------------------------

using Mimi

@defcomp ciam_protect begin
    # ---Time-related Parameters---
    tstep = Parameter()                     # Length of individual time-step (years)
    at = Parameter( index = [n])            # Array of time indices that mark starts of adaptation periods 
                                            # FLAG_1: This appears to work, but is it the correct approach?
                                            # FLAG_2: Mimi throws an error when you try to pass in non-float values
                                            #   (e.g. bools or ints)
    ntsteps = Parameter()                   # Number of time-steps     

    # ---Socioeconomic Parameters---
    pop_country = Parameter(index = [time]) # Population of country (million people) (from MERGE)
    refpopdens_country = Parameter()        # Reference population density of country (people / km^2)
    refpopdens_usa = Parameter()            # Reference population density of USA (people/km^2) 
    popdens1_seg = Parameter()              # Pop density of segment in time t = 1 (people/km^2)
    popdens = Variable(index = [time])      # Population density extrapolated forward in time (people / km^2)

    ypc_country = Parameter(index = [time]) # GDP per capita per country ($2010 per capita)
    ypc_usa = Parameter(index = [time])     # GDP per capita in USA; used as benchmark ($2010 per capita)
    ypc_seg = Variable(index = [time])      # GDP per capita by segment ($2010 per capita)
    ypc_scale = Variable()                  # Scaling factor to scale country to segment ypc


    # ---Land Parameters---  
    landinput = Parameter()                 # Set to 1.0 for FUND or 0.0 for GTAP
                                            #   FLAG_2 cont'd: Mimi throws error if I try to use Int value                     
    landdata = Variable()                   # Takes on value of either fundland or gtapland
    gtapland = Parameter()                  # GTAP land value in 2007 (million 2010$ / km^2)
    fundland = Variable()                   # FUND land value in 1995 (calculated in run_timestep) (million 2010$ / km^2), 
                                            #   Q maybe import directly? 
    dvbm = Parameter()                      # FUND value of OECD dryland per Darwin et al 1995 converted from $1995 ($2010M per sqkm) (5.376)
    land_appr = Variable(index = [time])    # Land appreciation rate (calculated as regression by Yohe ref Abraham and Hendershott) 
    coastland_scale = Variable(index = [time]) # Scaling factor to scale from interior land value to coast land value
    interior = Variable(index = [time])        # Value of interior land (function of land value and appreciation rate) ($2010M per sqkm)
    coastland = Variable(index = [time])       # Coastal land value (function of interior land and scaling factor) ($2010M per sqkm)
    landvalue = Variable(index = [time])    # Total endowment value of land ($2010M per sqkm)
    landrent = Variable(index = [time])     # Annual rental value of land (equal to landvalue * 0.04)

    ρ = Variable(index = [time])            # Country-wide resilience parameter (logistic function related to GDP)
    capital = Variable(index = [time])      # Total endowment value of capital stock (million $2010 / km^2)
    kgdp = Parameter()                      # Capital output ratio (per MERGE) (3 by default) 
    discountrate = Parameter()              # Discount rate (0.04 by default)
    discountfactor = Variable(index=[time]) # Discount factor (derived from discount rate)

    # ---Coastal Parameters---
    length = Parameter()                    # Segment length (km)

    # ---Protection Parameters---
    pc = Parameter()                        # Cost of protection per segment (million 2010$ / km / vertical m^2 )
    pcfixed = Parameter()                   # Fraction of protection cost that is fixed (not variable in height) (0.3)
    mc = Parameter()                        # Maintenance cost (Hillen et al, 2010) (2%/yr) 
    pc0 = Parameter()                       # Reference cost of protection (million 2010$ / km / vert m^2) (6.02 by default)

    # ---Surge Exposure Parameters---
    # Protection case
    pσ₀ = Parameter()
    pσ₀coef = Parameter()
    pσ₁ = Parameter()                       # psigA in GAMS code
    pσ₂ = Parameter()                       # psigB in GAMS code

    # ---Storm damage parameters---
    floodmortality = Parameter()            # Flood deaths as percent of exposed population; (Jonkman Vrijling 2008) (0.01) 
    vsl = Variable(index = [time])          # Value of statistica life (million 2010$)
                                           
    # ---Wetland Loss Parameters---
    wbvm = Parameter()                          # Annual value of wetland services (million 2010$ / km^2 / yr); (Brander et al 2006)  (0.376) 
    wetlandservice = Variable(index = [time])   # Annual value of wetland services TODO change doc
    wetlandarea = Parameter()                   # Initial wetland area in coastal segment (km^2)

    # ---Sea Level Rise Parameters---
    lslr = Parameter(index = [time])         # Local sea level rise (m) 
    lslrPlan = Parameter(index = [time])     # Planned-for sea level rise (m)
                                             # FLAG_3: Still confused about lslr vs lslrPlan                   
    slr10 = Parameter()                                           
    slr100 = Parameter()
    slr1000 = Parameter()
    slr10000 = Parameter()
 
 
    # ---Intermediate Variables---
    # Protection Costs 
    H10 = Variable(index= [time])
    H100 = Variable(index = [time])
    H1000 = Variable(index = [time])
    H10000 = Variable(index = [time])

    StormProtect10 = Variable(index = [time])       # Storm protection costs by wall height
    StormProtect100 = Variable(index = [time])
    StormProtect1000 = Variable(index = [time])
    StormProtect10000 = Variable(index = [time])
    WetlandProtect = Variable(index = [time])       # Value of wetland loss under protection 
    Protect10 = Variable(index = [time])
    Protect100 = Variable(index = [time])
    Protect1000 = Variable(index = [time])
    Protect10000 = Variable(index = [time])
    
    # ---Decision Variables---                                                    
    ProtectCost = Variable(index = [time])          # Total cost of protection for adaptation period (2010$) Q million?
    ProtectLevel = Variable(index = [time])         # Vector of level of protection chosen   

    # ---Outcome Variables---
    AdaptationOption = Variable(index = [time])     # Option chosen for adaptation period
    AdaptationCost = Variable(index = [time])       # Cost of option chosen for adaptation period (discretized according to tstep) 
    AdaptationLevel = Variable(index = [time])      # Level of protect or retreat (if chosen)


end

function run_timestep(s::ciam_protect, t::Int)
    p = s.Parameters
    v = s.Variables
        
    # Calculate intermediate variables if needed (t = 1, pre-adaptation and at start of each adaptation period)
    if t==1
        # Determine land input value (1.0 = FUND, 0.0 = GTAP)
        if p.landinput== 1.
            v.fundland = min(p.dvbm, max(0.005, p.dvbm * p.ypc_country[1] * p.refpopdens_country / (p.ypc_usa[1] * p.refpopdens_country)))
            v.landdata = v.fundland
        else
            v.landdata = p.gtapland
        end

        v.land_appr[t] = 1.
        v.wetlandservice[t] = p.wbvm * ((p.ypc_country[t] / p.ypc_usa[1])^1.16 * (p.refpopdens_country /27.59)^0.47) 
        v.popdens[t] = p.popdens1_seg 
        v.ypc_scale = max(0.9, (p.popdens1_seg/250)^0.05)
        v.ypc_seg[t] = p.ypc_country[t] * v.ypc_scale
        v.capital[t] = p.kgdp * v.ypc_seg[t] * v.popdens[t] * 1e-6
        v.ρ[t] = p.ypc_country[t] / (p.ypc_country[t] + p.ypc_usa[1])
        v.vsl[t] = 1e-6 * 216 * p.ypc_usa[t] * (p.ypc_country[t]/p.ypc_usa[t])^0.05
        v.interior[t] = v.land_appr[t] * v.landdata  
        v.coastland_scale[t] = max(0.5, log(1+v.popdens[t])/log(25))
        v.coastland[t] = v.coastland_scale[t] * v.interior[t]
        v.landvalue[t] = min(v.coastland[t], v.interior[t])
        v.landrent[t] = v.landvalue[t] * 0.04
        v.discountfactor[t] = 1/(1 + p.discountrate)^(p.tstep * (t-1))


    elseif ((t < p.at[1]) | (t in p.at))

        v.land_appr[t] = v.land_appr[t-1] * exp(0.565 * growthrate(p.ypc_country[t-1], p.ypc_country[t]) + 0.313 * growthrate(p.pop_country[t-1], p.pop_country[t]))
        v.wetlandservice[t] = v.land_appr[t] * v.wetlandservice[1]
        v.popdens[t] = v.popdens[t-1] * (1 + growthrate(p.pop_country[t-1], p.pop_country[t])) 
        v.ypc_seg[t] = p.ypc_country[t] * v.ypc_scale
        v.capital[t] = p.kgdp * v.ypc_seg[t] * v.popdens[t] * 1e-6        
        v.ρ[t] = p.ypc_country[t] / (p.ypc_country[t] + p.ypc_usa[1]) 
        v.vsl[t] = 1e-6 * 216 * p.ypc_usa[t] * (p.ypc_country[t]/p.ypc_usa[t])^0.05        
        v.interior[t] = v.land_appr[t] * v.landdata 
        v.coastland_scale[t] = max(0.5, log(1+v.popdens[t])/log(25))
        v.coastland[t] = v.coastland_scale[t] * v.interior[t]
        v.landvalue[t] = min(v.coastland[t], v.interior[t])
        v.landrent[t] = v.landvalue[t] * 0.04
        v.discountfactor[t] = 1/(1 + p.discountrate)^(p.tstep * (t-1))

        
    end

    # --- Calculate adaptation decisions ---
    # Case 1: Pre-adaptation period
    if t < p.at[1]
        v.AdaptationCost[t] = 0
        v.AdaptationOption[t] = -9  # Using -9 for "N/A" code for now
        v.AdaptationLevel[t] = -9 
               
    # Case 2: In the middle of adaptation period
    #   Cost already computed for time t, so do nothing
    elseif (t > p.at[1] && !(t in p.at))           

    # # Case 3: Start of new (or first) adaptation period
    else
        # Determine length of adaptation period ("atstep")
        f(x) = x == t
        at_index = find(f, p.at)[1]
        at_next = at_index + 1

        if at_next <= length(p.at)
            atstep = (p.at[at_next] - p.at[at_index])*p.tstep   # years
            next = Int(p.at[at_next])                           # index  
        else
            # Deal with special case of last adaptation period
            atstep = p.tstep*p.ntsteps - (p.at[at_index] * p.tstep) 
            next = p.ntsteps # Flag this assumes timesteps are indices not years (1:20 not 2010:2100)
        end

        if atstep==0
         # Cost in last period same as cost in previous period if last period starts new adaptation period
            v.ProtectLevel[t] = v.ProtectLevel[t-1]
            v.ProtectCost[t] = v.ProtectCost[t-1]
            v.AdaptationCost[t] = v.AdaptationCost[t-1]
            v.AdaptationOption[t] = v.AdaptationOption[t-1]
            v.AdaptationLevel[t] = v.AdaptationLevel[t-1] 

        elseif atstep < 0
            error("Adaptation period exceeds total time horizon")               
        
        else     
            # Begin optimization calculation 
            t_range = collect(t:(next-1))
            print(t_range)

            # Calculate intermediate values for all future time periods in atstep
            #   First time period was already calculated
            for i in collect((t+1):(next-1))
                v.land_appr[i] = v.land_appr[i-1] * exp(0.565 * growthrate(p.ypc_country[i-1], p.ypc_country[i]) + 0.313 * growthrate(p.pop_country[i-1], p.pop_country[i]))
                v.wetlandservice[i] = v.land_appr[i] * v.wetlandservice[1]
                v.popdens[i] = v.popdens[i-1] * (1 + growthrate(p.pop_country[i-1], p.pop_country[i]))
                v.ypc_seg[i] = p.ypc_country[i] * v.ypc_scale            
                v.capital[i] = p.kgdp * v.ypc_seg[i] * v.popdens[i] * 1e-6            
                v.ρ[i] = p.ypc_country[i] / (p.ypc_country[i] + p.ypc_usa[1])
                v.vsl[i] = 1e-6 * 216 * p.ypc_usa[i] * (p.ypc_country[i]/p.ypc_usa[i])^0.05            
                v.interior[i] = v.land_appr[i] * v.landdata 
                v.coastland_scale[i] = max(0.5, log(1+v.popdens[i])/log(25))
                v.coastland[i] = v.coastland_scale[i] * v.interior[i]
                v.landvalue[i] = min(v.coastland[i], v.interior[i])
                v.landrent[i] = v.landvalue[i] * 0.04
                v.discountfactor[i] = 1/(1 + p.discountrate)^(p.tstep * (i-1))
            end
        
            # ---Decision Variables---
            # Planned for sea level rise (slrPlan)
            # FLAG_3: do not know how this works in GAMS code; just guessing right now
            lslrPlan_at = sum([ p.lslrPlan[i] for i in t_range ])

        
            # Get previous construction height and retreat radius
            #   FLAG_4: Code appears to use H10_prev to account for previous adaptation, 
            #       which would mean no stock effects if they chose a different option in the previous period. 
            #       That seems wrong to me and I'd like to change it eventually but for now sticking with code. 

            # Calculate Construction Heights and Costs  
            # Only doing this for adaptation period, not individual timesteps       
            v.H10[t] = max(0, lslrPlan_at + p.slr10 / 2)
            v.H100[t] = max(0, lslrPlan_at + p.slr100)
            v.H1000[t] = max(0, lslrPlan_at + p.slr1000)
            v.H10000[t] = max(0, lslrPlan_at + p.slr10000)

            if at_index > 1
                at_prev = Int(p.at[at_index - 1])
                H10_prev = v.H10[at_prev] 
                H100_prev = v.H100[at_prev]
                H1000_prev = v.H1000[at_prev]
                H10000_prev = v.H10000[at_prev]
            
            else
                H10_prev = 0
                H100_prev = 0
                H1000_prev = 0
                H10000_prev = 0
            end

            # Calculate construction cost for timestep t in adaptation period (uniform for all timesteps) 
            # FLAG_5: Math: (H[t]^2 - H[t-1]^2) or (H[t]- H[t-1])^2 ? Going with GAMs verison for now
            Construct10 = (p.tstep / atstep) * sum([p.length * p.pc * (p.pcfixed + (1- p.pcfixed)*(v.H10[Int(p.at[at_index])]^2 - H10_prev^2) + p.mc*atstep*v.H10[Int(p.at[at_index])]) + p.length * 1.7 * v.H10[Int(p.at[at_index])] * v.landrent[i]/2*atstep for i in t_range])
            Construct100 = (p.tstep / atstep) * sum([p.length * p.pc * (p.pcfixed + (1- p.pcfixed)*(v.H100[Int(p.at[at_index])]^2 - H100_prev^2) + p.mc*atstep*v.H100[Int(p.at[at_index])]) + p.length * 1.7 * v.H100[Int(p.at[at_index])] * v.landrent[i]/2*atstep for i in t_range])
            Construct1000 = (p.tstep / atstep) * sum([p.length * p.pc * (p.pcfixed + (1- p.pcfixed)*(v.H1000[Int(p.at[at_index])]^2 - H1000_prev^2) + p.mc*atstep*v.H1000[Int(p.at[at_index])]) + p.length * 1.7 * v.H1000[Int(p.at[at_index])] * v.landrent[i]/2*atstep for i in t_range])
            Construct10000 = (p.tstep / atstep) * sum([p.length * p.pc * (p.pcfixed + (1- p.pcfixed)*(v.H10000[Int(p.at[at_index])]^2 - H10000_prev^2) + p.mc*atstep*v.H10000[Int(p.at[at_index])]) + p.length * 1.7 * v.H10000[Int(p.at[at_index])] * v.landrent[i]/2*atstep for i in t_range])
            
            # Calculate Storm Costs by Height and Wetland Costs for Protect Case (loop over adaptation period)
            # Also calculate total protection costs by period t 
                
            for i in t_range
                v.StormProtect10[i] = p.tstep * (1 - v.ρ[i]) * (p.pσ₀ + p.pσ₀coef * p.lslr[i]) / (1. + p.pσ₁ * exp(p.pσ₂ * max(0,(v.H10[Int(p.at[at_index])] - p.lslr[i])))) *
                    (v.capital[i] + v.popdens[i] * v.vsl[i] * p.floodmortality)
            
                v.StormProtect100[i] = p.tstep * (1 - v.ρ[i]) * (p.pσ₀ + p.pσ₀coef * p.lslr[i]) / (1. + p.pσ₁ * exp(p.pσ₂ * max(0,(v.H100[Int(p.at[at_index])] - p.lslr[i])))) *
                    (v.capital[i] + v.popdens[i] * v.vsl[i] * p.floodmortality)
          
                v.StormProtect1000[i] = p.tstep * (1 - v.ρ[i]) * (p.pσ₀ + p.pσ₀coef * p.lslr[i]) / (1. + p.pσ₁ * exp(p.pσ₂ * max(0,(v.H1000[Int(p.at[at_index])] - p.lslr[i])))) *
                    (v.capital[i] + v.popdens[i] * v.vsl[i] * p.floodmortality)
            
                v.StormProtect10000[i] = p.tstep * (1 - v.ρ[i]) * (p.pσ₀ + p.pσ₀coef * p.lslr[i]) / (1. + p.pσ₁ * exp(p.pσ₂ * max(0,(v.H10000[Int(p.at[at_index])] - p.lslr[i])))) *
                    (v.capital[i] + v.popdens[i] * v.vsl[i] * p.floodmortality)

                v.WetlandProtect[i] = p.tstep * p.wetlandarea * v.wetlandservice[i]
            
                v.Protect10[i] = Construct10 + v.StormProtect10[i] + v.WetlandProtect[i]
                v.Protect100[i] = Construct100 + v.StormProtect100[i] + v.WetlandProtect[i]
                v.Protect1000[i] = Construct1000 + v.StormProtect1000[i] + v.WetlandProtect[i]
                v.Protect10000[i] = Construct10000 + v.StormProtect10000[i] + v.WetlandProtect[i]
            
            end
    
            #  Calculate NPV of total protection cost by adaptation level
            NPVProtect10 = sum( [ v.discountfactor[i] * v.Protect10[i] for i in t_range] )
            NPVProtect100 = sum( [ v.discountfactor[i] * v.Protect100[i] for i in t_range] )
            NPVProtect1000 = sum( [ v.discountfactor[i] * v.Protect1000[i] for i in t_range] )
            NPVProtect10000 = sum( [ v.discountfactor[i] * v.Protect10000[i] for i in t_range] )
        
            # Store options and indicators in matrices
            # NPV Cost Matrix: [adaptationLevel adaptationOption NPV]
            # FLAG_6: I have no idea if this is the most computationally efficient way to do this; appreciate thoughts
            # Cannot assign strings to variables; using -1 to corespond to protect case, 
            #   -2 to Retreat, -3 to noadapt cases  
            NPVCostMatrix = [10 -1 NPVProtect10;
                      100 -1 NPVProtect100;
                      1000 -1 NPVProtect1000; 
                      10000 -1 NPVProtect10000
                      ]

            # Costs in individual years also stored in matrices indexed as (option, t)
            #   Row number corresponds to row in NPV cost matrix
            #   Column number corresponds to time period t
            TimestepCostMatrix = [ transpose(v.Protect10);
                            transpose(v.Protect100);
                            transpose(v.Protect1000);
                            transpose(v.Protect10000) 
                            ]

            # Choose cost-minimizing protection option 
            NPVProtectMatrix = NPVCostMatrix[find(x -> x == -1, NPVCostMatrix[:,2]),:] # Subset matrices to Protect Cases
            TimestepProtectMatrix = TimestepCostMatrix[find(x -> x == -1, NPVCostMatrix[:,2]),:]
            minProtectIndex = indmin(NPVProtectMatrix[:,3])
            v.ProtectLevel[t] = NPVProtectMatrix[minProtectIndex,1]
                             
            # Find overall cost-minimizing option
            minCostIndex = indmin(NPVCostMatrix[:,3])

            v.AdaptationOption[t] = NPVCostMatrix[minCostIndex, 2]
            v.AdaptationLevel[t] = NPVCostMatrix[minCostIndex, 1]

            # Look up costs in each timestep for option chosen from matrices 
            #   Does vary somewhat by t because of differing storm, wetland costs within at
            for i in t_range
                v.ProtectCost[i] = TimestepProtectMatrix[minProtectIndex, i]
                v.ProtectLevel[i] = v.ProtectLevel[t]
                v.AdaptationCost[i] = TimestepCostMatrix[minCostIndex, i]
                v.AdaptationOption[i] = v.AdaptationOption[t]
                v.AdaptationLevel[i] = v.AdaptationLevel[t]
            end
            
        end
    end
end

    
# Helper functions
function growthrate(x1, x2)
    epsilon = 1.e-9
    return (x2 / (x1 + epsilon) - 1.)
end
