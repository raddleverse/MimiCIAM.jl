### Downscale BRICK from GMSL to LSL
### 3/18/2020

# Retrieve BRICK fingerprints from NetCDF file 
function get_fingerprints()
    fp_file = datadep"BRICK fingerprints/FINGERPRINTS_SLANGEN_Bakker.nc"
    fplat = ncread(fp_file,"lat")
    fplon = ncread(fp_file,"lon")
    fpAIS = ncread(fp_file,"AIS")
    fpGSIC = ncread(fp_file,"GLAC")
    fpGIS = ncread(fp_file,"GIS")
    ncclose()

    return fplat,fplon,fpAIS,fpGSIC,fpGIS
end

# Get brick ensemble members for specified RCP from NetCDF file 
# Returns time x ens arrays for brick components 
function get_brickGMSL(gmslfile,rcp)
    
    brAIS= ncread(gmslfile,"AIS_RCP$(rcp)")
    brGSIC = ncread(gmslfile,"GSIC_RCP$(rcp)")
    brGIS = ncread(gmslfile,"GIS_RCP$(rcp)")
    brTE = ncread(gmslfile,"TE_RCP$(rcp)")
    brLWS = ncread(gmslfile,"LWS_RCP$(rcp)")
    brGMSL = ncread(gmslfile,"GlobalSeaLevel_RCP$(rcp)")
    btime = ncread(gmslfile,"time_proj")
    ncclose()

    return btime, brAIS, brGSIC, brGIS, brTE, brLWS, brGMSL

end

# Get CIAM lonlat tuples for specified segIDs
# segID order does not matter; will sort tuples alphabetically by segment name
function get_lonlat(segIDs)
    ciamlonlat = CSV.read("data/diva_segment_latlon.csv")
    #ciamlonlat = joinpath(@__DIR__,"..","data","diva_segment_latlon.csv")

    if segIDs==false
        filt = DataFrame(ciamlonlat)
    else
        filt = ciamlonlat |> @filter(_.segid in segIDs) |> DataFrame 
    end

    sort!(filt, [:segments])
    lons = filt.longi 
    lats = filt.lati 

    return collect(zip(lons,lats))
end

# Choose n ensemble members randomly within specified percentile range
# if percentile range is one number, will return that percentile (with respect to GMSL in specified year)
# time - BRICK time vector from netcdf 
# ens - BRICK GMSL matrix, time x num ensembles
# low - minimum percentile threshold (integer - e.g. 5 = 5th percentile)
# high - maximum percentile threshold
function choose_ensemble_members(time, ens, n, low, high, yend)
    if length(time)==size(ens)[1]
        end_year = findall(x -> x==yend, time)[1]
        
        val_low = percentile(ens[end_year,:],low)
        val_high = percentile(ens[end_year,:],high)

        ens_inds = findall(x -> x >= val_low && x <= val_high, ens[end_year,:])

        
        chosen_inds = ens_inds[sample(1:end,n,replace=false)]
        return chosen_inds 
      
    else
        println("Error: time dimension mismatch")
    end
end

# downscale_brick: downscale BRICK gmsl to lsl for all segments and ensembles of interest 
# Input: 
# brickcomps - BRICK components (time x ens matrices corresponding to brick gmsl components)
# lonlat - vector of (lon,lat) tuples, sorted corresp to segment name alphabetical order  
# Output: 
# lsl_out: ens x time x segment array of local sea levels, sorted in alphabetical order by segment name
# GMSL: global mean sea levels corresponding to local sea level arrays (time x ens)
function downscale_brick(brickcomps,lonlat, ensInds, ystart=2010, yend=2100, tstep=10)
    # To do - check with vectors of lat, lon 
    (fplat,fplon,fpAIS,fpGSIC,fpGIS) = get_fingerprints()
    (btime,AIS,GSIC,GIS,TE,LWS,GMSL) = brickcomps

    # Select indices of time of interest, with respect to timestep
    tinds = findall( x -> x .>= ystart && x .<=yend, btime)
    years = collect(ystart:yend)
    yinds = findall(x -> x % tstep==0, years)

    tdim=length(btime)
    
    if length(years)==length(tinds)
        tinds = tinds[yinds]
    else
        println("Error: years outside of bounds")
        return nothing
    end

    num_ens = length(ensInds)

    # Output matrix: ens x time x segment 
    lsl_out = zeros(num_ens, length(tinds), length(lonlat))
  
    # Trim component vectors to timesteps and ensembles. Assume interval is 1 year 
    if tdim==size(AIS)[1] # check that time dimension is 1
        AIS = AIS[tinds,ensInds]
        GSIC = GSIC[tinds,ensInds]
        GIS=GIS[tinds,ensInds]
        TE = TE[tinds,ensInds]
        LWS = LWS[tinds,ensInds]
        GMSL = GMSL[tinds,ensInds]
    else
        println("Error: time dimension is not 1 for brick components")
        return nothing
    end

    for f in 1:length(lonlat) # Loop through lonlat tuples 

        lon = lonlat[f][1]
        lat = lonlat[f][2]
        # Convert Longitude to degrees East
        # CIAM Lat is already in (-90,90) by default
        if lon <0
            lon = lon + 360
        end
    
        # Find fingerprint degrees nearest to lat,lon 
        ilat = findall(isequal(minimum(abs.(fplat.-lat))),abs.(fplat.-lat))
        ilon = findall(isequal(minimum(abs.(fplon.-lon))),abs.(fplon.-lon))
        

        # Take average of closest lat/lon values
        fpAIS_flat = collect(skipmissing(Iterators.flatten(fpAIS[ilon,ilat])))
        fpGSIC_flat = collect(skipmissing(Iterators.flatten(fpGSIC[ilon,ilat])))
        fpGIS_flat = collect(skipmissing(Iterators.flatten(fpGSIC[ilon,ilat])))
        
        fpAIS_loc = mean(fpAIS_flat[isnan.(fpAIS_flat).==false],dims=1)[1]
        fpGSIC_loc = mean(fpGSIC_flat[isnan.(fpGSIC_flat).==false],dims=1)[1]
        fpGIS_loc = mean(fpGIS_flat[isnan.(fpGIS_flat).==false],dims=1)[1]
        fpTE_loc = 1.0
        fpLWS_loc=1.0

        # Keep searching nearby lat/lon values if fingerprint value is NaN unless limit is hit 
        inc =1
        
        while isnan(fpAIS_loc) || isnan(fpGIS_loc) || isnan(fpGSIC_loc) && inc<5
       
            newlonStart = lon_subtractor.(fplon[ilon],inc)[1]
            newlatStart = lat_subtractor.(fplat[ilat],inc)[1]
            newlonEnd = lon_adder.(fplon[ilon],inc)[1]
            newlatEnd = lat_adder.(fplat[ilat],inc)[1]
            
            latInd1 = minimum(findall(isequal(minimum(abs.(fplat.-newlatStart))),abs.(fplat.-newlatStart)))
            #minimum(findall(x-> x in newlatStart,fplat))
            latInd2 = maximum(findall(isequal(minimum(abs.(fplat.-newlatEnd))),abs.(fplat.-newlatEnd)))
            #maximum(findall(x -> x in newlatEnd,fplat))

            lonInd1 = minimum(findall(isequal(minimum(abs.(fplon.-newlonStart))),abs.(fplon.-newlonStart)))
            #minimum(findall(x-> x in newlonStart,fplon))
            lonInd2 = maximum(findall(isequal(minimum(abs.(fplon.-newlonEnd))),abs.(fplon.-newlonEnd)))
            #maximum(findall(x -> x in newlonEnd,fplon))
        
            if latInd2 < latInd1
                latInds=[latInd1; 1:latInd2]
            else
                latInds=latInd1:latInd2
            end

            if lonInd2 < lonInd1
                lonInds=[lonInd1; 1:lonInd2]
            else
                lonInds = lonInd1:lonInd2
            end

            fpAIS_flat = collect(skipmissing(Iterators.flatten(fpAIS[lonInds,latInds])))
            fpGSIC_flat = collect(skipmissing(Iterators.flatten(fpGSIC[lonInds,latInds])))
            fpGIS_flat = collect(skipmissing(Iterators.flatten(fpGSIC[lonInds,latInds])))
          
            fpAIS_loc = mean(fpAIS_flat[isnan.(fpAIS_flat).==false],dims=1)[1]
            fpGSIC_loc = mean(fpGSIC_flat[isnan.(fpGSIC_flat).==false],dims=1)[1]
            fpGIS_loc = mean(fpGIS_flat[isnan.(fpGIS_flat).==false],dims=1)[1]

            inc = inc + 1 
        
        end
   
        # If still NaN, throw an error 
        if isnan(fpAIS_loc) || isnan(fpGIS_loc) || isnan(fpGSIC_loc)
            println("Error: no fingerprints found for ($(lon),$(lat))")
            return nothing 
        end

       # Multiply fingerprints by BRICK ensemble members 
        for n in 1:size(AIS)[2] # loop through ensemble members 
            lsl_out[n, :, f] = fpGIS_loc * GIS[:,n] + fpAIS_loc * AIS[:,n] + fpGSIC_loc * GSIC[:,n] + 
                fpTE_loc * TE[:,n] + fpLWS_loc * LWS[:,n]
        end

    end # End lonlat tuple 

    return lsl_out,GMSL 
end

# Driver function to downscale BRICK gmsl for specified segments 
function brick_lsl(rcp,segIDs,brickfile,n,low=5,high=95,ystart=2010,yend=2100,tstep=10)
    brickGMSL = get_brickGMSL(brickfile,rcp)
    brickEnsInds = choose_ensemble_members(brickGMSL[1],brickGMSL[7],n,low,high,yend)
    lonlat = get_lonlat(segIDs)

    (lsl,gmsl) = downscale_brick(brickGMSL, lonlat, brickEnsInds,ystart,yend,tstep)

    return lsl,gmsl
end

function brickCIAM_driver(rcp,brickfile,n,low=5,high=95,ystart=2010,yend=2100,tstep=10)

    # Get segIDs from initfile pre-specified subset
    ciamparams = MimiCIAM.init()
    if ciamparams["subset"][1]==false
        segIDs=false
    else
        subs = MimiCIAM.load_subset(ciamparams["subset"])
        sort!(subs)
        segIDs = MimiCIAM.segStr_to_segID(subs)
    end

    (lsl, gmsl) = brick_lsl(rcp,segIDs,brickfile,n,low,high,ystart,yend,tstep)
    num_ens = n

    if yend==2100
        t=10
    elseif yend==2200
        t=20
    elseif yend <2200 && yend >= 2010
        yend = round(yend/10)*10
        years = collect(ystart:tstep:yend)
        t = length(years)
    else
        println("Error: end year exceeds CIAM bounds.")
        return nothing 
    end

    # Set up CIAM using init.txt defaults for SSP and subsets 
    m = MimiCIAM.get_model(t=t) 
    globalNPV = zeros(n)
    for i in 1:n
        update_param!(m, :lslr, lsl[i,:,:])
        run(m)
        # store output 
        globalNPV[i] = m[:slrcost,:NPVOptimalTotal]
    end 

    return globalNPV
end



function adder(maxval)
    function y(point,n)
        if point + n > maxval
            return point + n - maxval
        else 
            return point + n
        end
    end
end

function subtractor(minval,maxval)
    function y(point,n)
        if point - n < minval
            return min(maxval,point - n + maxval)
        else
            return point - n
        end
    end
end

lon_subtractor = subtractor(1,360)
lon_adder = adder(360)
lat_adder = adder(180)
lat_subtractor = subtractor(1,180)