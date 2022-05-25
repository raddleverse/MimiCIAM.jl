using DataFrames
using CSVFiles
using Query

## Script to reproduce data processing doen by Lisa Rennels on 05/23/2022

################################################################################
## 1. Construction Cost Index (cci.csv)

# For construction cost indices, also attached the file (read in from cci tab) 
# plus the code used to fill in / limit the range. (source data file is 2011WB_ICP.xls)

# * Construction cost indices
# $CALL gdxxrw.exe construction_country_indices\2011WB_ICP.xls par=cci rng=cci!h8 rdim=1 cdim=0
# $gdxin 2011WB_ICP
# $load cci
# * correct for missing countries by assuming 1 and restrict factor to 0.5 to 2.5
# loop(country,
#         if(cci(country)=0, cci(country)=1; );
#         cci(country)=max(0.5,min(2.5,cci(country)));
# );
 
################################################################################
## 2. GTAP Land Value (gtapland.csv)

# the parameter is read in from the Agland tab cells E4:F191. Note the dollar year 
# adjustment in the code below. (source data file is GTAPagrent.xls)

# * Land values from GTAP
# $onecho > inputlist.txt
# par=gtapland rng=Agland!e4:f191 rdim=1 cdim=0
# par=countryarea rng=area!a2:b166 rdim=1  cdim=0
# $offecho
# $CALL gdxxrw.exe GTAPagrent.xls @inputlist.txt
# $gdxin GTAPagrent
# $load gtapland, countryarea
# * $2007M per sq km - I convert to $2010
# gtapland(country)=gtapland(country)/0.962;

################################################################################
## 3. SSP Data

# TODO

################################################################################
## 4. Reference Population Density (refpopdens.csv)

# The weighted average of reference population density of the segments in this 
# region, weighted by area 1 ie. area between 0 and 1 meter.

xsc = load(joinpath(@__DIR__, "input/xsc.csv")) |> DataFrame
data = load(joinpath(@__DIR__, "input/data.csv")) |> DataFrame

regions = load(joinpath(@__DIR__, "meta/rgnIDmap.csv")) |> DataFrame
df = DataFrame()
for r in regions.rgn
    segIDs = (xsc |> @filter(_.rgn == r) |> DataFrame).segID
    segNames = (xsc |> @filter(_.rgn == r) |> DataFrame).seg

    segData = data |> @filter(_.NA in segNames) |> DataFrame
    select!(segData, ["NA", "area1", "popdens"])
    insertcols!(segData, :weight => (segData.area1) / sum(segData.area1)) 
    value = max(1, sum(segData.popdens .* segData.weight))

    append!(df, DataFrame(:country => r, :refpopdens => value))
end

df |> save(joinpath(@__DIR__, "input/refpopdens.csv"))
