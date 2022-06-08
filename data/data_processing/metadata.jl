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
## 3. Population and GDP (pop.csv and ypcc.csv)

# The initial CIAM model relied on demographic inputs from EPRI’s global IAM called 
# MERGE (Model for Estimating the Regional and Global Effects of Greenhouse Gas 
# Reductions). Some minor adjustments for use in CIAM (i.e., aligning time indices, 
# extrapolating GDP growth after 2100 at a nominal rate, and converting to $2010)
# were made and exported into a GDX and then csv file MERGEdata.xls

# This Excel spreadsheet was reconfigured to create input files pop.csv and ypcc.csv

# NOTE that the MERGE model does not have data for PSE. For replication purposes
# we substitute in the SSP2 data through 2100, and repeat after that to match
# the process carried out with the MERGE data, sourced from models listed below 
# in section 4. This is only consequential for replication purposes.

# These two data files are ONLY used in replication exercises of the CIAM 
# (Diaz et al., 2016).

################################################################################
## 4. SSP Data (pop_IIASAGDP_SSPX_v9_130219 and ypcc_IIASAGDP_SSPX_v9_130219)

# The SSP Data were downloaded from to the file SspDb_country_data_2013-06-12.csv, 
# which was then filtered for two variables: Population and GDP|PPP. 

# Starting with pop.csv and ypcc.csv, any country available in the IIASA GDP model
# is replaced with data from the respective SSP/variable combination.

# An exception is made for PSE's GDP and GDP per Capita. Since this is unavailable 
# in MERGE, and GDP is not availabe in IIASA GDP, we use population from IIASA GDP
# but Per Capita GDP from OECD Env-Growth, the only model for which GDP is availabe.
# Note we calcualte Per Capita GDP from OECD Env-Growth's Population and GDP to be
# internally consistent, although upon inspection the two sets of population values
# match exactly except for a rounding difference.

# Country data added in 06/2022 are summarized in additional_country_ssp_data_sources.csv 
# and pulled directly from SspDb_country_data_2013-06-12.csv

################################################################################
## 5. Reference Population Density (refpopdens.csv)

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
