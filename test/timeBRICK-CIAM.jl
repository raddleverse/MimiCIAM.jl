## Test performance 
using MimiCIAM
using DataDeps

DataDeps.register(DataDeps.DataDep(
    "BRICK outputs",
    """
    Some BRICK output downloaded from https://download.scrim.psu.edu/Wong_etal_BRICK/BRICKms_Wong_etal_2017/output_model/BRICK-model_physical_control_31May2017.nc
    """,
    "https://download.scrim.psu.edu/Wong_etal_BRICK/BRICKms_Wong_etal_2017/output_model/BRICK-model_physical_control_31May2017.nc"))

brickfile = datadep"BRICK outputs/BRICK-model_physical_control_31May2017.nc"
n = 1
low=5
high=95
ystart=2010
yend=2100
tstep=10
rcp=85

# First run 
MimiCIAM.brickCIAM_driver(rcp,brickfile,n,low,high,ystart,yend,tstep)

## Time 1-run version 
@time MimiCIAM.brickCIAM_driver(rcp,brickfile,n,low,high,ystart,yend,tstep)

## Time 10-run version
n=10
@time MimiCIAM.brickCIAM_driver(rcp,brickfile,n,low,high,ystart,yend,tstep)