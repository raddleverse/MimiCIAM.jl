## Test performance 
using MimiCIAM
brickfile = "/Users/catherineledna/Downloads/BRICK-model_physical_control_31May2017.nc"
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