# Tests to check sensible output from the MimiCIAM model

# Original code: ~2020 Catherine Ledna
# Modified code: 1 Apr 2021 Tony Wong

#------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------

using MimiCIAM
using Test

@testset "MimiCIAM" begin
    m = MimiCIAM.get_model()
    run(m)

    isneg(x)=length(x[isless.(x,0)])
    
    vargroup2D = [:WetlandNoAdapt,:FloodNoAdapt,:StormCapitalNoAdapt,:StormPopNoAdapt,:RelocateNoAdapt,
                :NoAdaptCost,:OptimalCost,:OptimalStormCapital,:OptimalStormPop,:OptimalConstruct,
                :OptimalWetland,:OptimalFlood,:OptimalRelocate,:WetlandRetreat,:WetlandProtect]
    vargroup3D = [:Construct,:StormCapitalProtect,:StormPopProtect,:StormCapitalRetreat,
                :StormPopRetreat,:FloodRetreat,:RelocateRetreat,:RetreatCost,:ProtectCost]
    vargroup4 = [:OptimalH,:OptimalR,:DryLandLossOptimal]

    for v in vargroup2D
        var = m[:slrcost,v]
        @test isneg(var)==0
                    
    end

    for v in vargroup4
        var=m[:slrcost,v]
        tsteps=m[:slrcost,:ntsteps]
        for i in 2:tsteps
            @test sum(var[i,:].< var[i-1,:])==0

        end

    end
            
    for v in vargroup3D
        var=m[:slrcost,v]
        levels= size(var)[3]
        for k in 1:levels
            var2 = var[:,:,k]
            @test isneg(var2)==0
        end
    end

end

