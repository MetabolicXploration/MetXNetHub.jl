# activate env
using DevTests
@activate_testenv

# import stuff from `test/env/Project`
# `MetXOptim` depends on `MetXNetHub` but `MetXOptim` is needed for testing (fba)
# so it can't be a direct dependency of `MetXNetHub`. 
# `DevTests` (by running the test in an independent project) is avoiding the circular issue
using MetXNetHub
using MetXOptim 
import MetXOptim.Clp
using MetXBase
using Test

# test
@testset "MetXNetHub.jl" begin
    
    to_test = ["linear_net", "toy_net", "ecoli_core", "iJR904", "ECC2", "ECGS"]
    build_args = Dict()
    build_args["linear_net"] = (10,)
    
    # ------------------------------------------------
    let
        println()
        println("="^60)
        println("LOAD AND CACHE")
        println("."^60)
        println()

        MetXNetHub.clear_cache!()
        for id in to_test
            args = get(build_args, id, ())
            net0 = MetXNetHub.pull_net(id, args...; clear_cache = false)
            net1 = MetXNetHub.pull_net(id, args...; clear_cache = false)
            @info("Done", id)
            @test net0 == net1
        end
    end

    # ------------------------------------------------
    let
        println()
        println("="^60)
        println("FBA TEST")
        println("."^60)
        println()
    
        for id in to_test
        
            args = get(build_args, id, ())
            net = MetXNetHub.pull_net(id, args...)
            biom_id = get_extra(net, "BIOM")
            glcex_id = get_extra(net, "EX_GLC")
            lin_objective!(net, biom_id)

            opm = fba(net, Clp.Optimizer)
            objval = solution(opm, biom_id)
            glcval = solution(opm, glcex_id)

            @info("Done", id, objval, glcval)
            @test objval > 0
        end
    end
end
