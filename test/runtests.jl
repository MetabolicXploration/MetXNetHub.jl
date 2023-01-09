# activate env
using RunTestsEnv
@activate_testenv

# import stuff from `test/env/Project`.
# `MetXOptim` depends on `MetXNetHub` but `MetXOptim` is needed for testing (fba)
# so it can't be a direct dependency of `MetXNetHub`. 
# This is annoying while developing because `MetXOptim` will be `dev`ed
# `RunTestsEnv` (by running the test in an independent project) is avoiding the circular issue
using MetXNetHub
using MetXOptim 
import MetXOptim.Clp
using MetXBase
using Test

# test
@testset "MetXNetHub.jl" begin
    
    to_test = [
        "linear_net", "toy_net", "toy_net4D", 
        "ecoli_core", "ECC2", "ECGS", 
        "iJR904", "iJO1366", "Martinez_Monge_HEK293"
    ]
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
            biom_id = extras(net, "BIOM")
            glcex_id = extras(net, "EX_GLC")
            linear_coefficients!(net, biom_id, 1.0)

            opm = fba(net, Clp.Optimizer)
            objval = solution(opm, biom_id)
            glcval = solution(opm, glcex_id)

            @info("Done", id, objval, glcval)
            @test objval > 0
        end
    end
end
