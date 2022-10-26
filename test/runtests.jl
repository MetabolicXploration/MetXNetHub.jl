# activate env
using TestProject
@activate_testenv

# import
using MetXNetHub
using MetXOptim
import MetXOptim.Clp
using MetXBase
using Test

# test
@testset "MetXNetHub.jl" begin
    
    # to_test = ["linear_net", "toy_net", "ecoli_core", "iJR904", "ECC2", "ECGS"]
    to_test = ["linear_net"]
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
