# activate env
using RunTestsEnv
@activate_testenv

# import stuff from `test/env/Project`.
# `MetXOptim` depends on `MetXNetHub` but `MetXOptim` is needed for testing (fba)
# so it can't be a direct dependency of `MetXNetHub`. 
# This is annoying while developing because `MetXOptim` will be `dev`ed
# `RunTestsEnv` (by running the test in an independent project) is avoiding the circular issue
using MetXGEMs
using MetXNetHub
using MetXOptim 
import MetXOptim.Clp
using MetXBase
using Test

# test
@testset "MetXNetHub.jl" begin
    
    to_test = [
        
        # Toy models
        "linear_net", "toy_net", "toy_net4D", 

        # E coli
        "ecoli_core", "ECC2", "ECGS", 
        "iJR904", "iJO1366", 
        # "folsomPhysiologicalBiomassElemental2015", # TODO: make MetXCultureHub work
        
        # HEK
        "Martinez_Monge_HEK293",

        # CHO
        "iCHO2291",

        # Human
        "ENGRO1", 
        "SysBioChalmers_Human_GEM", 
        "SysBioChalmers_EnzymeConstrained_humanModels"
    ]

    # load args
    build_args = Dict{String, Vector{Any}}()
    build_args["linear_net"] = [(10,)]
    build_args["SysBioChalmers_EnzymeConstrained_humanModels"] = #...
        tuple.(MetXNetHub.__SysBioChalmers_EnzymeConstrained_humanModels_gemids)
    push!(build_args["SysBioChalmers_EnzymeConstrained_humanModels"], 
        (MetXNetHub.__SysBioChalmers_EnzymeConstrained_humanModels_gemids[1], "AGE1.HN.AAT"),
        (MetXNetHub.__SysBioChalmers_EnzymeConstrained_humanModels_gemids[1], "AGE1.HN"),
    )

    build_args["folsomPhysiologicalBiomassElemental2015"] = [
        ("ecoli_core", limid, Di) 
            for Di in 1:4 
                for limid in ["$(nut)_Limited" for nut in ["N", "C", "Fe"]]
    ]

    # ------------------------------------------------
    # TODO: No use case yet for caching. Loading is 'fast'
    # let

    #     println()
    #     println("="^60)
    #     println("LOAD AND CACHE")
    #     println("."^60)
    #     println()

    #     MetXNetHub.clear_cache!()
    #     for id in to_test
            
    #         nethub_status(id)
    #         println()

    #         argsv = get(build_args, id, [()])::Vector{<:Tuple}

    #         for args in argsv
    #             net0 = pull_net(id, args...; clear_cache = false)
    #             net1 = pull_net(id, args...; clear_cache = false)
    #             @test net0 == net1
    #         end
    #     end
    # end

    # ------------------------------------------------
    let
        println()
        println("="^60)
        println("FBA TEST")
        println("."^60)
        println()
    
        for id in to_test
        
            argsv = get(build_args, id, [()])

            for args in argsv
                
                @time begin
                    net = pull_net(id, args...)
                    obj_id = extras(net, "BIOM")
                    _size = size(net)

                    # Test biom is set
                    @test linear_coefficients(net, obj_id) == 1.0

                    opm = FBAOpModel(net, Clp.Optimizer)
                    net = nothing; GC.gc(); # free (just in case)
                    optimize!(opm)
                    objval = solution(opm, obj_id)

                    @info("DONE!!!", id, size = _size, obj_id, objval, args)
                    @test objval > 0
                end

            end # for args in argsv
        end # for id in to_test
    end

end
