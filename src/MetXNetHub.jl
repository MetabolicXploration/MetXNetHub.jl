# TODO: fix COBREXA.jl results on the tests
module MetXNetHub
    
    using MetXBase
    using MetXGEMs
    import MetXGEMs.COBREXA
    # using MetXCultureHub # TODO: makes MetXCultureHub works

    using Serialization
    using SparseArrays

    import Scratch

    export pull_net, register_network!, nethub_status

    # globals
    const NETS_REG = Dict{String, Dict{String, Any}}()
    NETS_DIR = ""

    #! include core
    include("core/api.jl")
    include("core/utils.jl")

    #! include nets
    include("nets/ECC2.jl")
    include("nets/ECC2comp.jl")
    include("nets/ECGS.jl")
    include("nets/ENGRO1.jl")
    include("nets/Martinez_Monge_HEK293.jl")
    include("nets/ecoli_core.jl")
    # include("nets/folsomPhysiologicalBiomassElemental2015.jl") # TODO: makes MetXCultureHub works
    include("nets/iCHO2291.jl")
    include("nets/iJO1366.jl")
    include("nets/iJR904.jl")
    include("nets/linear_net.jl")
    include("nets/toy_net.jl")
    include("nets/toy_net4D.jl")
    include("nets/toy_net_cost.jl")

    function __init__()
        # scratch
        global NETS_DIR = Scratch.get_scratch!("nets_dir")

        # register models
        empty!(NETS_REG)
        _register_ecoli_core()
        _register_linear_model()
        _register_toy_net()
        _register_toy_net4D()
        _register_iJR904()
        _register_ECC2()
        # _register_ECC2comp() # TODO: make it growth
        _register_ECGS()
        _register_iJO1366()
        _register_iCHO2291()
        _register_Martinez_Monge_HEK293()
        # _register_folsomPhysiologicalBiomassElemental2015()  # TODO: makes MetXCultureHub works
        _register_ENGRO1()

    end

end