module MetXNetHub

using Serialization
using MetXBase
using SparseArrays

import Scratch

export load_net, register_network!, registered_nets


# globals
const NETS_REG = Dict{String, Dict{String, Any}}()
NETS_DIR = ""

include("utils.jl")
include("api.jl")

include("nets/ecoli_core.jl")
include("nets/linear_net.jl")
include("nets/toy_net.jl")
include("nets/iJR904.jl")


function __init__()
    # scratch
    global NETS_DIR = Scratch.get_scratch!("nets_dir")

    # register models
    _register_ecoli_core()
    _register_linear_model()
    _register_toy_net()
    _register_iJR904()

end

end
