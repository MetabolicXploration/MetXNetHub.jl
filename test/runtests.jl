using MetXNetHub
using Test

@testset "MetXNetHub.jl" begin
    
    # tesr load and cache
    to_test = ["linear_net", "toy_net", "ecoli_core", "iJR904"]
    build_args = [(10,), (), ()]
    
    for (id, args) in zip(to_test, build_args)
        MetXNetHub.clear_cache!(id)
        net0 = MetXNetHub.load_net(id, args...; clear_cache = false)
        net1 = MetXNetHub.load_net(id, args...; clear_cache = false)
        @test net0 == net1
        MetXNetHub.clear_cache!(id)
    end

end
