# TODO: wait till network manipulation is ready at MetXBase

# function _toy_net_cost_builder(resp_cost = -0.1, E_demand = 5.0)
#     @assert resp_cost <= 0.0

#     net = load_net("toy_net")
#     net = _common_format(net)

#     # Atpm
#     lb!(net, "atpm", E_demand)

#     # cost
#     # Add a new metabolite simulating the cost penalazing 
#     # reaction fluxes. 
#     # A new balance equations is then added:
#     #     Σ(rᵢ*costᵢ) + tot_cost = 0
#     # Because the cost coefficients (costᵢ) < 0, the system must allocate 
#     # the fluxes (rᵢ) so that Σ(rᵢ*costᵢ) = tot_cost, and tot_cost
#     # are usually bounded [0.0, 1.0]

#     # We will add 1 met and 1 rxn
#     M, N = size(net)
#     net = expanded_model(net, M + 1, N + 1)
#     # Add cost met
#     set_met!(net, M + 1, Met("cost"; rxns = ["resp"], S = [resp_cost]))
#     set_rxn!(net, N + 1, Rxn("tot_cost"; mets = ["cost"], S = [1.0], lb = 0.0, ub = 1.0))

#     # extras
#     extras!(net, "BIOM", "biom")
#     extras!(net, "COST", "tot_cost")

#     return net
# end

# function _register__toy_net_cost()
#     register_network!("_toy_net_cost", _toy_net_cost_builder;
#         use_cache = false,
#         source = "Home made ;)", 
#         desc = "A model with M metabolites in a linear pathway:  -> Met_1 -> ... -> Met_M -> "
#     )
# end