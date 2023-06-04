# ## ------------------------------------------------------------------
# function _folsomPhysiologicalBiomassElemental2015_builder(netid::AbstractString, culid::AbstractString, Di::Int)
    
#     # net
#     net = pull_net(netid)
#     net = _common_format(net)
    
#     # medium from https://doi.org/10.1099/mic.0.000118 
#     cul = pull_cul("folsomPhysiologicalBiomassElemental2015")

#     X = queryfirst(cul, "X", culid, Di; extract = "val")
#     D = queryfirst(cul, "D", culid, Di; extract = "val")
#     c_glc = queryfirst(cul, "c_glc", culid, Di; extract = "val")
#     c_nh4 = queryfirst(cul, "c_nh4", culid, Di; extract = "val")

#     ex_glc_id = extras(net, "EX_GLC")
#     ex_nh4_id = extras(net, "EX_NH4")
#     bounds!(net, ex_glc_id; lb = -c_glc*D/X)
#     bounds!(net, ex_nh4_id; lb = -c_nh4*D/X)
    
#     # under bound biomass
#     objid = extras(net, "BIOM")
#     linear_coefficients!(net, objid, 1.0)
    
#     return net

# end

# ## ------------------------------------------------------------------
# function _register_folsomPhysiologicalBiomassElemental2015()
#     register_network!("folsomPhysiologicalBiomassElemental2015", 
#         _folsomPhysiologicalBiomassElemental2015_builder;
#         use_cache = false,
#         source = "Folsom, James Patrick, and Ross P. Carlson. “Physiological, Biomass Elemental Composition and Proteomic Analyses of Escherichia Coli Ammonium-Limited Chemostat Growth, and Comparison with Iron- and Glucose-Limited Chemostat Growth.” Microbiology (Reading, England) 161, no. 8 (August 2015): 1659–70. https://doi.org/10.1099/mic.0.000118.", 
#         desc = """A contextualized E coli network as described at https://doi.org/10.1099/mic.0.000118. The networks can be set to 3 limiting conditions ['N_Limited', 'C_Limited', 'Fe_Limited'] at 4 dilution rates 'Di ∈ 0:4', Ex: pull_net("folsomPhysiologicalBiomassElemental2015", "ecoli_core", "N_Limited", 1) uses the 'ecoli_cre' as basis and the data from the nitrogen limited culture at D[1]. It uses MetXCultureHub.pull_cul("folsomPhysiologicalBiomassElemental2015") for retriving the experimental data"""
#     )
# end