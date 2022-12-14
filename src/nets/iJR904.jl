function _iJR904_builder()

    # load
    net = _load_raw_model("iJR904.xml")
    net = MetXBase.dense_vecs(net)
    
    extras!(net, "BIOM", "R_BIOMASS_Ecoli")
    extras!(net, "EX_GLC", "R_EX_glc__D_e")
    extras!(net, "EX_NH4", "R_EX_nh4_e")
    extras!(net, "ATPM", "R_ATPM")
    
    linear_coefficients!(net, "R_BIOMASS_Ecoli", 1.0)

    return net
end

function _register_iJR904()
    register_network!("iJR904", _iJR904_builder;
        use_cache = true,
        source = "http://bigg.ucsd.edu/models/iJR904", 
        desc = "A genome-scale model of E. coli, just as shipped"
    )
end