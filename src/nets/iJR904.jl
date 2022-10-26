function _iJR904_builder()

    # load
    net = _load_raw_model("iJR904.xml")
    net = MetXBase.dense_vecs(net)
    
    set_extra!(net, "BIOM", "R_BIOMASS_Ecoli")
    set_extra!(net, "EX_GLC", "R_EX_glc__D_e")
    
    lin_objective!(net, "R_BIOMASS_Ecoli", 1.0)

    return net
end

function _register_iJR904()
    register_network!("iJR904", _iJR904_builder;
        use_cache = true,
        source = "http://bigg.ucsd.edu/models/iJR904", 
        desc = "A genome-scale model of E. coli, just as shipped"
    )
end