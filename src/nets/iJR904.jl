function _iJR904_builder()

    # load
    net = _load_raw_model("iJR904.json")
    net = _common_format(net)
    
    extras!(net, "BIOM", "BIOMASS_Ecoli")
    extras!(net, "EX_GLC", "EX_glc__D_e")
    extras!(net, "EX_NH4", "EX_nh4_e")
    extras!(net, "EX_GLU", "EX_glu__L_e")
    extras!(net, "EX_O2", "EX_o2_e")
    extras!(net, "EX_CO2", "EX_co2_e")
    extras!(net, "ATPM", "ATPM")
    
    linear_weights!(net, "BIOMASS_Ecoli", 1.0)

    return net
end

function _register_iJR904()
    register_network!("iJR904", _iJR904_builder;
        use_cache = false,
        source = "http://bigg.ucsd.edu/models/iJR904", 
        desc = "A genome-scale model of E. coli, just as shipped"
    )
end