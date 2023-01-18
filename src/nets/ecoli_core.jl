function _ecoli_core_builder()

    # load
    net = _load_raw_model("e_coli_core.mat")
    net = MetXBase.dense_vecs(net)
    
    extras!(net, "BIOM", "BIOMASS_Ecoli_core_w_GAM")
    extras!(net, "EX_GLC", "EX_glc__D_e")
    extras!(net, "EX_NH4", "EX_nh4_e")
    extras!(net, "EX_GLU", "EX_glu__L_e")
    extras!(net, "ATPM", "ATPM")

    linear_coefficients!(net, "BIOMASS_Ecoli_core_w_GAM", 1.0)

    return net
end

function _register_ecoli_core()
    register_network!("ecoli_core", _ecoli_core_builder;
        source = "http://bigg.ucsd.edu/models/e_coli_core", 
        desc = "The clasical core model of E. coli in a minimum glucose rich medium"
    )
end