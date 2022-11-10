function _ecoli_core_builder()

    # load
    net = _load_raw_model("e_coli_core.mat")
    net = MetXBase.dense_vecs(net)
    
    extras!(net, "BIOM", "BIOMASS_Ecoli_core_w_GAM")
    extras!(net, "EX_GLC", "EX_glc__D_e")

    lin_objective!(net, "BIOMASS_Ecoli_core_w_GAM", 1.0)

    return net
end

function _register_ecoli_core()
    register_network!("ecoli_core", _ecoli_core_builder;
        source = "http://bigg.ucsd.edu/models/e_coli_core", 
        desc = "The clasical core model of E. coli in a minimum glucose rich medium"
    )
end