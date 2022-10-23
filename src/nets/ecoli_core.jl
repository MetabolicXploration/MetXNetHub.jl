function _ecoli_core_builder()

    # load
    net = _load_raw_model("e_coli_core.mat")
    net = _dense_vecs(net)
    
    set_extra!(net, "BIOM", "BIOMASS_Ecoli_core_w_GAM")
    set_extra!(net, "GLC_EX", "EX_glc__D_e")

    return net
end

function _register_ecoli_core()
    register_network!("ecoli_core", _ecoli_core_builder;
        source = "http://bigg.ucsd.edu/models/e_coli_core", 
        desc = "The clasical core model of E. coli in a minimum glucose rich medium"
    )
end