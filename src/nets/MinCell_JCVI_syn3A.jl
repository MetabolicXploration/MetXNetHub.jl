# MinCell_JCVI_syn3A
function _MinCell_JCVI_syn3A_builder()

    # load
    # Hi, I'm trying to load a json model (download link: elife-36842-supp10-v2.json) from this paper.
    # TODO: Wait till COBREXA loading issue is resolved (PR: https://github.com/LCSB-BioCore/COBREXA.jl/pull/822)
    # to load the original file "elife-36842-supp10-v2.json"
    # net = _load_raw_model("elife-36842-supp10-v2.json")
    net = _load_raw_model("FORMATED-elife-36842-supp10-v2.json")
    net = _common_format(net)

    # TODO:define medium
    # open out close in
    # ex_rxns = filter(reactions(net)) do rxn
    #     contains(rxn, "_EX_")
    # end
    # bounds!(net, ex_rxns; lb = 0.0, ub = 1000.0)

    # # medium
    # lb!(net, "EX_glc__D_e", -10.0)
    # lb!(net, "EX_k_e", -1000.0)
    # lb!(net, "EX_cl_e", -1000.0)
    # lb!(net, "EX_ni2_e", -1000.0)
    # lb!(net, "EX_pi_e", -1000.0)
    # lb!(net, "EX_zn2_e", -1000.0)
    # lb!(net, "EX_nh4_e", -1000.0)
    # lb!(net, "EX_mn2_e", -1000.0)
    # lb!(net, "EX_o2_e", -1000.0)
    # lb!(net, "EX_fe2_e", -1000.0)
    # lb!(net, "EX_cu2_e", -1000.0)
    # lb!(net, "EX_ca2_e", -1000.0)
    # lb!(net, "EX_so4_e", -1000.0)
    # lb!(net, "EX_mobd_e", -1000.0)
    # lb!(net, "EX_mg2_e", -1000.0)
    # lb!(net, "EX_cobalt2_e", -1000.0)
    
    # extras!(net, "BIOM", "BIOMASS")
    # extras!(net, "EX_GLC", "EX_glc__D_e")
    # extras!(net, "EX_NH4", "EX_nh4_e")
    # extras!(net, "EX_GLU", "EX_glu__L_e")
    # extras!(net, "EX_O2", "EX_o2_e")
    # extras!(net, "EX_CO2", "EX_co2_e")
    # extras!(net, "ATPM", "ATPM")

    linear_weights!(net, "BIOMASS", 1.0)

    return net
end

function _register_MinCell_JCVI_syn3A()
    register_network!("MinCell_JCVI_syn3A", _MinCell_JCVI_syn3A_builder;
        use_cache = true,
        source = "Breuer, Marian, Tyler M Earnest, Chuck Merryman, Kim S Wise, Lijie Sun, Michaela R Lynott, Clyde A Hutchison, et al. “Essential Metabolism for a Minimal Cell.” eLife 8 (January 18, 2019): e36842. https://doi.org/10.7554/eLife.36842.", 
        download = "https://cdn.elifesciences.org/articles/36842/elife-36842-supp10-v3.json.zip",
        desc = "The raw JCVI-syn3A model of a minimal cell in a (TODO: define medium)"
    )
end