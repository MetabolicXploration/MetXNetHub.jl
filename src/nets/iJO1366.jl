function _iJO1366_builder()

    # load
    net = MetXNetHub._load_raw_model("iJO1366.xml")

    # open out close in
    ex_rxns = filter(reactions(net)) do rxn
        contains(rxn, "_EX_")
    end
    bounds!(net, ex_rxns; lb = 0.0, ub = 1000.0)

    # medium
    lb!(net, "R_EX_glc__D_e", -10.0)
    lb!(net, "R_EX_k_e", -1000.0)
    lb!(net, "R_EX_cl_e", -1000.0)
    lb!(net, "R_EX_ni2_e", -1000.0)
    lb!(net, "R_EX_pi_e", -1000.0)
    lb!(net, "R_EX_zn2_e", -1000.0)
    lb!(net, "R_EX_nh4_e", -1000.0)
    lb!(net, "R_EX_mn2_e", -1000.0)
    lb!(net, "R_EX_o2_e", -1000.0)
    lb!(net, "R_EX_fe2_e", -1000.0)
    lb!(net, "R_EX_cu2_e", -1000.0)
    lb!(net, "R_EX_ca2_e", -1000.0)
    lb!(net, "R_EX_so4_e", -1000.0)
    lb!(net, "R_EX_mobd_e", -1000.0)
    lb!(net, "R_EX_mg2_e", -1000.0)
    lb!(net, "R_EX_cobalt2_e", -1000.0)
    
    set_extra!(net, "BIOM", "R_BIOMASS_Ec_iJO1366_core_53p95M")
    set_extra!(net, "EX_GLC", "R_EX_glc__D_e")

    lin_objective!(net, "R_BIOMASS_Ec_iJO1366_core_53p95M", 1.0)

    return net
end

function _register_iJO1366()
    register_network!("iJO1366", _iJO1366_builder;
        use_cache = true,
        source = "Orth JD, Conrad TM, Na J, Lerman JA, Nam H, Feist AM, Palsson BØ. A comprehensive genome-scale reconstruction of Escherichia coli metabolism--2011. Mol Syst Biol. 2011 Oct 11;7:535. doi: 10.1038/msb.2011.65. PMID: 21988831; PMCID: PMC3261703.", 
        download = "http://bigg.ucsd.edu/models/iJO1366",
        desc = "The raw iJO1366 model of E Coli in a minimum glucose rich medium (max uptake: 10.0 mmol/gDW h)"
    )
end