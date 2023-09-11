function _iJO1366_builder()

    # load
    net = _load_raw_model("iJO1366.json")
    net = _common_format(net)

    # open out close in
    ex_rxns = filter(reactions(net)) do rxn
        contains(rxn, "_EX_")
    end
    bounds!(net, ex_rxns; lb = 0.0, ub = 1000.0)

    # medium
    lb!(net, "EX_glc__D_e", -10.0)
    lb!(net, "EX_k_e", -1000.0)
    lb!(net, "EX_cl_e", -1000.0)
    lb!(net, "EX_ni2_e", -1000.0)
    lb!(net, "EX_pi_e", -1000.0)
    lb!(net, "EX_zn2_e", -1000.0)
    lb!(net, "EX_nh4_e", -1000.0)
    lb!(net, "EX_mn2_e", -1000.0)
    lb!(net, "EX_o2_e", -1000.0)
    lb!(net, "EX_fe2_e", -1000.0)
    lb!(net, "EX_cu2_e", -1000.0)
    lb!(net, "EX_ca2_e", -1000.0)
    lb!(net, "EX_so4_e", -1000.0)
    lb!(net, "EX_mobd_e", -1000.0)
    lb!(net, "EX_mg2_e", -1000.0)
    lb!(net, "EX_cobalt2_e", -1000.0)
    
    extras!(net, "BIOM", "BIOMASS_Ec_iJO1366_core_53p95M")
    extras!(net, "EX_GLC", "EX_glc__D_e")
    extras!(net, "EX_NH4", "EX_nh4_e")
    extras!(net, "EX_GLU", "EX_glu__L_e")
    extras!(net, "EX_O2", "EX_o2_e")
    extras!(net, "EX_CO2", "EX_co2_e")
    extras!(net, "ATPM", "ATPM")

    linear_weights!(net, "BIOMASS_Ec_iJO1366_core_53p95M", 1.0)

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