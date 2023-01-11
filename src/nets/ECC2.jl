function _ECC2_builder()

    # load
    # 41598_2017_BFsrep39647_MOESM452_ESM
    net = _load_raw_model("ECC2.xml")
    net = MetXBase.dense_vecs(net)

    # elimiate external mets
    ex_mets = filter(metabolites(net)) do met
        endswith(met, "_ex")
    end
    empty_met!(net, ex_mets)  
    # and Biomass met
    empty_met!(net, "Biomass")  
    net = emptyless_model(net)
    
    # open out close in
    ex_rxns = filter(reactions(net)) do rxn
        contains(rxn, "_EX_")
    end
    bounds!(net, ex_rxns; lb = 0.0, ub = 1000.0)

    # minimum medium (open in)
    lb!(net, "R_EX_glc_LPAREN_e_RPAREN_", -10.0)
    lb!(net, "R_EX_ni2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_k_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_fe3_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_so4_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_cl_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_pi_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_mg2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_zn2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_nh4_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_mobd_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_cu2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_mn2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_ca2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_o2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_cobalt2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_fe2_LPAREN_e_RPAREN_", -1000.0)
    
    extras!(net, "BIOM", "R_Ec_biomass_iJO1366_core_53p95M")
    extras!(net, "EX_GLC", "R_EX_glc_LPAREN_e_RPAREN_")
    extras!(net, "EX_NH4", "R_EX_nh4_LPAREN_e_RPAREN_")
    extras!(net, "ATPM", "R_ATPM")

    linear_coefficients!(net, "R_Ec_biomass_iJO1366_core_53p95M", 1.0)

    return net
end

function _register_ECC2()
    register_network!("ECC2", _ECC2_builder;
        source = "Hädicke, O., Klamt, S. EColiCore2: a reference network model of the central metabolism of Escherichia coli and relationships to its genome-scale parent model. Sci Rep 7, 39647 (2017). https://doi.org/10.1038/srep39647", 
        desc = "The EColiCore2 model in a minimum glucose rich medium (max uptake: 10.0 mmol/gDW h), the raw model was extracted from 41598_2017_BFsrep39647_MOESM452_ESM.xls supplementary data"
    )
end