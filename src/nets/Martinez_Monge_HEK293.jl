# Taked from "Martínez-Monge, Iván, Joan Albiol, Martí Lecina, Leticia Liste-Calleja, Joan Miret, Carles Solà, and Jordi J. Cairó. “Metabolic Flux Balance Analysis during Lactate and Glucose Concomitant Consumption in HEK293 Cell Cultures.” Biotechnology and Bioengineering 116, no. 2 (2019): 388–404. https://doi.org/10.1002/bit.26858." supplementary materials. 

function _Martinez_Monge_HEK293_builder()

    rawfilename = "bit26858-sup-0001-S1_Model.xml"
    net = MetXNetHub._load_raw_model(rawfilename)

    # elimiate external mets
    ex_mets = filter(metabolites(net)) do met
        endswith(met, "_b")
    end
    empty_met!(net, ex_mets)  
    net = emptyless_model(net)

    # close intakes
    global exchs = filter(reactions(net)) do rxn
        startswith(rxn, "R_EX_")
    end
    lb!(net, exchs, 0.0)

    # add medium
    # TODO: find a realistic medium
    lb!(net, "R_EX_tyr_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_his_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_asn_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_arg_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_ile_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_phe_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_met_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_thr_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_lys_L_LPAREN_e_RPAREN_", -0.05)
    lb!(net, "R_EX_leu_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_val_L_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_trp_L_LPAREN_e_RPAREN_", -10_000)
    
    # lb!(net, "R_EX_glu_L_LPAREN_e_RPAREN_", -10_000)
    # lb!(net, "R_EX_gln_L_LPAREN_e_RPAREN_", -10_000)
    # lb!(net, "R_EX_glc_LPAREN_e_RPAREN_", -10_000)
    # lb!(net, "R_EX_lac_L_LPAREN_e_RPAREN_", -10_000)

    lb!(net, "R_EX_o2_LPAREN_e_RPAREN_", -10_000)

    lb!(net, "R_EX_h_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_pi_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_h2o_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_co2_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_nh4_LPAREN_e_RPAREN_", -10_000)
    lb!(net, "R_EX_hco3_LPAREN_e_RPAREN_", -10_000)

    # extras
    extras!(net, "BIOM", "R_Ex_Biomass")
    extras!(net, "EX_GLC", "R_EX_glc_LPAREN_e_RPAREN_")
    extras!(net, "ATPM", "R_DM_atp_c_")

    # objective
    lin_objective!(net, "R_Ex_Biomass", 1.0)

    return net

end

function _register_Martinez_Monge_HEK293()
    register_network!("Martinez_Monge_HEK293", _Martinez_Monge_HEK293_builder;
        use_cache = true,
        source = "https://doi.org/10.1002/bit.26858", 
        desc = "A genome-scale model of HEK293 (335, 354), derived from Recon2"
    )
end