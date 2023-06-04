# Taked from "Martínez-Monge, Iván, Joan Albiol, Martí Lecina, Leticia Liste-Calleja, Joan Miret, Carles Solà, and Jordi J. Cairó. “Metabolic Flux Balance Analysis during Lactate and Glucose Concomitant Consumption in HEK293 Cell Cultures.” Biotechnology and Bioengineering 116, no. 2 (2019): 388–404. https://doi.org/10.1002/bit.26858." supplementary materials. 

function _Martinez_Monge_HEK293_builder()

    # load 
    rawfilename = "bit26858-sup-0001-S1_Model.xml"
    rawfile = _raw_file(rawfilename)

    std_net = COBREXA.load_model(COBREXA.StandardModel, rawfile)

    # TODO: add gene associations
    # TODO: add notes to MetNet
    # "<notes>\n  <html:p>GENE_ASSOCIATION:  </html:p>\n  <html:p>PROTEIN_ASSOCIATION:  </html:p>\n  <html:p>PROTEIN_CLASS:  </html:p>\n  <html:p>SUBSYSTEM:  Transport, mitochondrial</html:p>\n  <html:p>CONFIDENCE_LEVEL: 2</html:p>\n  <html:p>PROTEIN_ASSOCIATION: </html:p>\n  <html:p>NOTES: NCD</html:p>\n</notes>"
    subSystems = String[]
    subSystems_reg = r"<html:p>SUBSYSTEM:(?<subsys>.*)</html:p>"
    for (rxn, dat) in std_net.reactions
        notes = dat.notes[""][1]
        m = match(subSystems_reg, string(notes))
        subsys = isnothing(m) ? "" : strip(m[:subsys])
        push!(subSystems, subsys)
    end
    
    net = convert(MetNet, std_net)
    net = MetNet(net; subSystems)
    net = _common_format(net)

    # elimiate external mets
    ex_mets = filter(metabolites(net)) do met
        endswith(met, "_b")
    end
    empty_met!(net, ex_mets)  
    net = emptyless_model(net)

    # close intakes
    exchs = filter(reactions(net)) do rxn
        startswith(rxn, "R_EX_")
    end
    lb!(net, exchs, 0.0)

    # add medium
    # TODO: find a realistic medium
    lb!(net, "R_EX_tyr_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_his_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_asn_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_arg_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_ile_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_phe_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_met_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_thr_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_lys_L_LPAREN_e_RPAREN_", -0.05)
    lb!(net, "R_EX_leu_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_val_L_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_trp_L_LPAREN_e_RPAREN_", -1000.0)
    
    # lb!(net, "R_EX_glu_L_LPAREN_e_RPAREN_", -1000.0)
    # lb!(net, "R_EX_gln_L_LPAREN_e_RPAREN_", -1000.0)
    # lb!(net, "R_EX_glc_LPAREN_e_RPAREN_", -1000.0)
    # lb!(net, "R_EX_lac_L_LPAREN_e_RPAREN_", -1000.0)

    lb!(net, "R_EX_o2_LPAREN_e_RPAREN_", -1000.0)

    lb!(net, "R_EX_h_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_pi_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_h2o_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_co2_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_nh4_LPAREN_e_RPAREN_", -1000.0)
    lb!(net, "R_EX_hco3_LPAREN_e_RPAREN_", -1000.0)

    # extras
    extras!(net, "BIOM", "R_Ex_Biomass")
    
    extras!(net, "EX_GLC", "R_EX_glc_LPAREN_e_RPAREN_")
    extras!(net, "EX_NH4", "R_EX_nh4_LPAREN_e_RPAREN_")

    extras!(net, "EX_TYR", "R_EX_tyr_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_HIS", "R_EX_his_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_ASN", "R_EX_asn_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_ARG", "R_EX_arg_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_ILE", "R_EX_ile_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_PHE", "R_EX_phe_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_MET", "R_EX_met_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_THR", "R_EX_thr_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_LYS", "R_EX_lys_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_LEU", "R_EX_leu_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_VAL", "R_EX_val_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_TRP", "R_EX_trp_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_GLU", "R_EX_glu_L_LPAREN_e_RPAREN_")
    extras!(net, "EX_GLN", "R_EX_gln_L_LPAREN_e_RPAREN_")
    
    extras!(net, "EX_LAC", "R_EX_lac_L_LPAREN_e_RPAREN_")
    
    extras!(net, "EX_CO2", "R_EX_co2_LPAREN_e_RPAREN_")
    extras!(net, "EX_O2", "R_EX_o2_LPAREN_e_RPAREN_")

    extras!(net, "ATPM", "R_DM_atp_c_")

    # objective
    linear_coefficients!(net, "R_Ex_Biomass", 1.0)

    return net

end

function _register_Martinez_Monge_HEK293()
    register_network!("Martinez_Monge_HEK293", _Martinez_Monge_HEK293_builder;
        use_cache = true,
        source = "https://doi.org/10.1002/bit.26858", 
        desc = "A genome-scale model of HEK293 (335, 354), derived from Recon2"
    )
end