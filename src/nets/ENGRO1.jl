function _ecoli_ENGRO1()

    # load
    # from https://doi.org/10.1371/journal.pcbi.1005758.s012
    net0 = _load_raw_model("ENGRO1.xml")

    # @show size(net0)

    # from https://doi.org/10.1371/journal.pcbi.1005758.s011
    src_ids = ["hexokinase", "G6P_isomerase", "PFK", "aldolase", "GAPDH", "PGK", "PGM", "enolase", "pyr_kinase", "LDH", "pyr_carboxylase", "pep_carboxykinase", "pyr_dehydrogenase", "cit_synthase", "aconitase", "IDH2_nadp", "IDH2_nad", "AKG_dehydrogenase", "SuCoA_synthase", "Succ_dehydrogenase", "fumarase", "malate_dehydrogenase", "malic_enzyme_NADP", "malic_enzyme_NAD", "ASPTA", "ALATA_L", "GLUDxm", "GLUDym", "Cit_lyase", "FFAsynthesis", "glutaminase", "gln_synthetase", "CPS1", "OTC", "ASS", "ASL", "ARG1", "ORNTArm", "G5SADrm", "P5CRm", "P5CRxm", "ORNDC", "ATP_ADP", "G6P_dehydrogenase", "gluconolactonase", "6PDG_dehydrogenase", "Ru5P_isomerase", "PRPPS", "SHMT1", "MTHFD1", "ser_synthesis", "RESP1", "RESP2", "Hcys_synthesis", "CYSTS", "CYSTGL", "GLUCYS", "GTHS", "GTHP", "GTHO", "glutamyltrasferase", "asparagine_synthetase", "CYSTA", "ASRGL1", "CAD", "SDS", "GATM", "EX_O2", "EX_Gln", "EX_Glc", "EX_Arg", "EX_THF", "EX_Met", "protein_synthesis", "biomass_synthesis", "DM_LACT", "DM_H2O2", "DM_NH3", "DM_Urea", "DM_PRPP", "DM_biomass", "DM_Putrescine", "DM_forTHF", "DM_GSSG"]
    lbs =  [0, 0, 0, 0, 0, 0, -1000, -1000, 0, 0, 0, 0, 0, 0, -1000, -1000, -1000, 0, -1000, -1000, -1000, -1000, 0, 0, -1000, -1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1000, 0, -1000, -1000, 0, 0, 0, 0, 0, 0, 0, -1000, -1000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -38, -40, -10, -20, -20, -20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ubs = [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]
    # @show size(src_ids)
    # @show size(lbs)
    # @show size(ubs)

    for (i, id) in enumerate(src_ids)
        id = string("R_", id)
        if id in net0.rxns  
            lb!(net0, id, lbs[i])
            ub!(net0, id, ubs[i])
            continue
        end
        @show id
    end

    # Eliminate the boundaries
    boundary_mets = findall(net0.mets) do id
        endswith(id, "b")
    end
    # @show length(boundary_mets)
    empty_met!(net0, boundary_mets)
    net = emptyless_model(net0)
    net = _common_format(net)

    # rename DMs
    net.rxns .= map(net.rxns) do rxn
        replace(rxn, r"\AR_DM_" => "R_EX_")
    end

    # Close intakes
    ex_rxns = findall(net.rxns) do id
        startswith(id, "R_EX_")
    end
    bounds!(net, ex_rxns, 0.0, 1000.0)
    # @show length(ex_rxns)
    
    # Open medium
    lb!(net, "R_EX_Arg", -1000.0)
    lb!(net, "R_EX_Met", -1000.0)
    lb!(net, "R_EX_Glc", -10.0)
    lb!(net, "R_EX_Gln", -10.0)
    lb!(net, "R_EX_O2", -1000.0)
    lb!(net, "R_EX_THF", -1000.0)
    
    # Extras
    extras!(net, "BIOM", "R_biomass_synthesis")
    extras!(net, "EX_GLC", "R_EX_Glc")
    extras!(net, "EX_GLN", "R_EX_Gln")
    extras!(net, "EX_NH3", "EX_nh3_e")
    extras!(net, "EX_O2", "R_EX_O2")
    extras!(net, "ATPM", "R_ATP_ADP")

    # default linear obj
    linear_coefficients!(net, "R_biomass_synthesis", 1.0)

    return net
end

function _register_ENGRO1()
    register_network!("ENGRO1", _ecoli_ENGRO1;
        use_cache = true,
        source = "https://doi.org/10.1371/journal.pcbi.1005758 S", 
        desc = "A core model of a human cell"
    )
end