function _iCHO2291_builder()

    # load
    net = _load_raw_model("iCHO2291.xml")
    net = _common_format(net)

    # Close intakes open productions
    for rxn in colids(net)
        startswith(rxn, "R_EX_") || continue
        bounds!(net, rxn, 0.0, 1000.0)
    end

    # Close sinks
    for rxn in colids(net)
        startswith(rxn, "R_SK_") || continue
        bounds!(net, rxn, 0.0, 0.0)
    end

    # Mminimum Medium
    bounds!(net, "R_EX_his_L__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_ile_L__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_leu_L__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_lys_L__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_met_L__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_phe_L__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_thr_L__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_trp_L__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_val_L__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_Tyr_ggn__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_o2__40__e__41__", -1000.0,  1000.0)
    bounds!(net, "R_EX_pi__40__e__41__", -1000.0,  1000.0)
    
    # bounds!(net, "R_EX_tyr_L__40__e__41__", -1.0, 1000.0)

    # extras
    bounds!(net, "R_EX_glc__40__e__41__", -1.0,  1000.0)
    bounds!(net, "R_EX_gln_L__40__e__41__", -1.0,  1000.0)

    # bounds!(net, "R_EX_chol__40__e__41__", -1000.0,  1000.0)
    # bounds!(net, "R_EX_hco3__40__e__41__", -1000.0,  1000.0)
    # bounds!(net, "R_EX_hxan__40__e__41__", -1000.0,  1000.0)
    # bounds!(net, "R_EX_arg_L__40__e__41__", -1000.0,  1000.0)
    # bounds!(net, "R_EX_asn_L__40__e__41__", -1000.0,  1000.0)
    # bounds!(net, "R_EX_asp_L__40__e__41__", -1000.0,  1000.0)
    # bounds!(net, "R_EX_cys_L__40__e__41__", -1000.0,  1000.0)
    # bounds!(net, "R_EX_pro_L__40__e__41__", -1000.0,  1000.0)
    # bounds!(net, "R_EX_h__40__e__41__", -1000.0,  1000.0)
    # bounds!(net, "R_EX_ser_L__40__e__41__", -1000.0,  1000.0)
    
    # bounds!(net, "R_EX_co2__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_h2o2__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_h2o__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_h2s__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_nh4__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_urea__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_lac_L__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_lac_D__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_ac__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_gly__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_for__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_ala_D__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_ala_L__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_ser_D__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_ser_L__40__e__41__", 0.0, 1000.0)
    # bounds!(net, "R_EX_pyr__40__e__41__", 0.0, 1000.0)
    
     # Demands
    bounds!(net, "R_DM_atp__91__c__93__", 0.0, 1000.0)
    bounds!(net, "R_DM_Ser__47__Thr__91__ly__93__", 0.0, 0.0)
    
    # sinks
    bounds!(net, "R_SK_Tyr_ggn__91__c__93__", -0.1, 1000.0)

    extras!(net, "BIOM", "R_biomass_cho")
    extras!(net, "EX_GLC", "R_EX_glc__40__e__41__")
    extras!(net, "EX_NH4", "R_EX_nh4__40__e__41__")
    extras!(net, "EX_GLU", "R_EX_glu_L__40__e__41__")
    extras!(net, "EX_O2", "R_EX_o2__40__e__41__")
    extras!(net, "EX_CO2", "R_EX_co2__40__e__41__")
    extras!(net, "ATPM", "R_DM_atp__91__c__93__")

    linear_weights!(net, "R_biomass_cho", 1.0)

    return net

end

function _register_iCHO2291()
    register_network!("iCHO2291", _iCHO2291_builder;
        use_cache = true,
        source = "Yeo, Hock Chuan, Jongkwang Hong, Meiyappan Lakshmanan, and Dong-Yup Lee. “Enzyme Capacity-Based Genome Scale Modelling of CHO Cells.” Metabolic Engineering 60 (July 1, 2020): 138-47. https://doi.org/10.1016/j.ymben.2020.04.005.", 
        download = "https://www.ebi.ac.uk/biomodels/MODEL1912180001",
        desc = "The iCHO2291 model was reconstructed by updating the previously published iCHO1766 model in a 6-step procedure: 1) identification and removal of duplicate metabolites and reactions, 2) replacement of lumped reactions into detailed steps and removal of biochemically inconsistent reactions, 3) update of GPR based on latest genome annotations (as on Nov 1, 2018), 4) correction of GPR and reaction compartment assignment based on subcellular localization, 5) inclusion of new reactions based on new genome annotations, and 6) metabolic gap identification and resolution."
    )
end