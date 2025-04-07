function _iML1515_builder()

    # load
    net = _load_raw_model("iML1515.json")
    net = _common_format(net)

    # close all intakes
    for rxn in reactions(net)
        startswith(rxn, "EX_") || continue
        bounds!(net,  rxn, 0.0, 1000.00)
    end

    # open sinks
    for rxn in reactions(net)
        startswith(rxn, "DM_") || continue
        bounds!(net, rxn, 0.0, 1000.0)
    end

    # minimum medium
    # #TODO: find why iML1515 is not capable of growing whitout oxygen
    # This is the problem you found in the pan network from 
    # monkGenomescaleMetabolicNetwork2022
    bounds!(net,  "EX_glc__D_e", -10.0, 1000.0) # limiting nutrient
    bounds!(net,  "EX_o2_e", -1000.0, 1000.0) # essential, do not know why
    bounds!(net,  "EX_nh4_e", -1000.0, 1000.0) # essential

    bounds!(net,  "EX_ca2_e", -1000.0, 1000.0)
    bounds!(net,  "EX_cl_e", -1000.0, 1000.0)
    bounds!(net,  "EX_co2_e", -1000.0, 1000.0)
    bounds!(net,  "EX_cobalt2_e", -1000.0, 1000.0)
    bounds!(net,  "EX_cu2_e", -1000.0, 1000.0)
    bounds!(net,  "EX_fe2_e", -1000.0, 1000.0)
    bounds!(net,  "EX_fe3_e", -1000.0, 1000.0)
    bounds!(net,  "EX_h2o_e", -1000.0, 1000.0)
    bounds!(net,  "EX_h_e", -1000.0, 1000.0)
    bounds!(net,  "EX_k_e", -1000.0, 1000.0)
    bounds!(net,  "EX_mg2_e", -1000.0, 1000.0)
    bounds!(net,  "EX_mn2_e", -1000.0, 1000.0)
    bounds!(net,  "EX_mobd_e", -1000.0, 1000.0)
    bounds!(net,  "EX_na1_e", -1000.0, 1000.0)
    bounds!(net,  "EX_ni2_e", -1000.0, 1000.0)
    bounds!(net,  "EX_pi_e", -1000.0, 1000.0)
    bounds!(net,  "EX_so4_e", -1000.0, 1000.0)
    bounds!(net,  "EX_zn2_e", -1000.0, 1000.0)
    bounds!(net,  "EX_cu_e", -1000.0, 1000.0)
    

    # This is required
    # rxn[888]: EX_adocbl_e (Adenosylcobalamin exchange)
    #  subsys: nothing
    #  lb: -1000.0, ub: 1000.0
    #  (-1.0) adocbl_e <==> 
    bounds!(net,  "EX_adocbl_e", -1000.0, 1000.0)


    # extras
    extras!(net, "BIOM", "BIOMASS_Ec_iML1515_WT_75p37M")
    extras!(net, "EX_GLC", "EX_glc__D_e")
    extras!(net, "EX_NH4", "EX_nh4_e")
    extras!(net, "EX_GLU", "EX_glu__L_e")
    extras!(net, "EX_O2", "EX_o2_e")
    extras!(net, "EX_CO2", "EX_co2_e")
    extras!(net, "ATPM", "ATPM")
    
    linear_weights!(net, "BIOMASS_Ec_iML1515_WT_75p37M", 1.0)

    return net
end

function _register_iML1515()
    register_network!("iML1515", _iML1515_builder;
        use_cache = false,
        source = "http://bigg.ucsd.edu/models/iML1515", 
        desc = "A genome-scale model of E. coli, just as shipped"
    )
end