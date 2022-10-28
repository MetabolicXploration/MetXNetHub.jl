function _ECCcomp_builder()

    # TODO: make it growth

    # load
    # 41598_2017_BFsrep39647_MOESM453_ESM
    net = _load_raw_model("ECC2comp.xml")

    # elimiate external mets
    ex_mets = filter(metabolites(net)) do met
        endswith(met, "_ex")
    end
    empty_met!(net, ex_mets)  
    # and Biomass met
    empty_met!(net, "Biomass")  
    net = emptyless_model(net)
    
    # open out close in
    # This net hase both reaction Ex and Up separated
    ex_rxns = filter(reactions(net)) do rxn
        endswith(rxn, "Up")
        endswith(rxn, "Ex")
    end
    bounds!(net, ex_rxns; lb = 0.0, ub = 1000.0)

    # minimum medium (open in)
    
    set_extra!(net, "BIOM", "Growth")
    set_extra!(net, "EX_GLC", "R_SuccUp")

    lin_objective!(net, "Growth", 1.0)

    return net
end

function _register_ECC2comp()
    register_network!("ECC2comp", _ECCcomp_builder;
        source = "Hädicke, O., Klamt, S. EColiCore2: a reference network model of the central metabolism of Escherichia coli and relationships to its genome-scale parent model. Sci Rep 7, 39647 (2017). https://doi.org/10.1038/srep39647", 
        desc = "The compressed EColiCore2 model (see source)  in a minimum glucose rich medium (max uptake: 10.0 mmol/gDW h), the raw model was extracted from 41598_2017_BFsrep39647_MOESM453_ESM.xls supplementary data"
    )
end