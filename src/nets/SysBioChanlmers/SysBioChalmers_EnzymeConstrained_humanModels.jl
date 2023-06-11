## ------------------------------------------------------------------
const __SysBioChalmers_EnzymeConstrained_humanModels_gemids = [
    "HOP62", "ec_HOP62",
    "HOP92", "ec_HOP92",
    "HS_578T", "ec_HS_578T",
    "HT29", "ec_HT29",
    "MALME_3M", "ec_MALME_3M",
    "MDMAMB_231", "ec_MDMAMB_231",
    "NCI_H226", "ec_NCI_H226",
    "O_786", "ec_O_786",
    "RPMI_8226", "ec_RPMI_8226",
    "SR", "ec_SR",
    "UO_31", "ec_UO_31",
]

## ------------------------------------------------------------------
function _SysBioChalmers_EnzymeConstrained_humanModels_tINIt_builder(gemid::String)
    # ------------------------------
    # cell line model
    
    net = _load_raw_model(
        "SysBioChalmers", "EnzymeConstrained_humanModels", "models", 
        gemid, 
        string(gemid, ".mat")
    )

    # fix size bug
    resize!(net.b, size(net, 1))
    
    # drug stuff gone
    empty_drug_reactions!(net)
    
    # Eliminate boundary metabolites
    for mi in eachindex(net.mets)
        met = net.mets[mi]
        endswith(met, 'x') || continue
        empty_met!(net, met; empty_void = false)
    end
    empty_void_iders!(net)
    @assert all(.!endswith.(net.mets, 'x'))

    # net
    net = emptyless_model(net)
    net = _common_format(net)
    
    # close exchanges
    for ri in eachindex(net.rxns)
        net.subSystems[ri] == "Exchange/demand reactions" || continue
        bounds!(net, ri, 0.0, 1000.0)
        # summary(net,net.rxns[ri])
    end

    # Open Ham medium
    # from: https://github.com/SysBioChalmers/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt GR	Growth on Ham's media (biomass production) task
    # medium_mets = ["MAM03135e","MAM01628e","MAM02982e","MAM02770e","MAM02159e","MAM01830e","MAM01401e","MAM02680e","MAM01513e","MAM02171e","MAM02583e","MAM02817e","MAM02842e","MAM02996e","MAM01361e","MAM02394e","MAM01965e","MAM02946e","MAM02387e","MAM02389e","MAM02630e","MAM02040e","MAM02833e","MAM02751e","MAM01821e","MAM01369e","MAM02184e","MAM02360e","MAM01327e","MAM01935e","MAM01365e","MAM02125e","MAM02426e","MAM02471e","MAM02724e","MAM03089e","MAM03101e","MAM01307e","MAM01986e","MAM02896e","MAM02993e","MAM01369e","MAM01974e","MAM01975e"]
    # medium_exchs = ["MAR09046", "MAR09065", "MAR09159", "MAR09068", "MAR09358", "MAR09146", "MAR09109", "MAR09145", "MAR09083", "MAR09361", "MAR09378", "MAR09144", "MAR09143", "MAR09423", "MAR09269", "MAR09167", "MAR09034", "MAR09074", "MAR09035", "MAR09036", "MAR09048", "MAR09047", "MAR09404", "MAR09072", "MAR09076", "MAR09062", "MAR09039", "MAR09040", "MAR09151", "MAR09153", "MAR09066", "MAR09038", "MAR09041", "MAR09042", "MAR09043", "MAR09045", "MAR09064", "MAR09061", "MAR09067", "MAR09069", "MAR09044", "MAR09062", "MAR09071", "MAR09063"]
    bounds!(net, "HMR_9034", -0.5, 1000.0) # glucose (metid: MAM01965e)
    bounds!(net, "HMR_9063", -0.5, 1000.0) # glutamine (metid: MAM01975e)
    bounds!(net, "HMR_9046", -0.1, 1000.0) # valine (metid: MAM03135e)
    bounds!(net, "HMR_9065", -0.1, 1000.0) # cysteine (metid: MAM01628e)
    bounds!(net, "HMR_9068", -0.1, 1000.0) # proline (metid: MAM02770e)
    bounds!(net, "HMR_9062", -0.1, 1000.0) # asparagine (metid: MAM01369e)
    bounds!(net, "HMR_9039", -0.1, 1000.0) # isoleucine (metid: MAM02184e)
    bounds!(net, "HMR_9040", -0.1, 1000.0) # leucine (metid: MAM02360e)
    bounds!(net, "HMR_9066", -0.1, 1000.0) # arginine (metid: MAM01365e)
    bounds!(net, "HMR_9038", -0.1, 1000.0) # histidine (metid: MAM02125e)
    bounds!(net, "HMR_9041", -0.1, 1000.0) # lysine (metid: MAM02426e)
    bounds!(net, "HMR_9042", -0.1, 1000.0) # methionine (metid: MAM02471e)
    bounds!(net, "HMR_9043", -0.1, 1000.0) # phenylalanine (metid: MAM02724e)
    bounds!(net, "HMR_9045", -0.1, 1000.0) # tryptophan (metid: MAM03089e)
    bounds!(net, "HMR_9064", -0.1, 1000.0) # tyrosine (metid: MAM03101e)
    bounds!(net, "HMR_9061", -0.1, 1000.0) # alanine (metid: MAM01307e)
    bounds!(net, "HMR_9067", -0.1, 1000.0) # glycine (metid: MAM01986e)
    bounds!(net, "HMR_9069", -0.1, 1000.0) # serine (metid: MAM02896e)
    bounds!(net, "HMR_9044", -0.1, 1000.0) # threonine (metid: MAM02993e)
    bounds!(net, "HMR_9062", -0.1, 1000.0) # asparagine (metid: MAM01369e)
    bounds!(net, "HMR_9071", -0.1, 1000.0) # glutamate (metid: MAM01974e)
    bounds!(net, "HMR_9109", -1000.0, 1000.0) # biotin (metid: MAM01401e)
    bounds!(net, "HMR_9145", -1000.0, 1000.0) # pantothenate (metid: MAM02680e)
    bounds!(net, "HMR_9167", -1000.0, 1000.0) # lipoic acid (metid: MAM02394e)
    bounds!(net, "HMR_9074", -1000.0, 1000.0) # sulfate (metid: MAM02946e)
    bounds!(net, "HMR_9404", -1000.0, 1000.0) # retinoate (metid: MAM02833e)
    bounds!(net, "HMR_9035", -1000.0, 1000.0) # linoleate (metid: MAM02387e)
    bounds!(net, "HMR_9036", -1000.0, 1000.0) # linolenate (metid: MAM02389e)
    bounds!(net, "HMR_9083", -1000.0, 1000.0) # choline (metid: MAM01513e)
    bounds!(net, "HMR_9361", -1000.0, 1000.0) # inositol (metid: MAM02171e)
    bounds!(net, "HMR_9143", -1000.0, 1000.0) # riboflavin (metid: MAM02842e)
    bounds!(net, "HMR_9378", -1000.0, 1000.0) # nicotinamide (metid: MAM02583e)
    bounds!(net, "HMR_9144", -1000.0, 1000.0) # pyridoxine (metid: MAM02817e)
    bounds!(net, "HMR_9423", -1000.0, 1000.0) # thymidine (metid: MAM02996e)
    bounds!(net, "HMR_9151", -1000.0, 1000.0) # alpha-tocopherol (metid: MAM01327e)
    bounds!(net, "HMR_9146", -1000.0, 1000.0) # folate (metid: MAM01830e)
    bounds!(net, "HMR_9358", -1000.0, 1000.0) # hypoxanthine (metid: MAM02159e)
    bounds!(net, "HMR_9269", -1000.0, 1000.0) # aquacob(III)alamin (metid: MAM01361e)
    bounds!(net, "HMR_9153", -1000.0, 1000.0) # gamma-tocopherol (metid: MAM01935e)
    bounds!(net, "HMR_9159", -1000.0, 1000.0) # thiamin (metid: MAM02982e)
    bounds!(net, "HMR_9048", -1000.0, 1000.0) # O2 (metid: MAM02630e)
    bounds!(net, "HMR_9047", -1000.0, 1000.0) # H2O (metid: MAM02040e)
    bounds!(net, "HMR_9072", -1000.0, 1000.0) # Pi (metid: MAM02751e)
    bounds!(net, "HMR_9076", -1000.0, 1000.0) # Fe2+ (metid: MAM01821e)

    # extras
    extras!(net, "BIOM", "biomass_human")
    extras!(net, "EX_GLC", "HMR_9034")
    extras!(net, "EX_GLN", "HMR_9063")
    extras!(net, "EX_NH4", "EX_nh4[e]")
    extras!(net, "EX_GLU", "HMR_9071")
    extras!(net, "EX_O2", "HMR_9048")
    extras!(net, "EX_CO2", "HMR_9058")
    extras!(net, "ATPM", "HMR_8474")

    return net
end

function _SysBioChalmers_EnzymeConstrained_humanModels_tINIT_EC_builder(gemid::String)

    ecnet = _load_raw_model(
        "SysBioChalmers", "EnzymeConstrained_humanModels", "models", 
        gemid, 
        "ecModel_batch.mat"
    )

    flx_inf = 1e6
    ecnet = _common_format(ecnet; flx_inf)
    resize!(ecnet.b, size(ecnet, 1))

    # drug stuff gone
    empty_drug_reactions!(ecnet)

    # Eliminate boundary metabolites
    for mi in eachindex(ecnet.mets)
        met = ecnet.mets[mi]
        endswith(met, 'x') || continue
        empty_met!(ecnet, met; empty_void = false)
    end
    empty_void_iders!(ecnet)
    @assert all(.!endswith.(ecnet.mets, 'x'))

    ecnet = emptyless_model(ecnet)

    # Exchanges: close input & open output
    for ri in eachindex(ecnet.rxns)
        ecnet.subSystems[ri] == "Exchange/demand reactions" || continue
        rxn = ecnet.rxns[ri]
        endswith(rxn, "_REV") && continue
        bounds!(ecnet, rxn, 0.0, flx_inf) # output
        bounds!(ecnet, string(rxn, "_REV"), 0.0, 0.0) # input
        # summary(ecnet, ecnet.rxns[ri])
    end

    # Open Ham medium
    # from: https://github.com/SysBioChalmers/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt GR	Growth on Ham's media (biomass production) task
    # medium_mets = ["MAM03135e","MAM01628e","MAM02982e","MAM02770e","MAM02159e","MAM01830e","MAM01401e","MAM02680e","MAM01513e","MAM02171e","MAM02583e","MAM02817e","MAM02842e","MAM02996e","MAM01361e","MAM02394e","MAM01965e","MAM02946e","MAM02387e","MAM02389e","MAM02630e","MAM02040e","MAM02833e","MAM02751e","MAM01821e","MAM01369e","MAM02184e","MAM02360e","MAM01327e","MAM01935e","MAM01365e","MAM02125e","MAM02426e","MAM02471e","MAM02724e","MAM03089e","MAM03101e","MAM01307e","MAM01986e","MAM02896e","MAM02993e","MAM01369e","MAM01974e","MAM01975e"]
    # medium_exchs = ["MAR09046", "MAR09065", "MAR09159", "MAR09068", "MAR09358", "MAR09146", "MAR09109", "MAR09145", "MAR09083", "MAR09361", "MAR09378", "MAR09144", "MAR09143", "MAR09423", "MAR09269", "MAR09167", "MAR09034", "MAR09074", "MAR09035", "MAR09036", "MAR09048", "MAR09047", "MAR09404", "MAR09072", "MAR09076", "MAR09062", "MAR09039", "MAR09040", "MAR09151", "MAR09153", "MAR09066", "MAR09038", "MAR09041", "MAR09042", "MAR09043", "MAR09045", "MAR09064", "MAR09061", "MAR09067", "MAR09069", "MAR09044", "MAR09062", "MAR09071", "MAR09063"]
    # medium_exchs = ["HMR_9034", "HMR_9063", "HMR_9046", "HMR_9065", "HMR_9068", "HMR_9062", "HMR_9039", "HMR_9040", "HMR_9066", "HMR_9038", "HMR_9041", "HMR_9042", "HMR_9043", "HMR_9045", "HMR_9064", "HMR_9061", "HMR_9067", "HMR_9069", "HMR_9044", "HMR_9062", "HMR_9071", "HMR_9109", "HMR_9145", "HMR_9167", "HMR_9074", "HMR_9404", "HMR_9035", "HMR_9036", "HMR_9083", "HMR_9361", "HMR_9143", "HMR_9378", "HMR_9144", "HMR_9423", "HMR_9151", "HMR_9146", "HMR_9358", "HMR_9269", "HMR_9153", "HMR_9159", "HMR_9048", "HMR_9047", "HMR_9072", "HMR_9076"]
    bounds!(ecnet, "HMR_9034_REV", 0.0, 0.5) # glucose (metid: MAM01965e)
    bounds!(ecnet, "HMR_9063_REV", 0.0, 0.5) # glutamine (metid: MAM01975e)
    bounds!(ecnet, "HMR_9046_REV", 0.0, 0.1) # valine (metid: MAM03135e)
    bounds!(ecnet, "HMR_9065_REV", 0.0, 0.1) # cysteine (metid: MAM01628e)
    bounds!(ecnet, "HMR_9068_REV", 0.0, 0.1) # proline (metid: MAM02770e)
    bounds!(ecnet, "HMR_9062_REV", 0.0, 0.1) # asparagine (metid: MAM01369e)
    bounds!(ecnet, "HMR_9039_REV", 0.0, 0.1) # isoleucine (metid: MAM02184e)
    bounds!(ecnet, "HMR_9040_REV", 0.0, 0.1) # leucine (metid: MAM02360e)
    bounds!(ecnet, "HMR_9066_REV", 0.0, 0.1) # arginine (metid: MAM01365e)
    bounds!(ecnet, "HMR_9038_REV", 0.0, 0.1) # histidine (metid: MAM02125e)
    bounds!(ecnet, "HMR_9041_REV", 0.0, 0.1) # lysine (metid: MAM02426e)
    bounds!(ecnet, "HMR_9042_REV", 0.0, 0.1) # methionine (metid: MAM02471e)
    bounds!(ecnet, "HMR_9043_REV", 0.0, 0.1) # phenylalanine (metid: MAM02724e)
    bounds!(ecnet, "HMR_9045_REV", 0.0, 0.1) # tryptophan (metid: MAM03089e)
    bounds!(ecnet, "HMR_9064_REV", 0.0, 0.1) # tyrosine (metid: MAM03101e)
    bounds!(ecnet, "HMR_9061_REV", 0.0, 0.1) # alanine (metid: MAM01307e)
    bounds!(ecnet, "HMR_9067_REV", 0.0, 0.1) # glycine (metid: MAM01986e)
    bounds!(ecnet, "HMR_9069_REV", 0.0, 0.1) # serine (metid: MAM02896e)
    bounds!(ecnet, "HMR_9044_REV", 0.0, 0.1) # threonine (metid: MAM02993e)
    bounds!(ecnet, "HMR_9062_REV", 0.0, 0.1) # asparagine (metid: MAM01369e)
    bounds!(ecnet, "HMR_9071_REV", 0.0, 0.1) # glutamate (metid: MAM01974e)
    bounds!(ecnet, "HMR_9109_REV", 0.0, flx_inf) # biotin (metid: MAM01401e)
    bounds!(ecnet, "HMR_9145_REV", 0.0, flx_inf) # pantothenate (metid: MAM02680e)
    bounds!(ecnet, "HMR_9167_REV", 0.0, flx_inf) # lipoic acid (metid: MAM02394e)
    bounds!(ecnet, "HMR_9074_REV", 0.0, flx_inf) # sulfate (metid: MAM02946e)
    bounds!(ecnet, "HMR_9404_REV", 0.0, flx_inf) # retinoate (metid: MAM02833e)
    bounds!(ecnet, "HMR_9035_REV", 0.0, flx_inf) # linoleate (metid: MAM02387e)
    bounds!(ecnet, "HMR_9036_REV", 0.0, flx_inf) # linolenate (metid: MAM02389e)
    bounds!(ecnet, "HMR_9083_REV", 0.0, flx_inf) # choline (metid: MAM01513e)
    bounds!(ecnet, "HMR_9361_REV", 0.0, flx_inf) # inositol (metid: MAM02171e)
    bounds!(ecnet, "HMR_9143_REV", 0.0, flx_inf) # riboflavin (metid: MAM02842e)
    bounds!(ecnet, "HMR_9378_REV", 0.0, flx_inf) # nicotinamide (metid: MAM02583e)
    bounds!(ecnet, "HMR_9144_REV", 0.0, flx_inf) # pyridoxine (metid: MAM02817e)
    bounds!(ecnet, "HMR_9423_REV", 0.0, flx_inf) # thymidine (metid: MAM02996e)
    bounds!(ecnet, "HMR_9151_REV", 0.0, flx_inf) # alpha0.0, tocopherol (metid: MAM01327e)
    bounds!(ecnet, "HMR_9146_REV", 0.0, flx_inf) # folate (metid: MAM01830e)
    bounds!(ecnet, "HMR_9358_REV", 0.0, flx_inf) # hypoxanthine (metid: MAM02159e)
    bounds!(ecnet, "HMR_9269_REV", 0.0, flx_inf) # aquacob(III)alamin (metid: MAM01361e)
    bounds!(ecnet, "HMR_9153_REV", 0.0, flx_inf) # gamma0.0, tocopherol (metid: MAM01935e)
    bounds!(ecnet, "HMR_9159_REV", 0.0, flx_inf) # thiamin (metid: MAM02982e)
    bounds!(ecnet, "HMR_9048_REV", 0.0, flx_inf) # O2 (metid: MAM02630e)
    bounds!(ecnet, "HMR_9047_REV", 0.0, flx_inf) # H2O (metid: MAM02040e)
    bounds!(ecnet, "HMR_9072_REV", 0.0, flx_inf) # Pi (metid: MAM02751e)
    bounds!(ecnet, "HMR_9076_REV", 0.0, flx_inf) # Fe2+ (metid: MAM01821e)
    
    # Protein pool (upper limit of the cost function)
    # C = 0.14825 # from 10.1126/scisignal.aaz1482
    # bounds!(ecnet, "prot_pool_exchange", 0.0, C) 
    bounds!(ecnet, "prot_pool_exchange", 0.0, flx_inf) # non restricted

    # extras
    extras!(ecnet, "BIOM", "biomass_human")
    extras!(ecnet, "EX_GLC", "HMR_9034_REV") # input
    extras!(ecnet, "EX_GLN", "HMR_9063_REV") # input
    extras!(ecnet, "EX_NH4", "EX_nh4[e]_REV") # input
    extras!(ecnet, "EX_GLU", "HMR_9071_REV")  # input
    extras!(ecnet, "EX_O2", "HMR_9048_REV") # input
    extras!(ecnet, "EX_CO2", "HMR_9058") # output
    extras!(ecnet, "ATPM", "arm_HMR_8474")
    extras!(ecnet, "PROT_POOL", "prot_pool_exchange")

    return ecnet
end

## ------------------------------------------------------------------
function _SysBioChalmers_EnzymeConstrained_humanModels_builder(gemid::String, biom_mod::String = "original")
    
    # check id
    ids = __SysBioChalmers_EnzymeConstrained_humanModels_gemids
    gemid in ids || error(
        "Unknown GEM id.\n", 
        " gemid: ", gemid, "\n", 
        " available: ", ids
    )

    if startswith(gemid, "ec_") # tINIT + EC models
        gemid = gemid[4:end] # remove 'ec_'
        net = _SysBioChalmers_EnzymeConstrained_humanModels_tINIT_EC_builder(gemid)
    else # tINIT models
        net = _SysBioChalmers_EnzymeConstrained_humanModels_tINIt_builder(gemid)
    end

    # biomass modification
    if biom_mod != "original"
        _SysBioChalmers_EnzymeConstrained_humanModels_Niklas_biomass!(net, biom_mod)
    end

    return net

end

function _register_SysBioChalmers_EnzymeConstrained_humanModels()
    register_network!("SysBioChalmers_EnzymeConstrained_humanModels", _SysBioChalmers_EnzymeConstrained_humanModels_builder;
        use_cache = true,
        source = "J. L. Robinson, P. Kocabas, H. Wang, P.-E. Cholley, et al. An atlas of human metabolism. Sci. Signal. 13, eaaz 1482 (2020). [doi:10.1126/scisignal.aaz1482](https://doi.org/10.1126/scisignal.aaz1482)", 
        download = """
            https://github.com/SysBioChalmers/EnzymeConstrained_humanModels, 
            commit a211fb69ba8816344bb06dce56d2a447df265550 (grafted, HEAD -> master, origin/master, origin/HEAD)
            Author: Iván Domenzain <ivand@chalmers.se>
            Date:   Mon Feb 3 17:40:11 2020 +0100
            
                Merge pull request #21 from SysBioChalmers/transfer2zenodo
                
                Transfer2zenodo
        """,
        desc = """
            The context-specific and enzymaticaly constraint models derived from Human-GEM with a few modifications (drug related reactions/metabolites removed) in Ham medium.
            Available lines: $(__SysBioChalmers_EnzymeConstrained_humanModels_gemids)
        """
    )
end