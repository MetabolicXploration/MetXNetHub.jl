## ------------------------------------------------------------------
const __SysBioChalmers_EnzymeConstrained_humanModels_gemids = ["HOP62", "HOP92", "HS_578T", "HT29", "MALME_3M", "MDMAMB_231", "NCI_H226", "O_786", "RPMI_8226", "SR", "UO_31"]
## ------------------------------------------------------------------
function _SysBioChalmers_EnzymeConstrained_humanModels_contextualized_builder(gemid::String)
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

    # @show size(net)

    # # test fba
    # biom_id = "biomass_human"
    # println("-"^50)
    # println("As shiped (no boundary)")
    # println()
    # opm = FBAOpModel(net, Gurobi.Optimizer)
    # optimize!(opm)
    # @show solution(opm, biom_id)
    
    # close exchanges
    for ri in eachindex(net.rxns)
        net.subSystems[ri] == "Exchange/demand reactions" || continue
        bounds!(net, ri, 0.0, 1000.0)
        # summary(net,net.rxns[ri])
    end

    # println("-"^50)
    # println("Intake closed")
    # println()
    # opm = FBAOpModel(net, Gurobi.Optimizer)
    # optimize!(opm)
    # @show solution(opm, biom_id)

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

    # println("-"^50)
    # println("Ham open")
    # println()
    # opm = FBAOpModel(net, Gurobi.Optimizer)
    # optimize!(opm)
    # @show solution(opm, biom_id)

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

function _SysBioChalmers_EnzymeConstrained_humanModels_contextualized_EP_builder(gemid::String)
    gemid = gemid[4:end] # remove 'ec_'
    
    # WIP: Move stuff from dev to here
end

## ------------------------------------------------------------------
function _SysBioChalmers_EnzymeConstrained_humanModels_builder(gemid::String)
    ids = __SysBioChalmers_EnzymeConstrained_humanModels_gemids
    
    if startswith(gemid, "ec_")
        ec_ids = string.(["ec_"], ids)
        gemid in ec_ids || error("Unknown GEM id. Available ", ec_ids)
        
        _SysBioChalmers_EnzymeConstrained_humanModels_contextualized_EP_builder(gemid)
    else
        gemid in ids || error("Unknown GEM id. Available ", ec_ids)
        _SysBioChalmers_EnzymeConstrained_humanModels_contextualized_builder(gemid)
    end

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