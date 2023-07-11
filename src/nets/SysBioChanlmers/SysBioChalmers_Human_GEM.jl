## ------------------------------------------------------------------
# Utils
# TODO: think about a Human_GEM.jl package for easy porting stuff from the MatLab original repo
function empty_drug_reactions!(net::MetNet; verbose = false)
    
    # from: https://github.com/SysBioChalmers/Human-GEM/blob/main/code/tINIT/removeDrugReactions.m
    drugNames = [
        "pravastatin","Gliclazide","atorvastatin","fluvastatin","fluvastain","fluvstatin",
        "simvastatin","cyclosporine","acetaminophen","cerivastatin","Tacrolimus","ibuprofen",
        "lovastatin","Losartan","nifedipine","pitavastatin","rosuvastatin","Torasemide","Midazolam", 
    ]

    drug_metis = Int[]
    for mi in eachindex(net.mets)
        metName = net.metNames[mi]
        for drugName in drugNames
            if contains(uppercase(metName), uppercase(drugName))
                verbose && println("-"^40)
                verbose && @show drugName metName
                push!(drug_metis, mi)
                break
            end
        end
    end
    

    drug_rxnis = Int[]
    for mi in drug_metis
        push!(drug_rxnis, met_rxns(net, mi)...)
    end

    # drugMets = contains(inModel.metNames,drugNames,'IgnoreCase',true);
    empty_met!(net, drug_metis; empty_void = false)
    empty_rxn!(net, drug_rxnis; empty_void = false)
    empty_void_iders!(net)
    
    return net

end

## ------------------------------------------------------------------
function _SysBioChanlmers_Human_GEM_builder()
    
    # load
    net = _load_raw_model("SysBioChalmers", "Human-GEM", "Human-GEM.mat")
    net = _common_format(net)

    # remove drug reactions
    net = empty_drug_reactions!(net)
    net = emptyless_model(net)

    # close intakes
    for i in eachindex(net.rxns)
        net.subSystems[i] == "Exchange/demand reactions" || continue
        bounds!(net, i, 0.0, 1000.0)
    end

    # open just Ham medium
    # from: https://github.com/SysBioChalmers/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt GR	Growth on Ham's media (biomass production) task
    # medium_mets = ["MAM03135e","MAM01628e","MAM02982e","MAM02770e","MAM02159e","MAM01830e","MAM01401e","MAM02680e","MAM01513e","MAM02171e","MAM02583e","MAM02817e","MAM02842e","MAM02996e","MAM01361e","MAM02394e","MAM01965e","MAM02946e","MAM02387e","MAM02389e","MAM02630e","MAM02040e","MAM02833e","MAM02751e","MAM01821e","MAM01369e","MAM02184e","MAM02360e","MAM01327e","MAM01935e","MAM01365e","MAM02125e","MAM02426e","MAM02471e","MAM02724e","MAM03089e","MAM03101e","MAM01307e","MAM01986e","MAM02896e","MAM02993e","MAM01369e","MAM01974e","MAM01975e"]
    # medium_exchs = ["MAR09046", "MAR09065", "MAR09159", "MAR09068", "MAR09358", "MAR09146", "MAR09109", "MAR09145", "MAR09083", "MAR09361", "MAR09378", "MAR09144", "MAR09143", "MAR09423", "MAR09269", "MAR09167", "MAR09034", "MAR09074", "MAR09035", "MAR09036", "MAR09048", "MAR09047", "MAR09404", "MAR09072", "MAR09076", "MAR09062", "MAR09039", "MAR09040", "MAR09151", "MAR09153", "MAR09066", "MAR09038", "MAR09041", "MAR09042", "MAR09043", "MAR09045", "MAR09064", "MAR09061", "MAR09067", "MAR09069", "MAR09044", "MAR09062", "MAR09071", "MAR09063"]
    # medium_exchs = ["HMR_9034", "HMR_9063", "HMR_9046", "HMR_9065", "HMR_9068", "HMR_9062", "HMR_9039", "HMR_9040", "HMR_9066", "HMR_9038", "HMR_9041", "HMR_9042", "HMR_9043", "HMR_9045", "HMR_9064", "HMR_9061", "HMR_9067", "HMR_9069", "HMR_9044", "HMR_9062", "HMR_9071", "HMR_9109", "HMR_9145", "HMR_9167", "HMR_9074", "HMR_9404", "HMR_9035", "HMR_9036", "HMR_9083", "HMR_9361", "HMR_9143", "HMR_9378", "HMR_9144", "HMR_9423", "HMR_9151", "HMR_9146", "HMR_9358", "HMR_9269", "HMR_9153", "HMR_9159", "HMR_9048", "HMR_9047", "HMR_9072", "HMR_9076"]
    bounds!(net, "MAR09034", -0.5, 1000.0) # glucose (metid: MAM01965e)
    bounds!(net, "MAR09063", -0.5, 1000.0) # glutamine (metid: MAM01975e)
    bounds!(net, "MAR09046", -0.1, 1000.0) # valine (metid: MAM03135e)
    bounds!(net, "MAR09065", -0.1, 1000.0) # cysteine (metid: MAM01628e)
    bounds!(net, "MAR09068", -0.1, 1000.0) # proline (metid: MAM02770e)
    bounds!(net, "MAR09062", -0.1, 1000.0) # asparagine (metid: MAM01369e)
    bounds!(net, "MAR09039", -0.1, 1000.0) # isoleucine (metid: MAM02184e)
    bounds!(net, "MAR09040", -0.1, 1000.0) # leucine (metid: MAM02360e)
    bounds!(net, "MAR09066", -0.1, 1000.0) # arginine (metid: MAM01365e)
    bounds!(net, "MAR09038", -0.1, 1000.0) # histidine (metid: MAM02125e)
    bounds!(net, "MAR09041", -0.1, 1000.0) # lysine (metid: MAM02426e)
    bounds!(net, "MAR09042", -0.1, 1000.0) # methionine (metid: MAM02471e)
    bounds!(net, "MAR09043", -0.1, 1000.0) # phenylalanine (metid: MAM02724e)
    bounds!(net, "MAR09045", -0.1, 1000.0) # tryptophan (metid: MAM03089e)
    bounds!(net, "MAR09064", -0.1, 1000.0) # tyrosine (metid: MAM03101e)
    bounds!(net, "MAR09061", -0.1, 1000.0) # alanine (metid: MAM01307e)
    bounds!(net, "MAR09067", -0.1, 1000.0) # glycine (metid: MAM01986e)
    bounds!(net, "MAR09069", -0.1, 1000.0) # serine (metid: MAM02896e)
    bounds!(net, "MAR09044", -0.1, 1000.0) # threonine (metid: MAM02993e)
    bounds!(net, "MAR09070", -0.1, 1000.0) # aspartate (metid: )
    bounds!(net, "MAR09071", -0.1, 1000.0) # glutamate (metid: MAM01974e)
    bounds!(net, "MAR09109", -1000.0, 1000.0) # biotin (metid: MAM01401e)
    bounds!(net, "MAR09145", -1000.0, 1000.0) # pantothenate (metid: MAM02680e)
    bounds!(net, "MAR09167", -1000.0, 1000.0) # lipoic acid (metid: MAM02394e)
    bounds!(net, "MAR09074", -1000.0, 1000.0) # sulfate (metid: MAM02946e)
    bounds!(net, "MAR09404", -1000.0, 1000.0) # retinoate (metid: MAM02833e)
    bounds!(net, "MAR09035", -1000.0, 1000.0) # linoleate (metid: MAM02387e)
    bounds!(net, "MAR09036", -1000.0, 1000.0) # linolenate (metid: MAM02389e)
    bounds!(net, "MAR09083", -1000.0, 1000.0) # choline (metid: MAM01513e)
    bounds!(net, "MAR09361", -1000.0, 1000.0) # inositol (metid: MAM02171e)
    bounds!(net, "MAR09143", -1000.0, 1000.0) # riboflavin (metid: MAM02842e)
    bounds!(net, "MAR09378", -1000.0, 1000.0) # nicotinamide (metid: MAM02583e)
    bounds!(net, "MAR09144", -1000.0, 1000.0) # pyridoxine (metid: MAM02817e)
    bounds!(net, "MAR09423", -1000.0, 1000.0) # thymidine (metid: MAM02996e)
    bounds!(net, "MAR09151", -1000.0, 1000.0) # alpha-tocopherol (metid: MAM01327e)
    bounds!(net, "MAR09146", -1000.0, 1000.0) # folate (metid: MAM01830e)
    bounds!(net, "MAR09358", -1000.0, 1000.0) # hypoxanthine (metid: MAM02159e)
    bounds!(net, "MAR09269", -1000.0, 1000.0) # aquacob(III)alamin (metid: MAM01361e)
    bounds!(net, "MAR09153", -1000.0, 1000.0) # gamma-tocopherol (metid: MAM01935e)
    bounds!(net, "MAR09159", -1000.0, 1000.0) # thiamin (metid: MAM02982e)
    bounds!(net, "MAR09048", -1000.0, 1000.0) # O2 (metid: MAM02630e)
    bounds!(net, "MAR09047", -1000.0, 1000.0) # H2O (metid: MAM02040e)
    bounds!(net, "MAR09072", -1000.0, 1000.0) # Pi (metid: MAM02751e)
    bounds!(net, "MAR09076", -1000.0, 1000.0) # Fe2+ (metid: MAM01821e)

    extras!(net, "BIOM", "MAR13082")
    extras!(net, "EX_GLC", "MAR09034")
    extras!(net, "EX_GLN", "MAR09063")
    extras!(net, "EX_NH4", "MAR11420")
    extras!(net, "EX_GLU", "MAR09071")
    extras!(net, "EX_O2", "MAR09048")
    extras!(net, "EX_CO2", "MAR09058")
    extras!(net, "ATPM", "MAR03964")
    
    return net

end

function _register_SysBioChanlmers_Human_GEM()
    register_network!("SysBioChalmers_Human_GEM", _SysBioChanlmers_Human_GEM_builder;
        use_cache = true,
        source = "J. L. Robinson, P. Kocabas, H. Wang, P.-E. Cholley, et al. An atlas of human metabolism. Sci. Signal. 13, eaaz 1482 (2020). [doi:10.1126/scisignal.aaz1482](https://doi.org/10.1126/scisignal.aaz1482)", 
        download = """https://github.com/SysBioChalmers/Human-GEM, 
            commit 194ebe5431c83e25f78df61caacad2fa485b5cb4 (grafted, HEAD -> main, tag: v1.15.0, origin/main, origin/HEAD)
            Author: Hao Wang <haowang.bioinfo@gmail.com>
            Date:   Tue May 9 10:07:23 2023 +0200
        """,
        desc = "The Human-GEM general model with a few modifications (drug related reactions/metabolites removed) in Ham medium"
    )
end