# read Beg et al., 2007 https://doi.org/10.1073/pnas.0609845104
# for better underestanding

function _ecoli_core_Beg2007_builder()

    # load
    net = pull_net("ecoli_core")

    # -------------------------------------------
    # Add nutrient uptake reations
    # NOTE: reactions are taken from iJO1366 with the help of Beg et al. 2007 fig 3...
    net = resize(net; 
        nmets = size(net, 1) + 11,
        nrxns = size(net, 2) + 13,
    )

    # -------------------------------------------
    # Galactose

    # rxn[935]: EX_gal_e (D-Galactose exchange)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) gal_e ==> 
    rxnid = "EX_gal_e"
    @assert !hascolid(net, rxnid)
    set_constraint!(net, rxnid;
        S = Dict("gal_e" => -1),
        lb = 0.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "EX_gal_e")
    net.rxnNames[rxnidx] = "D-Galactose exchange"
    net.subSystems[rxnidx] = "Extracellular exchange"
    metidx = metindex(net, "gal_e")
    net.metNames[metidx] = "D-Galactose"
    # println(rxnid, ": ", col_str(net, rxnid))

    # rxn[1310]: GALtex (D-galactose transport via diffusion (extracellular to periplasm))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) gal_e <==> (1.0) gal_p
    # Not apply

    # rxn[1308]: GALabcpp (D-galactose transport via ABC system (periplasm))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) atp_c + (-1.0) gal_p + (-1.0) h2o_c ==> (1.0) adp_c + (1.0) gal_c + (1.0) h_c + (1.0) pi_c
    rxnid = "GALabcpp"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "gal_e")
    @assert hasrowid(net, "atp_c")
    @assert hasrowid(net, "h2o_e")
    @assert hasrowid(net, "adp_c")
    @assert !hasrowid(net, "gal_c")
    @assert hasrowid(net, "h_c")
    @assert hasrowid(net, "pi_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "atp_c" => -1, "gal_e" => -1, "h2o_e" => -1, 
            "adp_c" => 1, "gal_c" => 1, "h_c" => 1, "pi_c" => 1, 
        ),
        lb = 0.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "GALabcpp")
    net.rxnNames[rxnidx] = "D-galactose transport via ABC system (periplasm)"
    net.subSystems[rxnidx] = "Transport, Inner Membrane"
    metidx = metindex(net, "gal_c")
    net.metNames[metidx] = "D-Galactose"
    # println(rxnid, ": ", col_str(net, rxnid))

    # rxn[1299]: GALKr (Galactokinase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) atp_c + (-1.0) gal_c <==> (1.0) adp_c + (1.0) gal1p_c + (1.0) h_c
    rxnid = "GALKr"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "gal_c")
    @assert !hasrowid(net, "gal1p_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "atp_c" => -1, "gal_c" => -1, 
            "adp_c" => 1, "gal1p_c" => 1, "h_c" => 1
        ),
        lb = -1000.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "GALKr")
    net.rxnNames[rxnidx] = "Galactokinase"
    net.subSystems[rxnidx] = "Alternate Carbon Metabolism"
    metidx = metindex(net, "gal1p_c")
    net.metNames[metidx] = "Alpha-D-Galactose 1-phosphate"
    # println(rxnid, ": ", col_str(net, rxnid))

    # rxn[2522]: UGLT (UDPglucose--hexose-1-phosphate uridylyltransferase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) gal1p_c + (-1.0) udpg_c <==> (1.0) g1p_c + (1.0) udpgal_c
    rxnid = "UGLT"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "gal1p_c")
    @assert !hasrowid(net, "udpg_c")
    @assert !hasrowid(net, "g1p_c")
    @assert !hasrowid(net, "udpgal_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "gal1p_c" => -1, "udpg_c" => -1, 
            "g1p_c" => 1, "udpgal_c" => 1
        ),
        lb = -1000.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "UGLT")
    net.rxnNames[rxnidx] = "UDPglucose--hexose-1-phosphate uridylyltransferase"
    net.subSystems[rxnidx] = "Alternate Carbon Metabolism"
    metidx = metindex(net, "udpg_c")
    net.metNames[metidx] = "UDPglucose"
    metidx = metindex(net, "g1p_c")
    net.metNames[metidx] = "D-Glucose 1-phosphate"
    metidx = metindex(net, "udpgal_c")
    net.metNames[metidx] = "UDPgalactose"
    # println(rxnid, ": ", col_str(net, rxnid))

    # rxn[2511]: UDPG4E (UDPglucose 4-epimerase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) udpg_c <==> (1.0) udpgal_c
    rxnid = "UDPG4E"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "udpg_c")
    @assert hasrowid(net, "udpgal_c")
    @assert hasrowid(net, "g1p_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "udpg_c" => -1, 
            "udpgal_c" => 1,
        ),
        lb = -1000.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "UDPG4E")
    net.rxnNames[rxnidx] = "UDPglucose 4-epimerase"
    net.subSystems[rxnidx] = "Alternate Carbon Metabolism"
    # println(rxnid, ": ", col_str(net, rxnid))

    # rxn[2082]: PGMT (Phosphoglucomutase)
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) g1p_c <==> (1.0) g6p_c
    rxnid = "PGMT"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "g6p_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "g1p_c" => -1, 
            "g6p_c" => 1,
        ),
        lb = -1000.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "PGMT")
    net.rxnNames[rxnidx] = "Phosphoglucomutase"
    net.subSystems[rxnidx] = "Alternate Carbon Metabolism"
    # println(rxnid, ": ", col_str(net, rxnid))

    # -------------------------------------------
    # Maltose

    # println()
    # println("-"^40)
    # println("Maltose")
    # println("-"^40)

    # rxn[1000]: EX_malt_e (Maltose exchange)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) malt_e ==> 
    rxnid = "EX_malt_e"
    @assert !hascolid(net, rxnid)
    @assert !hasrowid(net, "malt_e")
    set_constraint!(net, rxnid;
        S = Dict(
            "malt_e" => -1, 
        ),
        lb = -0.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "EX_malt_e")
    net.rxnNames[rxnidx] = "Maltose exchange"
    net.subSystems[rxnidx] = "Extracellular exchange"
    metidx = metindex(net, "malt_e")
    net.metNames[metidx] = "Maltose C12H22O11"
    # println(rxnid, ": ", col_str(net, rxnid))

    # rxn[1726]: MALTtexi (MaltoseMaltotriose transport via diffusion (extracellular to periplasm) irreversible)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) malt_e ==> (1.0) malt_p
    # Not apply

    # malEGFK
    # rxn[1724]: MALTabcpp (Maltose transport via ABC system (periplasm))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) atp_c + (-1.0) h2o_c + (-1.0) malt_p ==> (1.0) adp_c + (1.0) h_c + (1.0) malt_c + (1.0) pi_c
    rxnid = "MALTabcpp"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "malt_e")
    @assert !hasrowid(net, "malt_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "atp_c" => -1, "h2o_c" => -1, "malt_e" => -1, 
            "adp_c" => 1, "h_c" => 1, "malt_c" => 1, "pi_c" => 1
        ),
        lb = -0.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "MALTabcpp")
    net.rxnNames[rxnidx] = "Maltose transport via ABC system (periplasm)"
    net.subSystems[rxnidx] = "Transport, Inner Membrane"
    metidx = metindex(net, "malt_c")
    net.metNames[metidx] = "Maltose C12H22O11"
    # println(rxnid, ": ", col_str(net, rxnid))

    # NOTE
    # rxn[1725]: MALTptspp (Maltose transport via PEP:Pyr PTS (periplasm))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) malt_p + (-1.0) pep_c ==> (1.0) malt6p_c + (1.0) pyr_c
    # NOTE: malt6p_c is a dead end at iJO1366, so, using MALTabcpp

    # rxn[318]: AMALT1 (Amylomaltase (maltotriose))
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) malt_c + (-1.0) malttr_c ==> (1.0) glc__D_c + (1.0) maltttr_c
    # I will ignore maltttr_c
    rxnid = "AMALT1"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "glc__D_e")
    @assert hasrowid(net, "malt_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "malt_c" => -1,
            "glc__D_e" => 2,
        ),
        lb = -0.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "AMALT1")
    net.rxnNames[rxnidx] = "Amylomaltase (maltotriose)"
    net.subSystems[rxnidx] = "Alternate Carbon Metabolism"
    # println(rxnid, ": ", col_str(net, rxnid))

    # -------------------------------------------
    # Glycerol

    # println()
    # println("-"^40)
    # println("Glycerol")
    # println("-"^40)

    # rxn[958]: EX_glyc_e (Glycerol exchange)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) glyc_e ==> 
    rxnid = "EX_glyc_e"
    @assert !hascolid(net, rxnid)
    @assert !hasrowid(net, "glyc_e")
    set_constraint!(net, rxnid;
        S = Dict(
            "glyc_e" => -1,
        ),
        lb = -0.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "EX_glyc_e")
    net.rxnNames[rxnidx] = "Glycerol exchange"
    net.subSystems[rxnidx] = "Extracellular exchange"
    metidx = metindex(net, "glyc_e")
    net.metNames[metidx] = "Glycerol"
    # println(rxnid, ": ", col_str(net, rxnid))

    # rxn[1406]: GLYCtex (Glycerol transport via diffusion (extracellular to periplasm))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) glyc_e <==> (1.0) glyc_p
    # Not apply

    # rxn[1407]: GLYCtpp (Glycerol transport via channel (periplasm))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) glyc_c <==> (1.0) glyc_p
    rxnid = "GLYCtpp"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "glyc_e")
    @assert !hasrowid(net, "glyc_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "glyc_e" => -1,
            "glyc_c" => 1,
        ),
        lb = -1000.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "GLYCtpp")
    net.rxnNames[rxnidx] = "Glycerol transport via channel (periplasm)"
    net.subSystems[rxnidx] = "Transport, Inner Membrane"
    metidx = metindex(net, "glyc_c")
    net.metNames[metidx] = "Glycerol"
    # println(rxnid, ": ", col_str(net, rxnid))

    # rxn[1408]: GLYK (Glycerol kinase)
    # subsys: nothing
    # lb: 0.0, ub: 1000.0
    # (-1.0) atp_c + (-1.0) glyc_c ==> (1.0) adp_c + (1.0) glyc3p_c + (1.0) h_c
    rxnid = "GLYK"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "glyc_c")
    @assert !hasrowid(net, "glyc3p_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "atp_c" => -1, "glyc_c" => -1,
            "adp_c" => 1, "glyc3p_c" => 1, "h_c" => 1
        ),
        lb = 0.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "GLYK")
    net.rxnNames[rxnidx] = "Glycerol kinase"
    net.subSystems[rxnidx] = "Alternate Carbon Metabolism"
    metidx = metindex(net, "glyc3p_c")
    net.metNames[metidx] = "Glycerol 3-phosphate"
    # println(rxnid, ": ", col_str(net, rxnid))

    # rxn[1267]: G3PD2 (Glycerol-3-phosphate dehydrogenase (NADP))
    # subsys: nothing
    # lb: -1000.0, ub: 1000.0
    # (-1.0) glyc3p_c + (-1.0) nadp_c <==> (1.0) dhap_c + (1.0) h_c + (1.0) nadph_c
    rxnid = "G3PD2"
    @assert !hascolid(net, rxnid)
    @assert hasrowid(net, "glyc3p_c")
    @assert hasrowid(net, "dhap_c")
    @assert hasrowid(net, "nadp_c")
    @assert hasrowid(net, "nadph_c")
    set_constraint!(net, rxnid;
        S = Dict(
            "glyc3p_c" => -1, "nadp_c" => -1,
            "dhap_c" => 1, "h_c" => 1, "nadph_c" => 1
        ),
        lb = -1000.0, ub = 1000.0,
    )
    # meta
    rxnidx = rxnindex(net, "G3PD2")
    net.rxnNames[rxnidx] = "Glycerol-3-phosphate dehydrogenase (NADP)"
    net.subSystems[rxnidx] = "Alternate Carbon Metabolism"
    # println(rxnid, ": ", col_str(net, rxnid))

    # -------------------------------------------
    # trim
    net = emptyless_model(net)

    # -------------------------------------------
    # open complex medium

    # Bounds from Beg et al, 2007, fig 3
    # We just need the maximum rates. 
    # Original in (mmol/ min g)
    # 1 [mmol/ min g] * 60 = 1 [mmol/ h g]
    # adjustment growth factor (I aim to adjust the maximum growth to match glucose only regime ~0.8 1/h (Fig2, a))
    # At least I keep the relation between maximal nutrient intakes
    gf = 0.18
    lb!(net, "EX_glc__D_e", -0.9 * 60 * gf)
    lb!(net, "EX_lac__D_e", -1.0 * 60 * gf)
    lb!(net, "EX_malt_e", -0.1 * 60 * gf)
    lb!(net, "EX_gal_e", -0.2 * 60 * gf)
    lb!(net, "EX_glyc_e", -0.6 * 60 * gf)
    lb!(net, "EX_ac_e", -1.5 * 60 * gf)
    ub!(net, "EX_ac_e", 2.0 * 60 * gf)
    
    # extras
    extras!(net, "BIOM", "BIOMASS_Ecoli_core_w_GAM")
    extras!(net, "EX_GLC", "EX_glc__D_e")
    extras!(net, "EX_LAC", "EX_lac__D_e")
    extras!(net, "EX_MALT", "EX_malt_e")
    extras!(net, "EX_GAL", "EX_gal_e")
    extras!(net, "EX_GLYC", "EX_glyc_e")
    extras!(net, "EX_AC", "EX_ac_e")
    extras!(net, "EX_NH4", "EX_nh4_e")
    extras!(net, "EX_GLU", "EX_glu__L_e")
    extras!(net, "EX_O2", "EX_o2_e")
    extras!(net, "EX_CO2", "EX_co2_e")
    extras!(net, "ATPM", "ATPM")

    linear_weights!(net, "BIOMASS_Ecoli_core_w_GAM", 1.0)

    return net
end

function _register_ecoli_core_Beg2007()
    register_network!("ecoli_core_Beg2007", _ecoli_core_Beg2007_builder;
        use_cache = false,
        source = "Beg, Q. K., A. Vazquez, J. Ernst, M. A. de Menezes, Z. Bar-Joseph, A.-L. Barabási, and Z. N. Oltvai. 'Intracellular Crowding Defines the Mode and Sequence of Substrate Uptake by Escherichia Coli and Constrains Its Metabolic Activity.' Proceedings of the National Academy of Sciences 104, no. 31 (July 31, 2007): 12663-68. https://doi.org/10.1073/pnas.0609845104.", 
        download = "",
        desc = "The clasical core model of E. coli ('ecoli_core') with extra reactions to model Beg et al. 2007 nutrient intake patterns. The intake bounds try to match the max growth rate reported at the paper for each culture phase."
    )
end