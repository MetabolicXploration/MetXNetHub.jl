using JSON
using CSV
using DataFrames

function __PAMModel_alterProteomeRegulationPatterns2021_builder(
        options::Dict = Dict()
    )
    
    ## - --- -- . .- -.-. -.- .-.-. 
    # extract options
    f_UE0 = get(options, "f_UE0", 0.17)                 # [g / gCDW]
    u_ider = get(options, "u_ider", "EX_glc__D_e_b")    # nutrinet ider
    u_max = get(options, "u_max", 8.9)                  # [mmol/ h gCDW] max 
    f_T0 = get(options, "f_T0", 50.0)                   # [mg / gCDW]
    w_T = get(options, "w_T", 36.8)                     # [mg h / gCDW]
    f_Pc = get(options, "f_Pc", 0.26)                   # [g / gCDW]
    
    ## - --- -- . .- -.-. -.- .-.-. 
    # MARK: Load iML1515
    net0 = pull_net("iML1515")
    
    # Extract
    S0, b0, lb0, ub0, mets0, rxns0 = net0.S, net0.b, net0.lb, net0.ub, net0.mets, net0.rxns

    ## - --- -- . .- -.-. -.- .-.-. 
    # MARK: Load PaperSON
    # Load data from alterProteomeRegulationPatterns2021:
    # "Alter, Tobias B., Lars M. Blank, and Birgitta E. Ebert. “Proteome Regulation Patterns Determine Escherichia Coli Wild-Type and Mutant Phenotypes.” Edited by Joshua E. Elias. mSystems 6, no. 2 (April 27, 2021): 10.1128/msystems.00625-20. https://doi.org/10.1128/msystems.00625-20."
    paperSONDataDir = joinpath(pkgdir(MetXNetHub), "data/PaperSON")
    if (isdir(paperSONDataDir) == false)
        error("data/PaperSON directory not found. See https://github.com/MetabolicXploration/PaperSON")
    end
    alterDataDir = joinpath(paperSONDataDir, "data/alterProteomeRegulationPatterns2021")

    json_file = joinpath(alterDataDir, "raw.json")
    rawData = JSON.parsefile(json_file)
    file = rawData["data"]["msystems.00625-20-sd001.xlsx"]["PaperSON.file"]
    json_file = joinpath(alterDataDir, file)
    sd001Data = JSON.parsefile(json_file)

    file = sd001Data["data"]["Sheet 1"]["data"]["PaperSON.file"]
    csv_file = joinpath(alterDataDir, file)
    UEData = CSV.read(csv_file, DataFrame)

    ## - --- -- . .- -.-. -.- .-.-. 
    # MARK: split
    # We need to split all reactions because 
    # all variable `x` should have `wx = 0` constraints
    # Split reactions
    # S_f v_f + S_b v_b = 0
    # S_f = S0
    # S_b = -S0
    # v = v_f - v_b

    rxns1 = [
        [string(rxn, "_f") for rxn in rxns0]; 
        [string(rxn, "_b") for rxn in rxns0]; 
    ]
    mets1 = mets0
    S1 = [
        S0  -S0
    ]
    M, N = size(S1)
    b1 = b0
    c1 = zeros(N)
    lb1 = zeros(N)
    ub1 = [ub0; abs.(lb0)]

    # toSplit = Int[]
    # for rxn in UEData[!,"Reaction ID"]
    #     rxn = replace(rxn, r"_f$" => "")
    #     rxn = replace(rxn, r"_b$" => "")
    #     rxni = rxnindex(net0, rxn)
    #     push!(toSplit, rxni)
    # end
    # sort!(toSplit)
    # tsN = length(toSplit)
    # N = size(S0, 2)
    # notToSplit = setdiff(1:N, toSplit)

    # rxns1 = [
    #     [string(rxns0[rxni], "_f") for rxni in toSplit]; 
    #     [string(rxns0[rxni], "_b") for rxni in toSplit]; 
    #     [rxns0[rxni] for rxni in notToSplit]; 
    # ]
    # mets1 = mets0
    # S1 = [
    #     S0[:, toSplit]  -S0[:, toSplit]  S0[:, notToSplit]
    # ]
    # M, N = size(S1)
    # b1 = b0
    # c1 = zeros(N)
    # lb1 = [zeros(tsN); zeros(tsN); lb0[notToSplit]]
    # ub1 = [ub0[toSplit]; abs.(lb0[toSplit]); ub0[notToSplit]]

    ## - --- -- . .- -.-. -.- .-.-. 
    # MARK: f_AE
    # Add active enzyme contraints
    # An stash variable 'f_AE_stash_rxn' represents the
    # active enzymes proteome sector
    # 
    # p_e = v_e * M_e / kcat_e
    # ec_e = M_e / kcat_e    # [g s / mol]

    # enzyme cost coes
    # ec_e * v_e - p_e_stash = 0
    ec_dict = Dict() # enzymatic cost coe
    for (rxni, rxn) in enumerate(UEData[!,"Reaction ID"])
        M_e = UEData[rxni, "Enzyme molar mass [g/mol]"]
        M_e = M_e * 1e-3 # [g/mmol]
        kcat_e = UEData[rxni, "Turnover number (kcat) [1/s]"]
        kcat_e = kcat_e * 60 * 60 # [1 / h]
        
        if endswith(rxn, "_b")
            get!(ec_dict, rxn, 0)
            ec_dict[rxn] += M_e / kcat_e 
        elseif endswith(rxn, "_f")
            get!(ec_dict, rxn, 0)
            ec_dict[rxn] += M_e / kcat_e 
        else
            rxnf = string(rxn, "_f")
            get!(ec_dict, rxnf, 0)
            ec_dict[rxnf] += M_e / kcat_e 
            rxnb = string(rxn, "_b")
            get!(ec_dict, rxnb, 0)
            ec_dict[rxnb] += M_e / kcat_e 
        end
    end

    N = size(S1, 2)
    ec_coes = zeros(N)
    for (rxni, rxn) in enumerate(rxns1)
        ec_coes[rxni] = get(ec_dict, rxn, 0)
    end


    # f_AE stash variable
    rxns2 = [rxns1; "f_AE_stash_rxn"]
    mets2 = [mets1; "AE_stash_met"]
    lb2 = [lb1; 0]
    ub2 = [ub1; 1000.0]
    c2 = [c1; 0]
    b2 = [b1; 0]
    _fill0 = zeros(size(S1, 1))
    S2 = [
        S1          _fill0
        ec_coes'   -1
    ]

    ## - --- -- . .- -.-. -.- .-.-. 
    # MARK: f_UE
    # new balance
    # f_UE0 - w_UE * v_s - f_UE  = 0
    w_UE_coes = zeros(size(S2, 2))

    # ISSUE/QUESTION: 
    # - can this be extended to multiple nutrients at the same time?
    # - Maybe having a calibration curve for tot carbon
    Sidx = findfirst(rxns2 .== u_ider)
    # nutrient intake
    @assert u_max > 0
    w_UE_coes[Sidx] = f_UE0 / u_max            # [h / g mmol] 
    
    # f_AE stash variable
    rxns3 = [rxns2; "f_UE0_stash_rxn"; "f_UE_stash_rxn"]
    mets3 = [mets2; "UE_stash_met"]
    lb3 = [lb2; f_UE0; 0]
    ub3 = [ub2; f_UE0; 1000.0]
    c3 = [c2; 0; 0]
    b3 = [b2; 0]
    _col0 = zeros(size(S2, 1))
    S3 = [
        S2          _col0       _col0
        -w_UE_coes'  1           -1
    ]

    ## - --- -- . .- -.-. -.- .-.-. 
    # MARK: f_T
    # f_T0 + w_T z - f_T = 0
    f_T0 = f_T0 * 1e-3              # [g / gCDW]
    w_T = w_T * 1e-3                # [g / gCDW]

    zider = "BIOMASS_Ec_iML1515_WT_75p37M_f"    # biomass id
    zidx = findfirst(rxns3 .== zider)
    w_T_coes = zeros(size(S3, 2))
    w_T_coes[zidx] = w_T                        # [h / g mmol] 
    
    # f_AE stash variable
    rxns4 = [rxns3; "f_T0_stash_rxn"; "f_T_stash_rxn"]
    mets4 = [mets3; "T_stash_met"]
    lb4 = [lb3; f_T0; 0]
    ub4 = [ub3; f_T0; 1000.0]
    c4 = [c3; 0; 0]
    b4 = [b3; 0]
    _col0 = zeros(size(S3, 1))
    S4 = [
        S3          _col0       _col0
        w_T_coes'  1           -1
    ]

    ## - --- -- . .- -.-. -.- .-.-. 
    # MARK: f_Pc
    # TODO/ problem, the model is growth to little
    # See CODE_alter

    # protein allocation constraint
    # f_T + f_UE + f_AE - f_Pc = 0
    
    w_Pc_coes = zeros(size(S4, 2))
    fidx = findfirst(rxns4 .== "f_AE_stash_rxn")
    w_Pc_coes[fidx] = 1                             # [g / gCDW] 
    fidx = findfirst(rxns4 .== "f_UE_stash_rxn")
    w_Pc_coes[fidx] = 1                             # [g / gCDW] 
    fidx = findfirst(rxns4 .== "f_T_stash_rxn")
    w_Pc_coes[fidx] = 1                             # [g / gCDW] 

    # f_AE stash variable
    rxns5 = [rxns4; "f_Pc_stash_rxn"]
    mets5 = [mets4; "Pc_stash_met"]
    lb5 = [lb4; 0]
    ub5 = [ub4; f_Pc]
    c5 = [c4; 0]
    b5 = [b4; 0]
    _col0 = zeros(size(S4, 1))
    S5 = [
        S4          _col0
        w_Pc_coes'  -1
    ]

    # MARK: MetNet
    net0 = MetNet(;S=S5, b=b5, lb=lb5, ub=ub5, c=c5, rxns=rxns5, mets=mets5)

    extras = Dict(
        "f_UE0" => f_UE0,
        "u_ider" => u_ider,
        "u_max" => u_max,
        "f_T0" => f_T0,
        "w_T" => w_T,
        "f_Pc" => f_Pc,
    )

    return (net0, extras)
end

function _PAMModel_alterProteomeRegulationPatterns2021_builder(build_args...)

    # load
    net0, extras0 = __PAMModel_alterProteomeRegulationPatterns2021_builder(build_args...)
    net = _common_format(net0)
    
    extras!(net, "BIOM", "BIOMASS_Ec_iML1515_WT_75p37M_f")
    extras!(net, "EX_GLC", "EX_glc__D_e_b")
    extras!(net, "EX_NH4", "EX_nh4_e_b")
    extras!(net, "EX_GLU", "EX_glu__L_e_b")
    extras!(net, "EX_O2", "EX_o2_e_b")
    extras!(net, "EX_CO2", "EX_co2_e_b")
    extras!(net, "ATPM", "ATPM")
    for (key, val0) in extras0
        extras!(net, key, val0)
    end

    linear_weights!(net, "BIOMASS_Ec_iML1515_WT_75p37M_f", 1.0)

    return net
end

function _register_PAMModel_alterProteomeRegulationPatterns2021()
    register_network!(
        "PAMModel-alterProteomeRegulationPatterns2021", 
        _PAMModel_alterProteomeRegulationPatterns2021_builder;
        use_cache = false,
        source = "https://doi.org/10.1128/msystems.00625-20", 
        desc = "iML1515 model plus proteome allocation contraints"
    )
end