## ------------------------------------------------------------------
# 4D model
function _toy_net4D_builder()
    net = Dict()
    # TODO: make realistic S (make it growth like ecoli_core)
    net[:S] = 
    # rxns: gt   at  anap   ferm  resp  ldh   lt   biom    atpm  # mets
        [   1.0  0.0  0.0  -1.0   0.0   0.0   0.0  -16.2    0.0;  #  G
            0.0  1.0 -1.0   0.0   0.0   0.0   0.0  -9.53    0.0;  #  A
            0.0  0.0  0.0   2.0  18.0   0.0   0.0  -55.0   -5.0;  #  E
            0.0  0.0  1.5   2.0  -1.0  -1.0   0.0   0.0     0.0;  #  P
            0.0  0.0  0.0   0.0   0.0   1.0   1.0   0.0     0.0;  #  L
        ]
    
    net[:mets] = ["G", "A", "E", "P", "L"]
    net[:b] =    [0.0, 0.0, 0.0, 0.0, 0.0] # demand
    
    net[:metNames] = ["Glucose", "Aminoacids", "Energy", "Intermediate Product" , "Lactate"];
    
    net[:rxns] = [ "Ex_glc", "Ex_aa", "anap", "ferm", "resp", "ldh" , "Ex_lac" , "biom" , "atpm" ];
    net[:lb]   = [   0.0   ,   0.0  ,   0.0 ,   0.0 ,   0.0 ,   0.0 ,  -100.0  ,   0.0  ,   0.5  ];
    net[:ub]   = [  10.0   , 100.0  , 100.0 , 100.0 , 100.0 , 100.0 ,     0.0  , 100.0  ,  100.0 ];
    net[:c]    = [   0.0   ,   0.0  ,   0.0 ,   0.0 ,   0.0 ,   0.0 ,     0.0  ,   1.0  ,    0.0 ];
    net[:rxnNames] = ["Glucose transport", "Aminoacids transport", "Anaplerotic reactions", "Fermentation", "Respiration", 
        "Lactate DH", "Lactate transport", "Biomass production rate", "atp demand"];

    net = MetNet(;net...)

    # extras
    extras!(net, "BIOM", "biom")
    extras!(net, "EX_GLC", "Ex_glc")
    extras!(net, "EX_AA", "Ex_aa")
    
    linear_coefficients!(net, "biom", 1.0)
    
    return net
end

function _register_toy_net4D()
    register_network!("toy_net4D", _toy_net4D_builder;
        use_cache = false,
        source = "Home made ;)", 
        desc = "A small (5, 9) model with similar to `toy_net` but with 4 free variables."
    )
end