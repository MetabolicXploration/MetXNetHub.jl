function _toy_net_builder()
    net = Dict()
    net[:S] = 
    # rxns: gt    ferm  resp  ldh   lt   biom    atpm  # mets
        [   1.0  -1.0   0.0   0.0   0.0   0.0    0.0;  #  G
            0.0   2.0  18.0   0.0   0.0  -55.0  -5.0;  #  E
            0.0   2.0  -1.0  -1.0   0.0   0.0    0.0;  #  P
            0.0   0.0   0.0   1.0   1.0   0.0    0.0;  #  L
        ]
    
    net[:mets] = ["G", "E", "P", "L"]
    net[:b] =    [0.0, 0.0, 0.0, 0.0] # demand
    
    net[:metNames] = ["Glucose", "Energy", "Intermediate Product" , "Lactate"];
    
    net[:rxns] = ["Ex_glc"  ,"ferm" ,"resp" , "ldh" ,  "Ex_lac" , "biom", "atpm"];
    net[:lb] =   [0.0   , 0.0   , 0.0   ,  0.0  , -100.0,   0.0,     0.5];
    net[:ub] =   [10.0 , 100.0 , 100.0 , 100.0 ,    0.0, 100.0,    100.0];
    net[:c] =    [ 0.0 ,   0.0 ,   0.0 ,   0.0 ,    0.0,   1.0,      0.0];
    net[:rxnNames] = ["Glucose transport", "Fermentation", "Respiration", 
        "Lactate DH", "Lactate transport", "Biomass production rate", "atp demand"];

    net = MetNet(;net...)

    # extras
    set_extra!(net, "BIOM", "biom")
    
    return net
end

function _register_toy_net()
    register_network!("toy_net", _toy_net_builder;
        use_cache = false,
        source = "Home made ;)", 
        desc = "A small (4,7) model with mimicking glycolysis/fermentation + respiration pathways."
    )
end