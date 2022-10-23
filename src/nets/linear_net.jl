# A simple model
function _linear_S(M)
    N = M + 1
    # S
    S = spzeros(M,N)
    for i in 1:M
        S[i,i] = 1.0
        S[i,i+1] = -1.0
    end
    return S
end

function _linear_net_builder(M::Int = 5)

    net = Dict()
    net[:S] = _linear_S(M)
    M,N = size(net[:S])
    net[:c] = zeros(N)
    net[:b] = zeros(M)
    net[:lb] = zeros(N)
    net[:ub] = ones(N)
    net[:rxns] = [string("Rxn", i) for i in 1:N]
    net[:mets] = [string("Met", i) for i in 1:M]

    
    return MetNet(;net...)
end


function _register_linear_model()
    register_network!("linear_net", _linear_net_builder;
        use_cache = false,
        source = "Home made ;)", 
        desc = "A model with M metabolites in a linear pathway:  -> Met_1 -> ... -> Met_M -> "
    )
end