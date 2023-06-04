function _raw_file(f, fs...)
    rawfile = joinpath(pkgdir(MetXNetHub), "data", "raw", f, fs...) 
    isfile(rawfile) || error("Raw file missing, file: $rawfile")
    return rawfile
end

function _load_raw_model(f, fs...; kwargs...)
    return MetXGEMs.load_net(_raw_file(f, fs...); kwargs...)
end

function _common_format(net0::MetNet; flx_inf = 1000.0)
    # sparsity
    net = MetNet(net0; 
        S = _sparse(net0.S),
        lb = _dense(net0.lb),
        ub = _dense(net0.ub),
        c = _dense(net0.c),
        b = _dense(net0.b),
        rxnGeneMat = _sparse(net0.rxnGeneMat),
    )
    # bounds
    clampbounds!(net, -flx_inf, flx_inf)

    return net

end