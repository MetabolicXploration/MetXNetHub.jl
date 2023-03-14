function _load_raw_model(rawfilename; kwargs...)

    rawfile = joinpath(pkgdir(MetXNetHub), "data", "raw", basename(rawfilename)) 
    isfile(rawfile) || error("Raw file missing, file: $rawfile")

    return MetXGEMs.load_net(rawfile; kwargs...)
end