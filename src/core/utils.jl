function _raw_file(rawfilename)
    rawfile = joinpath(pkgdir(MetXNetHub), "data", "raw", basename(rawfilename)) 
    isfile(rawfile) || error("Raw file missing, file: $rawfile")
    return rawfile
end

function _load_raw_model(rawfilename; kwargs...)
    return MetXGEMs.load_net(_raw_file(rawfilename); kwargs...)
end