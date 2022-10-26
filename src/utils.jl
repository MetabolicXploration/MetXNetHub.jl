function _load_raw_model(rawfilename)

    rawfile = joinpath(pkgdir(MetXNetHub), "data", "raw", basename(rawfilename)) 
    isfile(rawfile) || error("Raw file missing, file: $rawfile")

    return MetXBase.load_net(rawfile)
end