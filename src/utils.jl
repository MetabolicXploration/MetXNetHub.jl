# make stuff dense
_dense_vecs(v) = (v isa SparseVector) ? Vector(v) : v

_dense_vecs(net::MetNet) = MetNet(net; 
    lb = _dense_vecs(net.lb),
    ub = _dense_vecs(net.ub),
    c = _dense_vecs(net.c),
    b = _dense_vecs(net.b),
)

function _load_raw_model(rawfilename)

    rawfile = joinpath(pkgdir(MetXNetHub), "data", "raw", basename(rawfilename)) 
    isfile(rawfile) || error("Raw file missing, file: $rawfile")

    return MetXBase.load_net(rawfile)
end