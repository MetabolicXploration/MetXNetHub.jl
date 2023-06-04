## ------------------------------------------------------------------
# CACHE
function _cache_file(id, build_args...)
    _hash = hash(hash.([id, build_args...]))
    return joinpath(NETS_DIR, string(basename(id), "_", _hash, ".jls"))
end

## ------------------------------------------------------------------
# PULL NET
function pull_net(id, build_args...; 
        clear_cache = false, 
    )::MetNet

    reg = get_reg(id)

    # cache
    # TODO: No use case yet for caching. Loading is 'fast'
    # cfile = _cache_file(id, build_args...)
    # clear_cache && rm(cfile; force = true)
    # use_cache = get(reg, "use_cache", false)
    
    # cache load
    # TODO: No use case yet for caching. Loading is 'fast'
    # if use_cache && isfile(cfile) 
    #     try; return deserialize(cfile)
    #         catch; rm(cfile; force = true)
    #     end
    # end

    # get net
    net = reg["builder"](build_args...)

    # cache write
    # TODO: No use case yet for caching. Loading is 'fast'
    # if use_cache
    #     mkpath(dirname(cfile))
    #     serialize(cfile, net)
    # end
    
    return net
end

## ------------------------------------------------------------------
function register_network!(id::String, builder::Function; 
        use_cache = true,
        desc...
    )
    haskey(NETS_REG, id) && error("Network '$id' already registered")

    NETS_REG[id] = Dict{String, Any}()

    # desc
    for (k, v) in desc
        NETS_REG[id][string(k)] = v
    end
    
    # use cache
    NETS_REG[id]["use_cache"] = use_cache

    # setup
    NETS_REG[id]["builder"] = builder

    return NETS_REG[id]
end

## ------------------------------------------------------------------
export nethub_status
function nethub_status(id)
    meta = get_reg(id)
    println("-"^40)
    println("id: ", id)
    println("source: ", get(meta, "source", ""))
    println("desc: ", get(meta, "desc", ""))
end

function nethub_status()
    for id in keys(NETS_REG)
        nethub_status(id)
    end
end

## ------------------------------------------------------------------
function clear_cache!(id, build_args...)
    cfile = _cache_file(id, build_args...)
    rm(cfile; force = true)
    return nothing
end

function clear_cache!()
    rm(NETS_DIR; recursive = true, force = true)
    mkpath(NETS_DIR)
    return nothing
end

function get_reg(id) 
    haskey(NETS_REG, id) || error("net '$id' not registered. See 'nethub_status()'")
    return NETS_REG[id]
end

nethubids() = collect(keys(NETS_REG))