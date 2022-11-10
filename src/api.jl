function pull_net(id, build_args...; 
        clear_cache = false, 
    )

    reg = get_reg(id)

    # cache
    cachefile = reg["cachefile"]
    clear_cache && clear_cache!(id)
    if isfile(cachefile) 
        try; return deserialize(cachefile)
            catch; clear_cache!(id)
        end
    end

    # get net
    net = reg["builder"](build_args...)

    # cache
    isempty(cachefile) || serialize(cachefile, net)
    
    return net
end

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
    
    # paths
    NETS_REG[id]["cachefile"] = use_cache ? joinpath(NETS_DIR, string(id, ".jls")) : ""
    
    # setup
    NETS_REG[id]["builder"] = builder

    return NETS_REG[id]
end

function hub_status()
    for (id, meta) in NETS_REG
        println("-"^40)
        println("id: ", id)
        println("source: ", get(meta, "source", ""))
        println("desc: ", get(meta, "desc", ""))
    end
end

function clear_cache!(id = nothing)

    if isnothing(id)
        rm(NETS_DIR; recursive = true, force = true)
        mkpath(NETS_DIR)
    else
        reg = get_reg(id)
        cachefile = get(reg, "cachefile", "")
        rm(cachefile; force = true)
    end
    return nothing
end

function get_reg(id) 
    haskey(NETS_REG, id) || error("net '$id' not registered. See 'hub_status()'")
    return NETS_REG[id]
end

nethubids() = collect(keys(NETS_REG))