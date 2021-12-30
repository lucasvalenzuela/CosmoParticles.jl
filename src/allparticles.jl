struct AllParticles{APC} <: AbstractParticles where {APC<:AbstractParticleCollection}
    particle_collection::APC
    props::Dict{Symbol,Any}

    AllParticles(pc::APC) where {APC} = new{APC}(pc, Dict{Symbol,Any}())
end


particle_collection(p::AllParticles) = p.particle_collection


function Base.getindex(ap::AllParticles, sym::Symbol)
    # get property index if property already saved or if key does not exist at all to raise Dict error
    if haskey(ap.props, sym) || !haskey(ap, sym)
        return ap.props[sym]
    end

    pc = ap.particle_collection

    # get array dimensions
    dims = 0 # 0 for vectors; for arrays: size(arr, 1)
    for key in keys(pc)
        p = pc[key]
        if haskey(p, sym)
            p[sym] isa AbstractVector && break
            if p[sym] isa AbstractMatrix
                dims = size(p[sym], 1)
                break
            end
        end
    end

    # create array of particle property arrays
    arrs = []

    # iterate through particles
    for p in values(pc)
        # fill array with missing if property does not exist
        if haskey(p, sym)
            # fill array with scalar value if property is not an array
            if p[sym] isa AbstractArray
                arr = p[sym]
            else
                if dims == 0
                    arr = Fill(p[sym], particle_number(p))
                else
                    arr = Fill(p[sym], dims, particle_number(p))
                end
            end
            push!(arrs, arr)
        else
            if dims == 0
                arr = Fill(missing, particle_number(p))
            else
                arr = Fill(missing, dims, particle_number(p))
            end
            push!(arrs, arr)
        end
    end

    cat_func = dims == 0 ? vcat : hcat
    lazy = ApplyArray(cat_func, arrs...)

    ap.props[sym] = lazy
    return lazy
end

function Base.setindex!(::AllParticles, _, ::Symbol)
    error(
        "Setting a property of an AllParticles object is not allowed. Try converting it via Particles(allparticles) first.",
    )
end

function Base.keys(ap::AllParticles)
    keys_arr = [keys(p) for p in values(ap.particle_collection)]
    length(keys_arr) == 1 && return keys_arr[1]
    return union(keys_arr...) |> collect
end

Base.haskey(p::AllParticles, key) = haskey(p.props, key) || key in keys(p)

Base.values(p::AllParticles) = Iterators.map(key -> getindex(p, key), keys(p))

function Base.copy(::AllParticles)
    error(
        "Copying an AllParticles object is not allowed. Try converting it via Particles(allparticles) instead.",
    )
    ## NOTE: potential way of making copy work
    # pnew = Particles(p.particle_collection)
    # copy!(pnew.props, p.props)
    # return pnew
end

function Base.copy!(::AllParticles, ::AllParticles, props=())
    error(
        "Copying particles to an AllParticles object is not allowed. Try converting it via Particles(allparticles) instead.",
    )
end

function Base.copy!(::AllParticles, ::AbstractParticles, props=())
    error(
        "Copying particles to an AllParticles object is not allowed. Try converting it via Particles(allparticles) instead.",
    )
end

function Base.copy!(dst::AbstractParticles, src::AllParticles, props=())
    empty!(dst)
    for key in keys(src)
        if isempty(props) || key in props
            val = src[key]
            dst[key] = Array(val)
        end
    end
    return dst
end

function Base.empty(::AllParticles)
    error(
        "Creating an empty AllParticles object is not supported. Try emptying the underlying particle collection instead.",
    )
end

function Base.empty!(::AllParticles)
    error(
        "Emptying an AllParticles object is not allowed. Try emptying the underlying particle collection instead.",
    )
end

function Base.isempty(p::AllParticles)
    for p in values(p.particle_collection)
        isempty(p) || return false
    end
    return true
end

Base.:(==)(p1::AllParticles, p2::AllParticles) = p1.particle_collection == p2.particle_collection

particle_name(::AllParticles) = "Particles"

function particle_number(ap::AllParticles)
    n = 0
    for p in values(ap.particle_collection)
        n += particle_number(p)
    end
    return n
end

function Base.propertynames(p::AllParticles; private=false)
    if private
        return [fieldnames(AllParticles) |> collect; keys(p) |> collect]
    else
        return keys(p) |> collect
    end
end

function Base.show(io::IO, mime::MIME"text/plain", p::AllParticles)
    printstyled(io, "all"; bold=true)
    print(io, ": ")
    show_properties(io, mime, p)
end


Particles(p::AllParticles, props=()) = copy!(Particles(:all), p, props)

function applyind!(::AllParticles, ::AbstractVector)
    error(
        "Indexing into an AllParticles object is not allowed. Try converting it via Particles(allparticles) first.",
    )
end

applyind(p::AllParticles, ind::AbstractVector) = applyind!(p, ind)
