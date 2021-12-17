_applyind(a::Number, _::Vector) = a
_applyind(a::Vector, ind::Vector) = a[ind]
_applyind(a::Matrix, ind::Vector) = a[:, ind]

function _applyind!(p::AbstractParticles, ind::Vector)
    for key in keys(p)
        p[key] = p[key][ind]
    end

    return p
end

function _applyind_to_copy(p::AbstractParticles, ind::Vector, affect=())
    # TODO: create Dict only for affected props instead of copying full Dict
    p = copy(p)

    if length(affect) == 0
        affect = keys(p)
    end

    for key in affect
        p[key] = _applyind(p[key], ind)
    end

    # remove keys in new Dict that weren't sorted
    # TODO: won't need this when Dict only created for affected
    for key in setdiff(keys(p), affect)
        delete!(get_props(p), key)
    end

    return p
end

function ind_of_a_in_set(a::AbstractVector, set::AbstractVector)
    if (issorted(a) && issorted(set))
        return ind_of_a_in_set_sorted(a, set)
    else
        return ind_of_a_in_set(a, Set(set))
    end
end

ind_of_a_in_set(a::AbstractVector, set::AbstractSet) = findall(in.(a, (set,)))

function ind_of_a_in_set_sorted(a::AbstractVector, set::AbstractVector)
    if isempty(a) || isempty(set)
        return Int64[]
    end

    na = length(a)
    nset = length(set)

    ind_all = Vector{Int64}(undef, min(na, nset))

    # NOTE: this is an alternative algorithm, which performs slightly worse than the below
    # iind = 0
    # iset = 1
    # for ia in 1:na
        # while a[ia] > set[iset] && iset < nset
            # iset += 1
        # end
        # if a[ia] == set[iset]
            # iind += 1
            # iset += 1
            # ind_all[iind] = ia
        # end
    # end

    iind = 0
    ia = 1
    iset = 1

    while true
        if a[ia] == set[iset]
            iind += 1
            ind_all[iind] = ia
            iset += 1
            ia += 1
            ((ia > na) || (iset > nset)) && break
        else
            if a[ia] < set[iset]
                ia += 1
                ia > na && break
            else
                iset += 1
                iset > nset && break
            end
        end
    end
    
    return resize!(ind_all, iind)
end
