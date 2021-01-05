#######
# Tag #
#######
struct DiffConfig{T}
    n_inputs :: Ref{Int}
    tag      :: T

    # Currently does not use tag functionality.
    function DiffConfig()::DiffConfig{Nothing}
        new{Nothing}(Ref(0), nothing)
    end
end

function new_dual(c::DiffConfig{T}, v::V) where {T, V}
    next_id = c.n_inputs[] + 1
    c.n_inputs[] = next_id
    Dual{T,V}(v, Partials{V}(Dict{Int, V}(next_id => one(V)), c.n_inputs))
end

function new_dual(c::DiffConfig{T}, v::AbstractArray{V}) where {T, V}
    duals = similar(v, Dual{T,V})
    next_id = c.n_inputs[]
    for i in eachindex(v)
        next_id += 1
        duals[i] = Dual{T,V}(v[i], Partials{V}(Dict{Int, V}(next_id => one(V)), c.n_inputs))
    end
    c.n_inputs[] = next_id
    return duals
end