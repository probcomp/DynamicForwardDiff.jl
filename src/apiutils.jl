#######
# Tag #
#######

const TAGCOUNT = Threads.Atomic{UInt}(0)

# each tag is assigned a unique number
# tags which depend on other tags will be larger
@generated function tagcount(::Type{Val{T}}) where {T}
    :($(Threads.atomic_add!(TAGCOUNT, UInt(1))))
end

function Tag()
    t = gensym()
    tagcount(Val{t}) # trigger generated function
    Val{t}()
end

@inline function ≺(::Type{Val{T1}}, ::Type{Val{T2}}) where {T1,T2}
    tagcount(Val{T1}) < tagcount(Val{T2})
end

struct DiffConfig{T}
    n_inputs :: Ref{Int}
    tag      :: T

    function DiffConfig()
        g = gensym()
        new{Val{g}}(Ref(0), Val{g}())
    end
end

function new_dual(c::DiffConfig{T}, v::V) where {T, V}
    next_id = c.n_inputs[] + 1
    c.n_inputs[] = next_id
    Dual{T,V}(v, Partials{V}(Dict{Int, V}(next_id => one(V)), c.n_inputs))
end