# Partial derivatives dy/dx_i for y of type V.
# This is a kind of sparse vector, whose length is
# the total number of inputs x_i encountered during
# the course of a computation. Length is stored as a reference;
# its value will be mutated whenever a new input is encountered.
struct Partials{V} <: AbstractVector{V}
    values::Dict{Int,V}
    length::Ref{Int}
end

const NO_LENGTH = Ref(0)

@inline valtype(::Partials{V}) where {V} = V
@inline valtype(::Type{Partials{V}}) where {V} = V

@inline npartials(p::Partials) = p.length[]
@inline npartials(p::Type{Partials{V}}) where {V} = p.length[]

@inline Base.length(p::Partials) = p.length[]
@inline Base.size(p::Partials) = (p.length[],)

@inline Base.@propagate_inbounds Base.getindex(partials::Partials{V}, i::Int) where {V} = get(partials.values, i, zero(V))

Base.iterate(partials::Partials) = length(partials) == 0 ? nothing : (partials[1], 2)
Base.iterate(partials::Partials, i) = length(partials) < i ? nothing : (partials[i], i + 1)

Base.IndexStyle(::Type{<:Partials}) = IndexLinear()


#####################
# Generic Functions #
#####################

# TODO: detect when all values are zero?
@inline Base.iszero(partials::Partials) = isempty(partials.values)

# TODO: how to handle length field here?
# Commenting out for now.
@inline Base.zero(partials::Partials) = zero(typeof(partials))
@inline Base.zero(::Type{Partials{V}}) where {V} = Partials{V}(Dict{Int, V}(), NO_LENGTH)

# Unclear what a "1" should mean, absent a specified length.
# @inline Base.one(partials::Partials) = one(typeof(partials))
# @inline Base.one(::Type{Partials{V}}) where {V} = Partials{V}(one_tuple(NTuple{N,V}))

# @inline Random.rand(partials::Partials) = rand(typeof(partials))
# @inline Random.rand(::Type{Partials{N,V}}) where {N,V} = Partials{N,V}(rand_tuple(NTuple{N,V}))
# @inline Random.rand(rng::AbstractRNG, partials::Partials) = rand(rng, typeof(partials))
# @inline Random.rand(rng::AbstractRNG, ::Type{Partials{N,V}}) where {N,V} = Partials{N,V}(rand_tuple(rng, NTuple{N,V}))

Base.isequal(a::Partials, b::Partials) = isequal(a.values, b.values)
Base.:(==)(a::Partials, b::Partials) = a.values == b.values

const PARTIALS_HASH = hash(Partials)

Base.hash(partials::Partials) = hash(partials.values, PARTIALS_HASH)
Base.hash(partials::Partials, hsh::UInt64) = hash(hash(partials), hsh)

# @inline Base.copy(partials::Partials) = partials
@inline Base.copy(partials::Partials) = Partials(copy(partials.values), partials.length)

# Not sure how to serialize:
# Base.read(io::IO, ::Type{Partials{V}}) where {V} = Partials{V}(ntuple(i->read(io, V), N))

# function Base.write(io::IO, partials::Partials)
#     for p in partials
#         write(io, p)
#     end
# end


########################
# Conversion/Promotion #
########################

Base.promote_rule(::Type{Partials{A}}, ::Type{Partials{B}}) where {A,B} = Partials{promote_type(A, B)}

Base.convert(::Type{Partials{V}}, partials::Partials) where {V} = Partials{V}(partials.values, partials.length)
Base.convert(::Type{Partials{V}}, partials::Partials{V}) where {V} = partials

########################
# Arithmetic Functions #
########################

function add_partial_dicts(a::Dict{Int, V}, b::Dict{Int, V}) where {V}
    Dict{Int, V}(k => get(a,k,zero(V)) + get(b,k,zero(V)) for k in union(keys(a), keys(b)))
end
function sub_partial_dicts(a::Dict{Int, V}, b::Dict{Int, V}) where {V}
    Dict{Int, V}(k => get(a,k,zero(V)) - get(b,k,zero(V)) for k in union(keys(a), keys(b)))
end
function mul_partial_dict_by_scalar(a::Real, b::Dict{Int, V}) where {V}
    Dict{Int, V}(k => v * a for (k, v) in b)
end
function unary_minus_partial_dict(a::Dict{Int, V}) where {V}
    Dict{Int, V}(k => -v for (k, v) in a)
end
function div_partial_dict_by_scalar(a::Dict{Int, V}, b::Real) where {V}
    Dict{Int, V}(k => v / b for (k, v) in a)
end
# Returns x_a * a + x_b * b.
function mul_dicts(a::Dict{Int, V}, b::Dict{Int, V}, x_a, x_b) where {V}
    return Dict{Int, V}(k => x_a * get(a, k, zero(V)) + x_b * get(b, k, zero(V)) for k in union(keys(a), keys(b)))
end

@inline _resolve_lengths(a::Ref{Int}, b::Ref{Int}) = a === NO_LENGTH ? b : a
@inline Base.:+(a::Partials, b::Partials) = Partials(add_partial_dicts(a.values, b.values), _resolve_lengths(a.length, b.length))
@inline Base.:-(a::Partials, b::Partials) = Partials(sub_partial_dicts(a.values, b.values), _resolve_lengths(a.length, b.length))
@inline Base.:-(partials::Partials) = Partials(unary_minus_partial_dict(partials.values), partials.length)
@inline Base.:*(x::Real, partials::Partials) = Partials(mul_partial_dict_by_scalar(x, partials.values), partials.length)
@inline Base.:*(partials::Partials, x::Real) = Partials(mul_partial_dict_by_scalar(x, partials.values), partials.length)
@inline Base.:/(partials::Partials, x::Real) = Partials(div_partial_dict_by_scalar(partials.values, x), partials.length)
@inline _mul_partials(a::Partials, b::Partials, x_a, x_b) = Partials(mul_dicts(a.values, b.values, x_a, x_b), _resolve_lengths(a.length, b.length))
@inline _div_partials(a::Partials, b::Partials, aval, bval) = _mul_partials(a, b, inv(bval), -(aval / (bval*bval)))


###################
# Pretty Printing #
###################

Base.show(io::IO, p::Partials) = print(io, "Partials: ", p.values)