"""
DynamicForwardDiff.can_dual(V::Type)
Determines whether the type V is allowed as the scalar type in a
Dual. By default, only `<:Real` types are allowed.
"""
can_dual(::Type{<:Real}) = true
can_dual(::Type) = false

# T is a "tag."
struct Dual{T,V} <: Real
    value::V
    partials::Partials{V}
    function Dual{T,V}(value::V, partials::Partials{V}) where {T,V}
        can_dual(V) || throw_cannot_dual(V)
        new{T,V}(value, partials)
    end
end


##############
# Exceptions #
##############

struct DualMismatchError{A,B} <: Exception
    a::A
    b::B
end

Base.showerror(io::IO, e::DualMismatchError{A,B}) where {A,B} =
print(io, "Cannot determine ordering of Dual tags $(e.a) and $(e.b)")

@noinline function throw_cannot_dual(V::Type)
    throw(ArgumentError("Cannot create a dual over scalar type $V." *
    " If the type behaves as a scalar, define FowardDiff.can_dual."))
end

"""
DynamicForwardDiff.≺(a, b)::Bool
Determines the order in which tagged `Dual` objects are composed. If true, then `Dual{b}`
objects will appear outside `Dual{a}` objects.

This is important when working with nested differentiation: currently, only the outermost
tag can be extracted, so it should be used in the _innermost_ function.
"""
≺(a,b) = throw(DualMismatchError(a, b))


################
# Constructors #
################

@inline Dual{T}(value::V, partials::Partials{V}) where {T,V} = Dual{T,V}(value, partials)

@inline function Dual{T}(value::A, partials::Partials{B}) where {T,A,B}
    C = promote_type(A, B)
    return Dual{T}(convert(C, value), convert(Partials{C}, partials))
end
@inline Dual(args...) = Dual{Nothing}(args...)

# we define these special cases so that the "constructor <--> convert" pun holds for `Dual`
@inline Dual{T,V}(x::Dual{T,V}) where {T,V} = x
@inline Dual{T,V}(x) where {T,V} = convert(Dual{T,V}, x)


##############################
# Utility/Accessor Functions #
##############################

@inline value(x) = x
@inline value(d::Dual) = d.value
@inline value(d::Array{<:Dual}) = map(value, d)

@inline value(::Type{T}, x) where T = x
@inline value(::Type{T}, d::Dual{T}) where T = value(d)
@inline function value(::Type{T}, d::Dual{S}) where {T,S}
    if S ≺ T
        d
    else
        throw(DualMismatchError(T, S))
    end
end

# Question: where to get the `length` tag?
@inline partials(x) = Partials{typeof(x)}(Dict{Int,typeof(x)}(), NO_LENGTH)
@inline partials(d::Dual) = d.partials
@inline partials(x, i...) = zero(x)
@inline Base.@propagate_inbounds partials(d::Dual, i) = d.partials[i]
@inline Base.@propagate_inbounds partials(d::Dual, i, j) = partials(d, i).partials[j]
@inline Base.@propagate_inbounds partials(d::Dual, i, j, k...) = partials(partials(d, i, j), k...)

@inline Base.@propagate_inbounds partials(::Type{T}, x, i...) where T = partials(x, i...)
@inline Base.@propagate_inbounds partials(::Type{T}, d::Dual{T}, i...) where T = partials(d, i...)
@inline function partials(::Type{T}, d::Dual{S}, i...) where {T,S}
    if S ≺ T
        zero(d)
    else
        throw(DualMismatchError(T, S))
    end
end

@inline npartials(d::Dual{T,V}) where {T,V} = npartials(d.partials)

@inline order(::Type{V}) where {V} = 0
@inline order(::Type{Dual{T,V}}) where {T,V} = 1 + order(V)

@inline valtype(::V) where {V} = V
@inline valtype(::Type{V}) where {V} = V
@inline valtype(::Dual{T,V}) where {T,V} = V
@inline valtype(::Type{Dual{T,V}}) where {T,V} = V

@inline tagtype(::V) where {V} = Nothing
@inline tagtype(::Type{V}) where {V} = Nothing
@inline tagtype(::Dual{T,V}) where {T,V} = T
@inline tagtype(::Type{Dual{T,V}}) where {T,V} = T


####################################
# N-ary Operation Definition Tools #
####################################


replace_match!(f, ismatch, x) = x

function replace_match!(f, ismatch, lines::AbstractArray)
    for i in eachindex(lines)
        line = lines[i]
        if ismatch(line)
            lines[i] = f(line)
        elseif isa(line, Expr)
            replace_match!(f, ismatch, line.args)
        end
    end
    return lines
end

# This is basically a workaround that allows one to use CommonSubexpressions.cse, but with
# qualified bindings. This is not guaranteed to be correct if the input expression contains
# field accesses.
function qualified_cse!(expr)
    placeholders = Dict{Symbol,Expr}()
    replace_match!(x -> isa(x, Expr) && x.head == :(.), expr.args) do x
        placeholder = Symbol("#$(hash(x))")
        placeholders[placeholder] = x
        placeholder
    end
    cse_expr = CommonSubexpressions.cse(expr, warn=false)
    replace_match!(x -> haskey(placeholders, x), cse_expr.args) do x
        placeholders[x]
    end
    return cse_expr
end

const AMBIGUOUS_TYPES = (AbstractFloat, Irrational, Integer, Rational, Real, RoundingMode)

macro define_binary_dual_op(f, xy_body, x_body, y_body)
    DFD = DynamicForwardDiff
    defs = quote
        @inline $(f)(x::$DFD.Dual{Txy}, y::$DFD.Dual{Txy}) where {Txy} = $xy_body
        @inline $(f)(x::$DFD.Dual{Tx}, y::$DFD.Dual{Ty}) where {Tx,Ty} = Ty ≺ Tx ? $x_body : $y_body
    end
    for R in AMBIGUOUS_TYPES
        expr = quote
            @inline $(f)(x::$DFD.Dual{Tx}, y::$R) where {Tx} = $x_body
            @inline $(f)(x::$R, y::$DFD.Dual{Ty}) where {Ty} = $y_body
        end
        append!(defs.args, expr.args)
    end
    return esc(defs)
end

macro define_ternary_dual_op(f, xyz_body, xy_body, xz_body, yz_body, x_body, y_body, z_body)
    DFD = DynamicForwardDiff
    defs = quote
        @inline $(f)(x::$DFD.Dual{Txyz}, y::$DFD.Dual{Txyz}, z::$DFD.Dual{Txyz}) where {Txyz} = $xyz_body
        @inline $(f)(x::$DFD.Dual{Txy}, y::$DFD.Dual{Txy}, z::$DFD.Dual{Tz}) where {Txy,Tz} = Tz ≺ Txy ? $xy_body : $z_body
        @inline $(f)(x::$DFD.Dual{Txz}, y::$DFD.Dual{Ty}, z::$DFD.Dual{Txz}) where {Txz,Ty} = Ty ≺ Txz ? $xz_body : $y_body
        @inline $(f)(x::$DFD.Dual{Tx}, y::$DFD.Dual{Tyz}, z::$DFD.Dual{Tyz}) where {Tx,Tyz} = Tyz ≺ Tx ? $x_body  : $yz_body
        @inline function $(f)(x::$DFD.Dual{Tx}, y::$DFD.Dual{Ty}, z::$DFD.Dual{Tz}) where {Tx,Ty,Tz}
            if Tz ≺ Tx && Ty ≺ Tx
                $x_body
            elseif Tz ≺ Ty
                $y_body
            else
                $z_body
            end
        end
    end
    for R in AMBIGUOUS_TYPES
        expr = quote
            @inline $(f)(x::$DFD.Dual{Txy}, y::$DFD.Dual{Txy}, z::$R) where {Txy} = $xy_body
            @inline $(f)(x::$DFD.Dual{Tx}, y::$DFD.Dual{Ty}, z::$R)  where {Tx,Ty} = Ty ≺ Tx ? $x_body : $y_body
            @inline $(f)(x::$DFD.Dual{Txz}, y::$R, z::$DFD.Dual{Txz}) where {Txz} = $xz_body
            @inline $(f)(x::$DFD.Dual{Tx}, y::$R, z::$DFD.Dual{Tz}) where {Tx,Tz} = Tz ≺ Tx ? $x_body : $z_body
            @inline $(f)(x::$R, y::$DFD.Dual{Tyz}, z::$DFD.Dual{Tyz}) where {Tyz} = $yz_body
            @inline $(f)(x::$R, y::$DFD.Dual{Ty}, z::$DFD.Dual{Tz}) where {Ty,Tz} = Tz ≺ Ty ? $y_body : $z_body
        end
        append!(defs.args, expr.args)
        for Q in AMBIGUOUS_TYPES
            Q === R && continue
            expr = quote
                @inline $(f)(x::$DFD.Dual{Tx}, y::$R, z::$Q) where {Tx} = $x_body
                @inline $(f)(x::$R, y::$DFD.Dual{Ty}, z::$Q) where {Ty} = $y_body
                @inline $(f)(x::$R, y::$Q, z::$DFD.Dual{Tz}) where {Tz} = $z_body
            end
            append!(defs.args, expr.args)
        end
        expr = quote
            @inline $(f)(x::$DFD.Dual{Tx}, y::$R, z::$R) where {Tx} = $x_body
            @inline $(f)(x::$R, y::$DFD.Dual{Ty}, z::$R) where {Ty} = $y_body
            @inline $(f)(x::$R, y::$R, z::$DFD.Dual{Tz}) where {Tz} = $z_body
        end
        append!(defs.args, expr.args)
    end
    return esc(defs)
end

function unary_dual_definition(M, f)
    DFD = DynamicForwardDiff
    Mf = M == :Base ? f : :($M.$f)
    work = qualified_cse!(quote
        val = $Mf(x)
        deriv = $(DiffRules.diffrule(M, f, :x))
    end)
    return quote
        @inline function $M.$f(d::$DFD.Dual{T}) where T
            x = $DFD.value(d)
            $work
            return $DFD.Dual{T}(val, deriv * $DFD.partials(d))
        end
    end
end

function binary_dual_definition(M, f)
    FD = DynamicForwardDiff
    dvx, dvy = DiffRules.diffrule(M, f, :vx, :vy)
    Mf = M == :Base ? f : :($M.$f)
    xy_work = qualified_cse!(quote
        val = $Mf(vx, vy)
        dvx = $dvx
        dvy = $dvy
    end)
    dvx, _ = DiffRules.diffrule(M, f, :vx, :y)
    x_work = qualified_cse!(quote
        val = $Mf(vx, y)
        dvx = $dvx
    end)
    _, dvy = DiffRules.diffrule(M, f, :x, :vy)
    y_work = qualified_cse!(quote
        val = $Mf(x, vy)
        dvy = $dvy
    end)
    expr = quote
        @define_binary_dual_op(
        $M.$f,
        begin
            vx, vy = $FD.value(x), $FD.value(y)
            $xy_work
            return $FD.Dual{Txy}(val, $FD._mul_partials($FD.partials(x), $FD.partials(y), dvx, dvy))
        end,
        begin
            vx = $FD.value(x)
            $x_work
            return $FD.Dual{Tx}(val, dvx * $FD.partials(x))
        end,
        begin
            vy = $FD.value(y)
            $y_work
            return $FD.Dual{Ty}(val, dvy * $FD.partials(y))
        end
        )
    end
    return expr
end


#####################
# Generic Functions #
#####################

Base.copy(d::Dual{T,V}) where {T, V} = Dual{T,V}(d.value, copy(d.partials))

Base.eps(d::Dual) = eps(value(d))
Base.eps(::Type{D}) where {D <: Dual} = eps(valtype(D))

function Base.nextfloat(d::DynamicForwardDiff.Dual{T,V}) where {T,V}
    DynamicForwardDiff.Dual{T}(nextfloat(d.value), d.partials)
end

function Base.prevfloat(d::DynamicForwardDiff.Dual{T,V}) where {T,V}
    DynamicForwardDiff.Dual{T}(prevfloat(d.value), d.partials)
end

Base.rtoldefault(::Type{D}) where {D <: Dual} = Base.rtoldefault(valtype(D))

Base.floor(::Type{R}, d::Dual) where {R <: Real} = floor(R, value(d))
Base.floor(d::Dual) = floor(value(d))

Base.ceil(::Type{R}, d::Dual) where {R <: Real} = ceil(R, value(d))
Base.ceil(d::Dual) = ceil(value(d))

Base.trunc(::Type{R}, d::Dual) where {R <: Real} = trunc(R, value(d))
Base.trunc(d::Dual) = trunc(value(d))

Base.round(::Type{R}, d::Dual) where {R <: Real} = round(R, value(d))
Base.round(d::Dual) = round(value(d))

Base.hash(d::Dual) = hash(value(d))
Base.hash(d::Dual, hsh::UInt) = hash(value(d), hsh)

# function Base.read(io::IO, ::Type{Dual{T,V,N}}) where {T,V,N}
#     value = read(io, V)
#     partials = read(io, Partials{N,V})
#     return Dual{T,V,N}(value, partials)
# end

# function Base.write(io::IO, d::Dual)
#     write(io, value(d))
#     write(io, partials(d))
# end

@inline Base.zero(d::Dual) = zero(typeof(d))
@inline Base.zero(::Type{Dual{T,V}}) where {T,V} = Dual{T}(zero(V), zero(Partials{V}))

@inline Base.one(d::Dual) = one(typeof(d))
@inline Base.one(::Type{Dual{T,V}}) where {T,V} = Dual{T}(one(V), zero(Partials{V}))

@inline Random.rand(rng::AbstractRNG, d::Dual) = rand(rng, value(d))
@inline Random.rand(::Type{Dual{T,V}}) where {T,V} = Dual{T}(rand(V), zero(Partials{V}))
@inline Random.rand(rng::AbstractRNG, ::Type{Dual{T,V}}) where {T,V} = Dual{T}(rand(rng, V), zero(Partials{V}))
@inline Random.randn(::Type{Dual{T,V}}) where {T,V} = Dual{T}(randn(V), zero(Partials{V}))
@inline Random.randn(rng::AbstractRNG, ::Type{Dual{T,V}}) where {T,V} = Dual{T}(randn(rng, V), zero(Partials{V}))
@inline Random.randexp(::Type{Dual{T,V}}) where {T,V} = Dual{T}(randexp(V), zero(Partials{V}))
@inline Random.randexp(rng::AbstractRNG, ::Type{Dual{T,V}}) where {T,V} = Dual{T}(randexp(rng, V), zero(Partials{V}))

# ------------#
# Predicates #
# ------------#

isconstant(d::Dual) = iszero(partials(d))

const UNARY_PREDICATES = Symbol[:isinf, :isnan, :isfinite, :iseven, :isodd, :isreal, :isinteger]

const BINARY_PREDICATES = Symbol[:isequal, :isless, :<, :>, :(==), :(!=), :(<=), :(>=)]

for pred in UNARY_PREDICATES
    @eval Base.$(pred)(d::Dual) = $(pred)(value(d))
end

for pred in BINARY_PREDICATES
    @eval begin
        @define_binary_dual_op(
        Base.$(pred),
        $(pred)(value(x), value(y)),
        $(pred)(value(x), y),
        $(pred)(x, value(y))
        )
    end
end


########################
# Promotion/Conversion #
########################

Base.@pure function Base.promote_rule(::Type{Dual{T1,V1}},
    ::Type{Dual{T2,V2}}) where {T1,V1,T2,V2}
    # V1 and V2 might themselves be Dual types
    if T2 ≺ T1
        Dual{T1,promote_type(V1, Dual{T2,V2})}
    else
        Dual{T2,promote_type(V2, Dual{T1,V1})}
    end
end

function Base.promote_rule(::Type{Dual{T,A}},
    ::Type{Dual{T,B}}) where {T,A,B}
    return Dual{T,promote_type(A, B)}
end

for R in (Irrational, Real, BigFloat, Bool)
    if isconcretetype(R) # issue #322 in ForwardDiff.jl repo
        @eval begin
            Base.promote_rule(::Type{$R}, ::Type{Dual{T,V}}) where {T,V} = Dual{T,promote_type($R, V)}
            Base.promote_rule(::Type{Dual{T,V}}, ::Type{$R}) where {T,V} = Dual{T,promote_type(V, $R)}
        end
    else
        @eval begin
            Base.promote_rule(::Type{R}, ::Type{Dual{T,V}}) where {R <: $R,T,V} = Dual{T,promote_type(R, V)}
            Base.promote_rule(::Type{Dual{T,V}}, ::Type{R}) where {T,V,R <: $R} = Dual{T,promote_type(V, R)}
        end
    end
end

Base.convert(::Type{Dual{T,V}}, d::Dual{T}) where {T,V} = Dual{T}(convert(V, value(d)), convert(Partials{V}, partials(d)))
Base.convert(::Type{Dual{T,V}}, x) where {T,V} = Dual{T}(convert(V, x), zero(Partials{V}))
Base.convert(::Type{Dual{T,V}}, x::Number) where {T,V} = Dual{T}(convert(V, x), zero(Partials{V}))
Base.convert(::Type{D}, d::D) where {D <: Dual} = d

Base.float(d::Dual{T,V}) where {T,V} = convert(Dual{T,promote_type(V, Float16)}, d)
Base.AbstractFloat(d::Dual{T,V}) where {T,V} = convert(Dual{T,promote_type(V, Float16)}, d)

###################################
# General Mathematical Operations #
###################################

@inline Base.conj(d::Dual) = d

for (M, f, arity) in DiffRules.diffrules()
    in((M, f), ((:Base, :^), (:Base, :/), (:Base, :+), (:Base, :-))) && continue
    M == :NaNMath && continue
    if arity == 1
        eval(unary_dual_definition(M, f))
    elseif arity == 2
        eval(binary_dual_definition(M, f))
    else
        error("DynamicForwardDiff currently only knows how to autogenerate Dual definitions for unary and binary primitives.")
    end
end

#################
# Special Cases #
#################

# +/- #
# -----#

@define_binary_dual_op(
Base.:+,
begin
    vx, vy = value(x), value(y)
    Dual{Txy}(vx + vy, partials(x) + partials(y))
end,
Dual{Tx}(value(x) + y, partials(x)),
Dual{Ty}(x + value(y), partials(y))
)

@define_binary_dual_op(
Base.:-,
begin
    vx, vy = value(x), value(y)
    Dual{Txy}(vx - vy, partials(x) - partials(y))
end,
Dual{Tx}(value(x) - y, partials(x)),
Dual{Ty}(x - value(y), -partials(y))
)

@inline Base.:-(d::Dual{T}) where {T} = Dual{T}(-value(d), -partials(d))

# * #
# ---#

@inline Base.:*(d::Dual, x::Bool) = x ? d : (signbit(value(d)) == 0 ? zero(d) : -zero(d))
@inline Base.:*(x::Bool, d::Dual) = d * x

# / #
# ---#

# We can't use the normal diffrule autogeneration for this because (x/y) === (x * (1/y))
# doesn't generally hold true for floating point; see issue #264
@define_binary_dual_op(
Base.:/,
begin
    vx, vy = value(x), value(y)
    Dual{Txy}(vx / vy, _div_partials(partials(x), partials(y), vx, vy))
end,
Dual{Tx}(value(x) / y, partials(x) / y),
begin
    v = value(y)
    divv = x / v
    Dual{Ty}(divv, -(divv / v) * partials(y))
end
)


# exponentiation #
# ----------------#

for f in (:(Base.:^),)
    @eval begin
        @define_binary_dual_op(
        $f,
        begin
            vx, vy = value(x), value(y)
            expv = ($f)(vx, vy)
            powval = vy * ($f)(vx, vy - 1)
            if isconstant(y)
                logval = one(expv)
            elseif iszero(vx) && vy > 0
                logval = zero(vx)
            else
                logval = expv * log(vx)
            end
            new_partials = _mul_partials(partials(x), partials(y), powval, logval)
            return Dual{Txy}(expv, new_partials)
        end,
        begin
            v = value(x)
            expv = ($f)(v, y)
            if y == zero(y)
                new_partials = zero(partials(x))
            else
                new_partials = partials(x) * y * ($f)(v, y - 1)
            end
            return Dual{Tx}(expv, new_partials)
        end,
        begin
            v = value(y)
            expv = ($f)(x, v)
            deriv = (iszero(x) && v > 0) ? zero(expv) : expv * log(x)
            return Dual{Ty}(expv, deriv * partials(y))
        end
        )
    end
end


@inline Base.literal_pow(::typeof(^), x::Dual{T}, ::Val{0}) where {T} =
Dual{T}(one(value(x)), zero(partials(x)))

for y in 1:3
    @eval @inline function Base.literal_pow(::typeof(^), x::Dual{T}, ::Val{$y}) where {T}
        v = value(x)
        expv = v^$y
        deriv = $y * v^$(y - 1)
        return Dual{T}(expv, deriv * partials(x))
    end
end


# hypot #
# -------#

@inline function calc_hypot(x, y, z, ::Type{T}) where T
    vx = value(x)
    vy = value(y)
    vz = value(z)
    h = hypot(vx, vy, vz)
    p = (vx / h) * partials(x) + (vy / h) * partials(y) + (vz / h) * partials(z)
    return Dual{T}(h, p)
end

@define_ternary_dual_op(
Base.hypot,
calc_hypot(x, y, z, Txyz),
calc_hypot(x, y, z, Txy),
calc_hypot(x, y, z, Txz),
calc_hypot(x, y, z, Tyz),
calc_hypot(x, y, z, Tx),
calc_hypot(x, y, z, Ty),
calc_hypot(x, y, z, Tz),
)

# Skipping fma and muladd for now

# sincos #
#--------#

@inline sincos(x) = (sin(x), cos(x))

@inline function sincos(d::Dual{T}) where T
    sd, cd = sincos(value(d))
    return (Dual{T}(sd, cd * partials(d)), Dual{T}(cd, -sd * partials(d)))
end


###################
# Pretty Printing #
###################

function Base.show(io::IO, d::Dual{T,V}) where {T,V}
    print(io, "Dual{$(repr(T))}(", value(d))
    print(io, ", ", partials(d))
    print(io, ")")
end

function Base.typemin(::Type{DynamicForwardDiff.Dual{T,V}}) where {T,V}
    DynamicForwardDiff.Dual{T,V}(typemin(V), zero(Partials{V}))
end

function Base.typemax(::Type{DynamicForwardDiff.Dual{T,V}}) where {T,V}
    DynamicForwardDiff.Dual{T,V}(typemax(V), zero(Partials{V}))
end