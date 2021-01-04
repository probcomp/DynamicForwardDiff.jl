module DynamicForwardDiff

using DiffRules, DiffResults
using DiffResults: DiffResult, MutableDiffResult, ImmutableDiffResult
using StaticArrays
using Random
import SpecialFunctions
import CommonSubexpressions

include("partials.jl")
include("dual.jl")
include("apiutils.jl")

end # module
