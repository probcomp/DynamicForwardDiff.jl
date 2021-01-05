# DynamicForwardDiff.jl

Much of this code is copied (sometimes with slight modifications) from [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl). 
The modifications are designed to enable differentiation with respect to a dynamically
growing set of input variables.
For example, a probabilistic programming library might run a stochastic computation and compute derivatives of the output with respect to the random numbers sampled during execution (whose values, or even number, cannot be known ahead of time).

The package may also be useful for computing Jacobians when each output depends only on a small number of inputs: we use sparse memory representation for dual numbers, i.e., intermediate values do not explicitly store partial derivatives that are equal to 0.

## Usage

Begin by creating a DiffConfig:

```julia
cfg = DynamicForwardDiff.DiffConfig()
```

Then, at any point during a computation, call `DynamicForwardDiff.new_dual(cfg, v)` on a value `v` of type `<: Real` or `AbstractArray{<: Real}` to begin tracking it as an input. (In the case of an `AbstractArray`, each scalar entry will be registered as separate input.) It is assumed that for two scalar inputs `u` and `v` registered during the course of a computation, `du/dv = 0`.

The final result of the computation will be a value of type `DynamicForwardDiff.Dual`. You can use `value(result)` to extract its value, or `partials(result)` to extract the gradient of `result` with respect to all inputs added by `cfg`, in the order they were added.

For example:

```julia
cfg = DynamicForwardDiff.DiffConfig()

# Stochastic computation with either one or two random choices
x = DynamicForwardDiff.new_dual(cfg, rand())
if x > 0.5
  x = x^2
  y = DynamicForwardDiff.new_dual(cfg, rand())
else
  y = 7
end
result = x*y

# Gradient of `result` with respect to the random choices
# (will be either length 2 or length 1).
println(DynamicForwardDiff.partials(result))
```

