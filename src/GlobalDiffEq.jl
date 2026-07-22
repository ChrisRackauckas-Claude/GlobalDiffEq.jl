module GlobalDiffEq

import OrdinaryDiffEq, Richardson, SciMLBase
import SciMLBase: __solve
using PrecompileTools: @setup_workload, @compile_workload
using SciMLBase: AbstractDAEProblem, AbstractODEAlgorithm, AbstractODEProblem, ODEProblem,
    solve

abstract type GlobalDiffEqAlgorithm <: AbstractODEAlgorithm end

"""
    GlobalRichardson(alg)

Wrap the fixed-step ODE algorithm `alg` with global Richardson extrapolation.

`GlobalRichardson` solves the problem at successively refined fixed step sizes and uses
Richardson extrapolation to improve the solution at the requested output times. The wrapped
algorithm must be an `AbstractODEAlgorithm`; adaptivity is disabled for the inner solves.

# Arguments

- `alg`: An ODE algorithm to run at refined fixed step sizes, such as `Tsit5()` or `SSPRK33()`.

# Fields

- `alg`: The wrapped fixed-step ODE algorithm.

# Example

```jldoctest
julia> using GlobalDiffEq, OrdinaryDiffEq, SciMLBase

julia> prob = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0));

julia> sol = solve(prob, GlobalRichardson(Tsit5()); dt = 0.1);

julia> length(sol.u) > 1
true
```
"""
struct GlobalRichardson{A <: AbstractODEAlgorithm} <: GlobalDiffEqAlgorithm
    alg::A
end

# Forward algorithm traits to the wrapped algorithm
# This allows GlobalRichardson to inherit capabilities from the inner algorithm
SciMLBase.allows_arbitrary_number_types(alg::GlobalRichardson) =
    SciMLBase.allows_arbitrary_number_types(alg.alg)
SciMLBase.allowscomplex(alg::GlobalRichardson) =
    SciMLBase.allowscomplex(alg.alg)
SciMLBase.isautodifferentiable(alg::GlobalRichardson) =
    SciMLBase.isautodifferentiable(alg.alg)

function __solve(
        prob::Union{AbstractODEProblem, AbstractDAEProblem},
        alg::GlobalRichardson, args...;
        dt, kwargs...
    )
    opt = Dict(kwargs)
    otheropts = delete!(copy(opt), :dt)
    tstops = get(opt, :tstops, range(prob.tspan[1], stop = prob.tspan[2], step = dt))
    local sol
    val,
        err = Richardson.extrapolate(
        dt, rtol = get(opt, :reltol, 1.0e-3),
        atol = get(opt, :abstol, 1.0e-6), contract = 0.5
    ) do _dt
        sol = solve(prob, alg.alg, args...; dt = _dt, adaptive = false, otheropts...)
        # Convert Vector{Vector{T}} to Matrix{T} for Richardson.jl compatibility
        reduce(hcat, sol.(tstops))
    end
    return sol
end

export GlobalRichardson

@setup_workload begin
    # Simple test ODE: exponential decay du/dt = -u
    function f!(du, u, p, t)
        du[1] = -u[1]
    end
    u0 = [1.0]
    tspan = (0.0, 1.0)
    prob = ODEProblem(f!, u0, tspan)

    @compile_workload begin
        solve(
            prob, GlobalRichardson(OrdinaryDiffEq.Tsit5()),
            dt = 0.1, reltol = 1.0e-3, abstol = 1.0e-6
        )
    end
end

end
