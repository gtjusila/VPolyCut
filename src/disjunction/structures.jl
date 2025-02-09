#
# structures.jl
# Contains the data structures that are needed to represent disjunctions
# along with some getter and setter functions to modify them
#

using SCIP

"""
    DisjunctiveTerm

A disjunctive term is a set of local bounds on the variables. It have two member variables `lower_bounds` and `upper_bounds` which are dictionaries of variable pointers to the lower and upper bounds respectively. 
"""
@kwdef struct DisjunctiveTerm
    lower_bounds::Dict{Ptr{SCIP.SCIP_VAR},SCIP.SCIP_Real} = Dict()
    upper_bounds::Dict{Ptr{SCIP.SCIP_VAR},SCIP.SCIP_Real} = Dict()
end

"""
    add_lower_bound!(term::DisjunctiveTerm, var::Ptr{SCIP.SCIP_VAR}, lower::SCIP.SCIP_Real)

Add a lower bound to the disjunctive term
"""
function set_lower_bound!(
    term::DisjunctiveTerm,
    var::Ptr{SCIP.SCIP_VAR},
    lower::SCIP.SCIP_Real
)
    term.lower_bounds[var] = lower
end

"""
    add_upper_bound!(term::DisjunctiveTerm, var::Ptr{SCIP.SCIP_VAR}, upper::SCIP.SCIP_Real)

Add an upper bound to the disjunctive term
"""
function set_upper_bound!(
    term::DisjunctiveTerm,
    var::Ptr{SCIP.SCIP_VAR},
    upper::SCIP.SCIP_Real
)
    term.upper_bounds[var] = upper
end

"""
    add_bound!(term::DisjunctiveTerm, var::Ptr{SCIP.SCIP_VAR}, lower::SCIP.SCIP_Real, upper::SCIP.SCIP_Real)

Add both lower and upper bounds to the disjunctive term
"""
function add_bound!(
    term::DisjunctiveTerm,
    var::Ptr{SCIP.SCIP_VAR},
    lower::SCIP.SCIP_Real,
    upper::SCIP.SCIP_Real
)
    set_lower_bound!(term, var, lower)
    add_upper_bound!(term, var, upper)
end

@enum DisjunctiveTermBoundStatus UNBOUNDED LOWER_BOUNDED UPPER_BOUNED BOUNDED

"""
    get_disjunctive_term_bound_status(term::DisjunctiveTerm, var::Ptr{SCIP.SCIP_VAR})

Get the bound status of the variable in the disjunctive term
"""
function get_disjunctive_term_bound_status(
    term::DisjunctiveTerm,
    var::Ptr{SCIP.SCIP_VAR}
)
    if haskey(term.lower_bounds, var) && haskey(term.upper_bounds, var)
        return BOUNDED
    elseif haskey(term.lower_bounds, var)
        return LOWER_BOUNDED
    elseif haskey(term.upper_bounds, var)
        return UPPER_BOUNDED
    else
        return UNBOUNDED
    end
end

"""
    get_lower_bound(term::DisjunctiveTerm, var::Ptr{SCIP.SCIP_VAR})

Get the lower bound of the variable in the disjunctive term. Defaults to typemin(Int64) if the variable is not bounded
"""
function get_lower_bound(
    term::DisjunctiveTerm,
    var::Ptr{SCIP.SCIP_VAR}
)
    if haskey(term.lower_bounds, var)
        return term.lower_bounds[var]
    end
    return typemin(Int64)
end

"""
    get_upper_bound(term::DisjunctiveTerm, var::Ptr{SCIP.SCIP_VAR})

Get the upper bound of the variable in the disjunctive term. Default to typemax(Int64) if the variable is not bounded
"""
function get_upper_bound(
    term::DisjunctiveTerm,
    var::Ptr{SCIP.SCIP_VAR}
)
    if haskey(term.upper_bounds, var)
        return term.upper_bounds[var]
    end
    return typemax(Int64)
end

"""
    get_bounds(term::DisjunctiveTerm, var::Ptr{SCIP.SCIP_VAR})

Get the lower and upper bounds of the variable in the disjunctive term.
"""
function get_bounds(
    term::DisjunctiveTerm,
    var::Ptr{SCIP.SCIP_VAR}
)
    return get_lower_bound(term, var), get_upper_bound(term, var)
end

"""
    apply_bound_changes(term::DisjunctiveTerm, scip::SCIP.SCIPData)

scip should be in probing mode, this function will apply the bounds in the disjunctive term to the scip data
"""
function apply_bound_changes(term::DisjunctiveTerm, scip::SCIP.SCIPData)
    @assert is_true(SCIP.SCIPinProbing(scip)) "SCIP should be in probing mode"
    for (var, lower) in term.lower_bounds
        SCIP.SCIPchgVarLbProbing(scip, var, lower)
    end
    for (var, upper) in term.upper_bounds
        SCIP.SCIPchgVarUbProbing(scip, var, upper)
    end
end

"""
    Disjunction

A disjunction is a set of disjunctive terms. It have one member variable `terms` which is a vector of disjunctive terms.
It is iterable so if `d` is a disjunction then `for term in d` will iterate over the disjunctive terms.
"""
@kwdef struct Disjunction
    terms::Vector{DisjunctiveTerm} = []
end

function Base.iterate(d::Disjunction, state=1)
    return state > length(d.terms) ? nothing : (d.terms[state], state + 1)
end

function Base.length(d::Disjunction)
    return length(d.terms)
end

function add_term!(d::Disjunction, term::DisjunctiveTerm)
    push!(d.terms, term)
end