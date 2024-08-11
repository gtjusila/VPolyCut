using SCIP
using LoggingExtras
using JuMP
using MathOptInterface
using Lazy

mutable struct Ray <: AbstractVector{SCIP.SCIP_Real}
    coefficients::Vector{SCIP.SCIP_Real}
    generating_variable::Variable
end

function get_coefficients(ray::Ray)
    return ray.coefficients
end

function get_generating_variable(ray::Ray)
    return ray.generating_variable
end

function set_coefficients!(ray::Ray, coefficients::Vector{SCIP.SCIP_Real})
    ray.coefficients = coefficients
end

@forward Ray.coefficients Base.size, Base.getindex, Base.setindex!

const Point = Vector{SCIP.SCIP_Real}

function get_logfolder_path()
    logfolder_path = joinpath(dirname(Base.active_project()), "logs")
    if !isdir(logfolder_path)
        mkdir(logfolder_path)
    end
    return logfolder_path
end

function setup_file_logger(filename::String)
    return FormatLogger(open(filename, "w")) do io, args
        # Write the module, level and message only
        println(io, "[", args.level, "] ", args.message)
    end
end

# A Hack to allow calls to SCIP for julia model
function create_subscip_model()
    inner = MOI.Bridges.full_bridge_optimizer(SCIP.Optimizer(), Float64)
    model = direct_generic_model(Float64, inner)
    set_attribute(model, "display/verblevel", 0)
    set_attribute(model, "presolving/maxrounds", 0)
    backend = unsafe_backend(model)
    scip = backend.inner
    SCIP.@SCIP_CALL SCIP.SCIPsetSubscipsOff(scip, SCIP.SCIP_Bool(true))
    return model
end

# Add Row alpha * x <= beta to scip 
function add_sepa_row!(
    scip::SCIP.SCIPData,
    sepa::T,
    alpha::Vector{Float64},
    xs::Vector{Ptr{SCIP.SCIP_Var}},
    beta::Float64;
    valid_globaly::Bool=true,
    modifiable::Bool=false,
    removable::Bool=false
) where {T<:SCIP.AbstractSeparator}
    @assert length(alpha) == length(xs)

    new_row = Ref{Ptr{SCIP.SCIP_ROW}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPcreateEmptyRowSepa(
        scip,
        new_row,
        scip.sepas[sepa],
        "",
        -SCIP.SCIPinfinity(scip),
        beta,
        !valid_globaly,
        modifiable,
        removable
    )

    infeasible = Ref{SCIP.SCIP_Bool}(0)
    for (var, sol) in zip(xs, alpha)
        @assert !is_infinity(scip, abs(sol))
        if is_non_zero(scip, sol)
            SCIP.@SCIP_CALL SCIP.SCIPaddVarToRow(scip, new_row[], var, sol)
        end
    end

    #SCIP.@SCIP_CALL SCIP.SCIPprintRow(scip, new_row[], C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPaddRow(scip, new_row[], true, infeasible)
    SCIP.@SCIP_CALL SCIP.SCIPreleaseRow(scip, new_row)
end