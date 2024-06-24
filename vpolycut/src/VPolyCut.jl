using SCIP
using Random
using LinearAlgebra
using Printf

include("constants.jl")
include("VPolySeparator.jl")
include("debug_utils.jl")
include("CornerPolyhedron.jl")
include("BranchAndBound.jl")

function printCut(coef, b)
    str = ""
    for (index, i) in enumerate(coef)
        str = str * " $i" * "x$index"
    end
    str = str * " ≥ $b"
    return str
end

function SCIP.exec_lp(sepa::VPolySeparator)
    # Only limit 1 call 
    if sepa.called >= CALL_LIMIT
        return SCIP.SCIP_DIDNOTRUN
    end
    sepa.called += 1

    # Precondition
    if SCIP.SCIPisLPSolBasic(sepa.scipd) == false
        @warn "Solution is non basic"
        return SCIP.SCIP_DIDNOTRUN
    elseif SCIP.SCIPgetLPSolstat(sepa.scipd) != SCIP.SCIP_LPSOLSTAT_OPTIMAL
        @warn "LP is not optimal"
        return SCIP.SCIP_DIDNOTRUN
    elseif !(SCIP.SCIPgetStage(sepa.scipd) == SCIP.SCIP_STAGE_SOLVING)
        @warn "SCIP not solving"
        return SCIP.SCIP_DIDNOTRUN
    end

    # STEP 1: Get Corner Polyhedron of current LP solution
    lp_sol = nothing
    lp_rays = nothing
    lp_sol, lp_rays = get_corner_polyhedron(sepa.scipd)

    # [CLEAN] For Debug Write LP of current node file THIS STEP MUST BE DONE AFTER getting the corner polyhedron since otherwise LP SOLSTAT will be affected
    name = string(sepa.called)
    SCIP.@SCIP_CALL SCIP.SCIPwriteMIP(sepa.scipd, joinpath(WRITE_PATH, name * ".mip"), true, false, true)

    if DEBUG_PRINT_ORIGINAL_CORNER_POLYHEDRON
        println("====================")
        println("Separator Called")
        println("LP Solution is " * string(lp_sol))
        println("LP Rays are " * string(lp_rays))
        println("====================")
    end

    # [CLEAN] Might actually not need this
    dim = length(lp_sol)
    if length(lp_rays) != dim
        @warn "Not enough ray"
        return SCIP.SCIP_DIDNOTRUN
    end

    #STEP 2: get a v-polyhedral description of the disjunctive hull
    points_collection = []
    rays_collection = []
    fixed = Set()

    SCIP.@SCIP_CALL SCIP.SCIPstartProbing(sepa.scipd)
    get_point_ray_collection(sepa.scipd, points_collection, rays_collection, fixed)
    SCIP.@SCIP_CALL SCIP.SCIPendProbing(sepa.scipd)

    @info "Finish collecting"
    @info points_collection
    @info rays_collection

    #Start Constructing PRLP
    lpi = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPlpiCreate(lpi, C_NULL, "Separating LP", SCIP.SCIP_OBJSEN_MINIMIZE)

    # Get the basis matrix for the basic space
    ray_matrix = stack(lp_rays)

    # Add Variables
    SCIP.@SCIP_CALL SCIP.SCIPlpiAddCols(lpi[], dim, ones(dim), -SCIP.SCIPinfinity(sepa.scipd) * ones(dim), SCIP.SCIPinfinity(sepa.scipd) * ones(dim), C_NULL, 0, C_NULL, C_NULL, C_NULL)
    for point in points_collection
        lhs = [1.0]
        rhs = [SCIP.SCIPinfinity(sepa.scipd)]
        point = point - lp_sol
        point = ray_matrix \ point
        buffer = zeros(dim)
        ind = zeros(Cint, dim)
        nnonz = 0
        j = 1
        for i in 1:dim
            if abs(point[i]) > 10e-6
                nnonz += 1
                buffer[j] = point[i]
                ind[j] = i - 1
                j += 1
            end
        end
        beg = [0]
        SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(lpi[], 1, pointer(lhs), pointer(rhs), C_NULL, nnonz, pointer(beg), pointer(ind), pointer(buffer))
    end
    for ray in rays_collection
        lhs = [0.0]
        rhs = [SCIP.SCIPinfinity(sepa.scipd)]
        ray = ray_matrix \ ray
        @info "Ray" ray
        buffer = zeros(dim)
        ind = zeros(Cint, dim)
        nnonz = 0
        j = 1
        for i in 1:dim
            if abs(ray[i]) > 10e-6
                nnonz += 1
                buffer[j] = ray[i]
                ind[j] = i - 1
                j += 1
            end
        end
        beg = [0]
        SCIP.@SCIP_CALL SCIP.SCIPlpiAddRows(lpi[], 1, pointer(lhs), pointer(rhs), C_NULL, nnonz, pointer(beg), pointer(ind), pointer(buffer))
    end
    SCIP.@SCIP_CALL SCIP.SCIPlpiWriteLP(lpi[], "./temp/Cut$name.lp")
    # Solve The LP
    if SCIP.SCIPlpiHasPrimalSolve() == 0
        @warn "LP Solver Does Not Support Primal Solve"
        return SCIP.SCIP_DIDNOTRUN
    end
    SCIP.@SCIP_CALL SCIP.SCIPlpiSolvePrimal(lpi[])
    if SCIP.SCIPlpiIsPrimalFeasible(lpi[]) == 0
        @warn "LP Solver Failed TO Solve Cut Generation LP"
        return SCIP.SCIP_DIDNOTRUN
    end
    obj = Ref{SCIP.SCIP_Real}(0.0)
    sol = zeros(SCIP.SCIP_Real, dim)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(lpi[], obj, pointer(sol), C_NULL, C_NULL, C_NULL)
    println("CGLP Solution is " * string(sol))
    ray_m = ray_matrix'
    sol = ray_m \ sol
    b = dot(sol, lp_sol) + 1
    @info printCut(sol, b)
    nnonz = 0
    for s in sol
        if abs(s) >= EPSILON
            nnonz += 1
        end
    end
    row = Ref{Ptr{SCIP.SCIP_ROW}}(C_NULL)
    cols = SCIP.LibSCIP.SCIPgetLPCols(sepa.scipd)
    cols = unsafe_wrap(Vector{Ptr{SCIP.LibSCIP.SCIP_COL}}, cols, dim)
    cols__ = Array{Ptr{SCIP.SCIP_ROW}}(undef, nnonz)
    vals__ = Array{SCIP.SCIP_Real}(undef, nnonz)
    infeasible = Ref{SCIP.SCIP_Bool}(0)
    cnt = 1
    for (index, s) in enumerate(sol)
        if abs(s) >= EPSILON
            cols__[cnt] = cols[index]
            vals__[cnt] = s
            cnt += 1
        end
    end

    SCIP.@SCIP_CALL SCIP.SCIPcreateRowSepa(sepa.scipd, row, sepa.scipd.sepas[sepa], "", nnonz, pointer(cols__), pointer(vals__), b, SCIP.SCIPinfinity(sepa.scipd), true, false, false)
    SCIP.@SCIP_CALL SCIP.SCIPaddRow(sepa.scipd, row[], true, infeasible)

    @warn "Row ADDED"

    return SCIP.SCIP_SEPARATED
end
