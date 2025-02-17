using SCIP

"""
    construct_prlp(point_ray_collection::PointRayCollection, non_basic_space::NonBasicSpace)

Construct a PRLP object from a given point ray collection
"""
function construct_prlp(
    scip::SCIP.SCIPData,
    point_ray_collection::PointRayCollection;
    prlp_solve_method::PRLPsolveAlgorithm = PRIMAL_SIMPLEX,
    prlp_allow_warm_start::Bool = true
)
    problem_dimension = dimension(point_ray_collection)
    prlp = PRLP(;
        dimension = problem_dimension,
        scip = scip,
        allow_warm_start = prlp_allow_warm_start)

    PRLPsetSolvingAlgorithm(prlp, prlp_solve_method)

    first = true
    for point in get_points(point_ray_collection)
        PRLPaddPoint(prlp, point)
    end
    for ray in get_rays(point_ray_collection)
        PRLPaddRay(prlp, ray)
    end

    PRLPconstructLP(prlp)
    @debug "PRLP Constructed"

    return prlp
end