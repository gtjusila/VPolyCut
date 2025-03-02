using SCIP

"""
    construct_prlp(point_ray_collection::PointRayCollection, non_basic_space::NonBasicSpace)

Construct a PRLP object from a given point ray collection
"""
function construct_prlp(
    sepa::VPCSeparator
)
    scip = sepa.shared_data.scipd
    point_ray_collection = sepa.shared_data.point_ray_collection
    prlp_solve_method = PRLPsolveAlgorithm(sepa.parameters.prlp_solve_method)
    prlp_allow_warm_start = sepa.parameters.prlp_allow_warm_start

    problem_dimension = dimension(point_ray_collection)

    if sepa.parameters.apply_beta_scaling
        if is_GE(sepa.shared_data.lp_obj_nonbasic, 1.0)
            beta = sepa.shared_data.lp_obj_nonbasic
        end
    end

    prlp = PRLP(;
        dimension = problem_dimension,
        scip = scip,
        allow_warm_start = prlp_allow_warm_start,
        beta = beta)

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