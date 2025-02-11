using SCIP

"""
    construct_prlp(point_ray_collection::PointRayCollection, non_basic_space::NonBasicSpace)

Construct a PRLP object from a given point ray collection
"""
function construct_prlp(
    point_ray_collection::PointRayCollection
)
    @debug "Constructing PRLP"
    problem_dimension = dimension(point_ray_collection)
    prlp = PRLP(problem_dimension)

    for point in get_points(point_ray_collection)
        PRLPaddPoint(prlp, get_point(point))
    end
    for ray in get_rays(point_ray_collection)
        PRLPaddRay(prlp, ray)
    end

    PRLPconstructLP(prlp)
    @debug "PRLP Constructed"

    return prlp
end