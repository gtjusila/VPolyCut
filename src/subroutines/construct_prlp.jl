using SCIP

"""
    construct_prlp(point_ray_collection::PointRayCollection, non_basic_space::NonBasicSpace)

Construct a PRLP object from a given point ray collection
"""
function construct_prlp(
    point_ray_collection::PointRayCollection, non_basic_space::NonBasicSpace;
    scip::SCIP.SCIPData = SCIP.Optimizer().inner
)
    # Step 1: Project points and ray
    projected_points = [
        project_point_to_nonbasic_space(non_basic_space, get_point(point))
        for point in get_points(point_ray_collection)
    ]
    projected_rays = [
        project_ray_to_nonbasic_space(non_basic_space, ray)
        for ray in get_rays(point_ray_collection)
    ]

    #Step 2: Apply zeroing step
    points_zeroed = map(projected_points) do point
        return [!is_zero(scip, p) ? p : 0.0 for p in point]
    end

    rays_zeroed = map(projected_rays) do ray
        return [!is_zero(scip, p) ? p : 0.0 for p in get_coefficients(ray)]
    end

    #Step 3: Construct PRLP
    @debug "Constructing PRLP"
    problem_dimension = length(points_zeroed[1])
    prlp = PRLP(problem_dimension)

    for point in points_zeroed
        PRLPaddPoint(prlp, point)
    end
    for ray in rays_zeroed
        PRLPaddRay(prlp, ray)
    end

    PRLPconstructLP(prlp)
    @debug "PRLP Constructed"

    return prlp
end