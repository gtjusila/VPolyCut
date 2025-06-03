function estimate_prlp_nonzero(sepa::VPCSeparator)
    shared = sepa.shared_data

    estimate::Int64 = 0
    corner_polyhedron = CornerPolyhedron(shared.scipd, shared.nonbasic_space)
    estimate = SparseArrays.nnz(corner_polyhedron.lp_sol)
    @debug "Estimated number of nonzeros in PRLP: $estimate"
    for ray in get_lp_rays(corner_polyhedron)
        estimate += SparseArrays.nnz(ray)
    end
    @debug "Estimated number of nonzeros in PRLP after rays: $estimate"
    estimate = estimate * sepa.parameters.n_leaves
    @debug "Estimated number of nonzeros in PRLP after leaves: $estimate"
    sepa.statistics.prlp_nnonzeros_estimate = estimate
end