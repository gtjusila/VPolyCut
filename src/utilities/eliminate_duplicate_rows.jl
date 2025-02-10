# Eliminate Duplicate rows
# Based on Tomlin, L. A., and James S. Welch. "Finding duplicate rows in a linear programming model." Operations Research Letters 5, no. 1 (1986): 7-11.
# Two rows are duplicate if they are positive scalar multiples of each other
# One can implement this more generally but we only need the positive scalar case.
using SCIP

const INF = typemax(Int) - 1

function get_unique_row_indices(
    original_matrix::Matrix{SCIP.SCIP_Real},
    scip::SCIP.SCIPData)::Vector{Int}

    # Original matrix is of size n x m
    n = size(original_matrix, 1)
    m = size(original_matrix, 2)

    f = zeros(Float64, n) # Row Factors
    set = zeros(Int, n) # an integer array indicating the set number to which a row belong

    # 5.1 Row Pass
    next_set = 1
    v = n # The number of rows not in the infinity set. All rows are initially considered

    # 5.2 Matrix Pass
    for j in 1:m
        # Iterate through each column j
        n_new_rows = 0
        current_set = next_set
        next_set = next_set + 1
        r0 = -1

        # Iterate through each row r with a_{rj} != 0
        non_zero_indices = findall(z -> !is_zero(z), original_matrix[:, j])
        for (k, r) in enumerate(non_zero_indices)
            # Divide into cases based on set[r]
            if set[r] == 0
                # 5.2.a
                r0 = r
                f[r] = original_matrix[r, j]
                set[r] = current_set
                n_new_rows += 1

            elseif set[r] < current_set
                # 5.2.b
                satisfied = false # Condition (1)

                # 5.2.b.i 
                for l in (k + 1):length(non_zero_indices)
                    i = non_zero_indices[l]
                    if set[i] == set[r] &&
                        is_EQ(
                        scip,
                        (original_matrix[r, j] * f[i]) / (f[r] * original_matrix[i, j]),
                        1.0
                    )
                        set[r] = next_set
                        set[i] = next_set
                        satisfied = true
                    end
                end # for i>r

                # 5.2.b.ii
                if satisfied
                    next_set = next_set + 1
                else
                    set[r] = INF
                    v = v - 1
                    if v == 0
                        @goto TERMINATE_ALGORITHM
                    end
                end

            elseif set[r] == INF || set[r] > current_set
                # Nothing to do here you can ignore this row safely
            else
                error("You should not be here. r is $r and set[r] is $(set[r])")
            end # if set[r] == 0 or set[r] < current_set
        end # for r with a_rj != 0

        # 5.2.c
        if n_new_rows == 1
            set[r0] = INF
            v = v - 1
            if v == 0
                @goto TERMINATE_ALGORITHM
            end
        end

        if n <= 1
            # current_set is empty
        end
    end # for each nonfixed column j

    @label TERMINATE_ALGORITHM
    unique_row_indices = []
    for r in 1:n
        if set[r] == 0
            # Row is empty 5.3.a
        elseif set[r] == INF
            # Row is unique 5.3.b
            push!(unique_row_indices, r)
        elseif set[r] > 0
            # Row belong to the set[r]-th set of duplicates
            # Row is the first row in the set
            for i in (r + 1):n
                # Go over all other elements of set[r]
                if set[i] == set[r]
                    if is_positive(scip, f[i]) == is_positive(scip, f[r])
                        # f[i] and f[r] have the same sign
                        set[i] = -1
                    end
                end
            end

            # Push row as the prototypicall row of type set[r]
            push!(unique_row_indices, r)
        elseif set[r] == -1
            # Row has been marked as duplicate
        else
            error("You should not be here")
        end
    end
    return unique_row_indices
end

function get_unique_row_indices(
    set_of_rows::Vector{Vector{Float64}},
    scip::SCIP.SCIPData
)::Vector{Int}
    matrix = reduce(vcat, transpose.(set_of_rows))
    unique_rows = get_unique_row_indices(matrix, scip)
    return unique_rows
end