#
# check_lpi_conventions.jl
#
# A small script with a bunch of asserts to ensure that the LPI used
# by SCIP conforms to the standard mentioned in https://scipopt.org/doc/html/group__LPIS.php
#
using JuMP
using VPolyhedralCut
using VPolyhedralCut.SCIPJLUtils
using SCIP
using Printf
using LinearAlgebra

### Global variables ###
const counter = Ref(1)

### Helper functions ###
function printmatrix(matrix::T, title) where {T<:AbstractMatrix}
    println(title)
    nrows, ncolumns = size(matrix)
    for i in 1:nrows
        for j in 1:ncolumns
            print(@sprintf("%5.2f ", matrix[i, j]))
        end
        println()
    end
end

function isequalwithinscipepsilon(
    scip::SCIP.SCIPData, matrix1::T, matrix2::T
) where {T<:AbstractMatrix}
    if size(matrix1) != size(matrix2)
        return false
    end
    nrows, ncolumns = size(matrix1)
    for i in 1:nrows
        for j in 1:ncolumns
            if (SCIP.SCIPisEQ(scip, matrix1[i, j], matrix2[i, j]) == 0)
                @info "Differ On Entry $(i) $(j)"
                return false
            end
        end
    end
    return true
end

### Test Function ###

function run_test_with_model(model)
    println("Running Test $(counter[])")
    println(model)
    counter[] = counter[] + 1
    data = Dict(
        "scip" => get_scip_data_from_model(model)
    )

    # Turn off display and presolving
    set_attribute(model, "display/verblevel", 0)
    set_attribute(model, "presolving/maxrounds", 0)

    lp_test_function = function (data)
        scip = data["scip"]
        nvars = SCIP.SCIPgetNLPCols(scip)
        nrows = SCIP.SCIPgetNLPRows(scip)
        constraint_matrix = zeros(SCIP.SCIP_Real, nrows, nvars + nrows)

        lp_cols = unsafe_wrap(Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIPgetLPCols(scip), nvars)
        lp_rows = unsafe_wrap(Vector{Ptr{SCIP.SCIP_ROW}}, SCIP.SCIPgetLPRows(scip), nrows)

        ###  Get constraint matrix ###
        for i in 1:nrows
            row = lp_rows[i]
            nnonz = SCIP.SCIProwGetNNonz(row)
            nonzero_cols = unsafe_wrap(
                Vector{Ptr{SCIP.SCIP_COL}}, SCIP.SCIProwGetCols(row), nnonz
            )
            nonzero_vals = unsafe_wrap(
                Vector{SCIP.SCIP_Real}, SCIP.SCIProwGetVals(row), nnonz
            )

            for j in 1:nnonz
                for k in 1:nvars
                    if lp_cols[k] == nonzero_cols[j]
                        constraint_matrix[i, k] = nonzero_vals[j]
                    end
                end
            end

            constraint_matrix[i, nvars + i] = 1 # Slack Variable
        end

        printmatrix(constraint_matrix, "Constraint Matrix")

        ### Get basis indices ###

        basis_indices = zeros(Cint, nrows)
        SCIP.@SCIP_CALL SCIP.SCIPgetLPBasisInd(scip, pointer(basis_indices))

        # Convert to julia indexing and rows to slack variable
        for i in 1:nrows
            # C to Julia index (+1)
            if basis_indices[i] >= 0
                basis_indices[i] += 1
            else
                basis_indices[i] = nvars - basis_indices[i]
            end
        end

        println("Basis Indices")
        println(basis_indices)

        ### Get Basis Matrix ###
        basis_matrix = constraint_matrix[:, basis_indices]
        printmatrix(basis_matrix, "Basis Matrix")

        ### Get Basis Inverse By Row ###
        basis_inverse_by_row = zeros(SCIP.SCIP_Real, nrows, nrows)
        for i in 1:nrows
            buffer = zeros(SCIP.SCIP_Real, nrows)
            SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvRow(
                scip, i - 1, pointer(buffer), C_NULL, C_NULL
            )
            basis_inverse_by_row[i, :] = buffer
        end
        printmatrix(basis_inverse_by_row, "Basis Inverse By Row")

        ### Basis Inverse By Row * Basis  (Should be I) ###
        basis_inverse_row_times_basis = basis_inverse_by_row * basis_matrix
        printmatrix(basis_inverse_row_times_basis, "Basis Inverse By Row * Basis Matrix")

        @assert (isequalwithinscipepsilon(
            scip, basis_inverse_row_times_basis, Matrix{Float64}(I, nrows, nrows)
        ))
        ### Get Basis Inverse By Column ###
        basis_inverse_by_col = zeros(SCIP.SCIP_Real, nrows, nrows)
        for i in 1:nrows
            buffer = zeros(SCIP.SCIP_Real, nrows)
            SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvCol(
                scip, i - 1, pointer(buffer), C_NULL, C_NULL
            )
            basis_inverse_by_col[:, i] = buffer
        end
        printmatrix(basis_inverse_by_col, "Basis Inverse By Col")

        ### Basis Inverse By Column * Basis ###
        basis_inverse_by_col_times_basis = basis_inverse_by_col * basis_matrix
        printmatrix(basis_inverse_by_col_times_basis, "Basis Inverse By Col * Basis Matrix")

        @assert (isequalwithinscipepsilon(
            scip, basis_inverse_by_col_times_basis, Matrix{Float64}(I, nrows, nrows)
        ))

        basis_inverse_times_constraint = basis_inverse_by_col * constraint_matrix
        printmatrix(
            basis_inverse_times_constraint,
            "Basis inverse times constraint (solution, using basis inverse by col)"
        )

        ### Basis Inverse times Constraint by column ###
        tableau_by_column = zeros(SCIP.SCIP_Real, nrows, nrows + nvars)
        for i in 1:nvars
            buffer = zeros(SCIP.SCIP_Real, nrows)
            SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvACol(scip, i - 1, buffer, C_NULL, C_NULL)
            tableau_by_column[:, i] = buffer
        end
        tableau_by_column[:, (nvars + 1):end] = basis_inverse_by_col
        printmatrix(tableau_by_column, "Basis inverse times constraint (by column)")
        @assert (isequalwithinscipepsilon(
            scip, basis_inverse_times_constraint, tableau_by_column)
)

        ### Basis Inverse times Constraint by row ###
        tableau_by_row = zeros(SCIP.SCIP_Real, nrows, nrows + nvars)
        for i in 1:nrows
            buffer = zeros(SCIP.SCIP_Real, nrows + nvars)
            SCIP.@SCIP_CALL SCIP.SCIPgetLPBInvARow(
                scip, i - 1, C_NULL, buffer, C_NULL, C_NULL
            )
            tableau_by_row[i, :] = buffer
        end
        tableau_by_row[:, (nvars + 1):end] = basis_inverse_by_col
        printmatrix(tableau_by_row, "Basis inverse times constraint (by row)")
        @assert (isequalwithinscipepsilon(
            scip, basis_inverse_times_constraint, tableau_by_row
        ))
    end

    lp_test = initiate_callback(
        lp_test_function,
        model,
        data
    )

    SCIP.@SCIP_CALL SCIP.SCIPtransformProb(data["scip"])

    register_callback(
        model,
        SCIP.SCIP_EVENTTYPE_FIRSTLPSOLVED,
        lp_test
    )

    optimize!(model)
end

# Taken from  Example from Chapter 2 of the book Vanderbei - Linear Programming: Foundations and Extensions

### Test case 1 ###
model = direct_model(SCIP.Optimizer())

@variable(model, x[1:3] >= 0);
@objective(model, Max, 5 * x[1] + 4 * x[2] + 3 * x[3]);
@constraint(model, c1, 2 * x[1] + x[2] + x[3] <= 5);
@constraint(model, c2, 4 * x[1] + x[2] + 2 * x[3] <= 11);
@constraint(model, c3, 3 * x[1] + 4 * x[2] + 2 * x[3] <= 8);

run_test_with_model(model)

### Test case 2 ###
model = direct_model(SCIP.Optimizer())

@variable(model, x[1:3] >= 0);
@objective(model, Min, -5 * x[1] - 4 * x[2] - 3 * x[3]);
@constraint(model, c1, 2 * x[1] + x[2] + x[3] <= 5);
@constraint(model, c2, 4 * x[1] + x[2] + 2 * x[3] <= 11);
@constraint(model, c3, 3 * x[1] + 4 * x[2] + 2 * x[3] <= 8);

run_test_with_model(model)

### Test case 3 ###
model = direct_model(SCIP.Optimizer())

@variable(model, x[1:3] >= 0);
@objective(model, Max, 5 * x[1] + 4 * x[2] + 3 * x[3]);
@constraint(model, c1, -2 * x[1] - x[2] - x[3] >= -5);
@constraint(model, c2, -4 * x[1] - x[2] - 2 * x[3] >= -11);
@constraint(model, c3, -3 * x[1] - 4 * x[2] - 2 * x[3] >= -8);

run_test_with_model(model)

### Test case 4 ###
model = direct_model(SCIP.Optimizer())

@variable(model, x[1:3] >= 0);
@objective(model, Min, -5 * x[1] - 4 * x[2] - 3 * x[3]);
@constraint(model, c1, -2 * x[1] - x[2] - x[3] >= -5);
@constraint(model, c2, -4 * x[1] - x[2] - 2 * x[3] >= -11);
@constraint(model, c3, -3 * x[1] - 4 * x[2] - 2 * x[3] >= -8);

run_test_with_model(model)

### Test case 5 ###
model = direct_model(SCIP.Optimizer())

@variable(model, x[1:3] >= 0);
@objective(model, Min, -5 * x[1] - 4 * x[2] - 3 * x[3]);
@constraint(model, c1, -2 * x[1] - x[2] - x[3] >= -5);
@constraint(model, c2, -4 * x[1] - x[2] - 2 * x[3] >= -11);
@constraint(model, c3, 3 * x[1] + 4 * x[2] + 2 * x[3] <= 8);

run_test_with_model(model)

### Test case 6 ###
model = direct_model(SCIP.Optimizer())

@variable(model, x[1:3] >= 0);
@objective(model, Max, 5 * x[1] + 4 * x[2] + 3 * x[3]);
@constraint(model, c1, -2 * x[1] - x[2] - x[3] >= -5);
@constraint(model, c2, -4 * x[1] - x[2] - 2 * x[3] >= -11);
@constraint(model, c3, 3 * x[1] + 4 * x[2] + 2 * x[3] <= 8);

run_test_with_model(model)