using SCIP
"""
DenseTableau

A dense tableau data structure
"""
struct DenseTableau
    tableau_data::Matrix{SCIP.SCIP_Real}

end

"""
size(obj::DenseTableau)
"""
function size(obj::DenseTableau)
    return size(obj.tableau_data)
end

const Tableau = Matrix{SCIP.SCIP_Real}

function create_tableau(row_count::UInt, column_count::UInt)
    return Tableau(undef, row_count, column_count)
end

function get_tableau_entry(tableau::Tableau, row::UInt, column::UInt)::SCIP.SCIP_Real
    return tableau[row, column]
end

function set_tableau_entry(tableau::Tableau, row::UInt, column::UInt, value::SCIP.SCIP_Real)::SCIP.SCIP_Real
    tableau[row, column] = value
end

function copy_vector_to_tableau_row(tableau::Tableau, tableau_row::UInt, vector::Vector)
    tableau[tableau_row, :] = vector
end