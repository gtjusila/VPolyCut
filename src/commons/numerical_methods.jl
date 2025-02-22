# wrapper around SCIP Numerics Library To Prevent ugly SCIP.SCIPisInfinity(scip, x) == 1 calls
using SCIP

const global numeric_scip = Ref{Union{SCIP.SCIPData,Nothing}}(nothing)
const global scip_epsilon = Ref{SCIP.SCIP_Real}(0.0)
function is_numeric_scip_set()
    return !isnothing(numeric_scip[])
end
function check_numeric_scip()
    if !is_numeric_scip_set()
        throw(
            "numeric_scip is not set. before using the wrapped numerics function set numeric_scip by calling set_numeric_scip(scip)"
        )
    end
end

function set_numeric_scip(scip::SCIP.SCIPData)
    global numeric_scip[] = scip
    global scip_epsilon
    scip_epsilon[] = SCIP.SCIPepsilon(scip)
end

function is_infinity(x::SCIP.SCIP_Real)
    return SCIP.SCIPisInfinity(numeric_scip[], x) == 1
end

function is_zero(x::SCIP.SCIP_Real)
    global scip_epsilon
    return abs(x) <= scip_epsilon[]
end

function is_non_zero(x::SCIP.SCIP_Real)
    return SCIP.SCIPisZero(numeric_scip[], x) == 0
end

function is_positive(x::SCIP.SCIP_Real)
    return SCIP.SCIPisPositive(numeric_scipp[], x) == 1
end

function is_negative(x::SCIP.SCIP_Real)
    return SCIP.SCIPisNegative(numeric_scip[], x) == 1
end

function is_EQ(x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisEQ(numeric_scip[], x, y) == 1
end

function is_GT(x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisGT(numeric_scip[], x, y) == 1
end

function is_GE(x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisGE(numeric_scip[], x, y) == 1
end

function is_LT(x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisLT(numeric_scip[], x, y) == 1
end

function is_LE(x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisLE(numeric_scip[], x, y) == 1
end

function is_integral(x::SCIP.SCIP_Real)
    return SCIP.SCIPisFeasIntegral(numeric_scip[], x) == 1
end

function is_true(i)
    return i == 1
end

function is_false(i)
    return i == 0
end
