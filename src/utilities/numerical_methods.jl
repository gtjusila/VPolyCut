# wrapper around SCIP Numerics Library To Prevent ugly SCIP.SCIPisInfinity(scip, x) == 1 calls
using SCIP

function is_infinity(scip::SCIP.SCIPData, x::SCIP.SCIP_Real)
    return SCIP.SCIPisInfinity(scip, x) == 1
end

function is_zero(scip::SCIP.SCIPData, x::SCIP.SCIP_Real)
    return SCIP.SCIPisZero(scip, x) == 1
end

function is_non_zero(scip::SCIP.SCIPData, x::SCIP.SCIP_Real)
    return SCIP.SCIPisZero(scip, x) == 0
end

function is_positive(scip::SCIP.SCIPData, x::SCIP.SCIP_Real)
    return SCIP.SCIPisPositive(scip, x) == 1
end

function is_negative(scip::SCIP.SCIPData, x::SCIP.SCIP_Real)
    return SCIP.SCIPisNegative(scip, x) == 1
end

function is_EQ(scip::SCIP.SCIPData, x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisEQ(scip, x, y) == 1
end

function is_GT(scip::SCIP.SCIPData, x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisGT(scip, x, y) == 1
end

function is_GE(scip::SCIP.SCIPData, x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisGE(scip, x, y) == 1
end

function is_LT(scip::SCIP.SCIPData, x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisLT(scip, x, y) == 1
end

function is_LE(scip::SCIP.SCIPData, x::SCIP.SCIP_Real, y::SCIP.SCIP_Real)
    return SCIP.SCIPisLE(scip, x, y) == 1
end

function is_integral(scip::SCIP.SCIPData, x::SCIP.SCIP_Real)
    return SCIP.SCIPisFeasIntegral(scip, x) == 1
end
