using SCIP

"""
check_time_limit_exceeded_and_raise

Check if time limit is exceeded and raise an exception if it is
"""
function check_time_limit_exceeded_and_raise(scip::SCIP.SCIPData)
    time_limit = Ref{SCIP.SCIP_Real}(0.0)
    SCIP.@SCIP_CALL SCIP.SCIPgetRealParam(scip, "limits/time", time_limit)
    solvingtime = SCIP.SCIPgetSolvingTime(scip)
    if is_infinity(time_limit) == false
        if solvingtime >= time_limit[]
            throw(TimeLimiExceeded())
        end
    end
end
