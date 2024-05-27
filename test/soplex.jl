using Test
import Polyhedra
import SCIP
import MathOptInterface as MOI

# Helper Functions 

function print_lpi_information(lpi)

    println("====================")
    println("Printing LP Information") 

    col_num = Ref{Cint}(0)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetNCols(lpi[],col_num)
    row_num = Ref{Cint}(0)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetNRows(lpi[],row_num)
    col_num = Int(col_num[])
    row_num = Int(row_num[])
    println()
    println("LP have $col_num columns and $row_num rows")

    primal = zeros(SCIP.SCIP_Real, col_num)
    dual =  zeros(SCIP.SCIP_Real, row_num)
    activity =  zeros(SCIP.SCIP_Real, row_num)
    red_cost = zeros(SCIP.SCIP_Real, col_num)
    SCIP.@SCIP_CALL SCIP.SCIPlpiGetSol(lpi[],C_NULL, primal, dual, activity, red_cost)
    
    println()
    println("Primal Solution")
    println("----------")
    println(primal)

    println()
    println("Dual Solution")
    println("----------")
    println(dual)
    
    println()
    println("Row Activity")
    println("----------")
    println(activity)
    
    println()
    println("Reduced Cost")
    println("----------")
    println(red_cost)
    
    cstat = Cint(-1)*ones(Cint,col_num)
    rstat = Cint(-1)*ones(Cint,row_num)

    SCIP.@SCIP_CALL SCIP.SCIPlpiGetBase(lpi[],cstat,rstat)
    
    println()
    println("Column Basis Status")
    println("----------")

    for (i,col) in enumerate(cstat)
        @assert(col >= 0)
        col = SCIP.SCIP_BaseStat(Int(col))
        println("Column "*string(i)*" status "*string(col))
    end

    println()
    println("Row Basis Status")
    println("----------")
    for (i,row) in enumerate(rstat)
        row = SCIP.SCIP_BaseStat(Int(row))
        println("Row "*string(i)*" status "*string(row))
    end

    println("====================")
end 

function run_test_soplex() 
    # Create a new model
    optimizer = SCIP.Optimizer()
    inner = optimizer.inner
    scip = inner.scip[]
    
    lpi = Ref{Ptr{SCIP.SCIP_LPI}}(C_NULL)
    SCIP.@SCIP_CALL SCIP.SCIPlpiCreate(lpi, C_NULL, "Test", SCIP.SCIP_OBJSEN_MINIMIZE)
    SCIP.@SCIP_CALL SCIP.SCIPlpiReadLP(lpi[],joinpath(@__DIR__,"data","soplex_test.lp"))
    SCIP.@SCIP_CALL SCIP.SCIPlpiSolvePrimal(lpi[])

    print_lpi_information(lpi)
end

run_test_soplex()