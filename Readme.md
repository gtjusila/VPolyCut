# Implementing VPolyhedral Cut with SCIP

## Structure
---

The main code of the module which combines all codes using include is the `VPolyhedralCut.jl`. To implement a cutting plane in SCIP.jl, one should implement a subtype `T` of `SCIP.AbstractSparator` and a corresponding method `SCIP.exec_lp(data::T)`. The subtype is defined in `VPCStructures.jl` and the method is defined in `VPCexeclp.jl`.

The code is divided into four main folders
- `branchandbound` defines a structure `BranchAndBound` which holds the data for a branch and bound tree
- `disjunction` defines a structure `Disjunction` which holds the data for a disjunction. There is an adaptor that takes leaves of a `BranchAndBound` tree into a `Disjunction` object
- `pointray` defines a structure `PointRayCollection`