# Implementing VPolyhedral Cut with SCIP

## Structure
---

The main code of the module which combines all codes using include is the `VPolyhedralCut.jl`. To implement a cutting plane in SCIP.jl, one should implement a subtype `T` of `SCIP.AbstractSparator` and a corresponding method `SCIP.exec_lp(data::T)`. The subtype is defined in `VPCStructures.jl` and the method is defined in `VPCexeclp.jl`.

The code is divided into six main folders
- `nonbasicspace` defines a structure NonbasicSpace which assist in converting points and rays from the structural problem of the variable into the nonbasic space
- `branchandbound` defines a structure `BranchAndBound` which holds the data for a branch and bound tree
- `disjunction` defines a structure `Disjunction` which holds the data for a disjunction. There is an adaptor that takes leaves of a `BranchAndBound` tree into a `Disjunction` object
- `pointray` defines the structure `PointRayCollection` an object which can be created from a disjunction. It will collect points and rays from a disjunction.
- `prlp` defines the structure `PRLP` which allows us to construct the separating LP and manage its solution.
- `cutpool` defines the structure `CutPool` which manages the addition of cuts into SCIP.
- `commons` contains various utility function definitions such as numerics, log helpers and helper functions for working with SCIP tableau.