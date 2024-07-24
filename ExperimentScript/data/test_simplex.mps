* SCIP STATISTICS
*   Problem name     : 
*   Variables        : 2 (0 binary, 2 integer, 0 implicit integer, 0 continuous)
*   Constraints      : 3
NAME          
OBJSENSE
  MAX
ROWS
 N  Obj 
 L  c0 
 L  c1 
 L  c2 
COLUMNS
    INTSTART  'MARKER'                            'INTORG'                           
    x0        Obj                              1  c0                               1 
    x0        c1                              -1 
    x1        Obj                              1  c0                               1 
    x1        c2                            -0.4 
    INTEND    'MARKER'                            'INTEND'                           
RHS
    RHS       c0                             2.5  c1                               0 
    RHS       c2                               0 
BOUNDS
 FR Bound     x0                                 
 FR Bound     x1                                 
ENDATA