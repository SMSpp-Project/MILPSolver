* min x1 +3x2 +3x3 +x4
* s.t.
* 3x1 +4x2 −3x3 +x4 = 2,
* 3x1 −2x2 +6x3 −x4 = 1,
* 6x1 +4x2      +x4 = 4
*
* x1, x2, x3, x4 ≥ 0
*
NAME          test1
OBJSENSE
  MIN
ROWS
 N  Obj
 E  cs_0_0_0
 E  cs_0_0_1
 E  cs_0_0_2
COLUMNS
    xs_0_0_0  Obj                              1  cs_0_0_2                         6
    xs_0_0_0  cs_0_0_1                         3  cs_0_0_0                         3
    xs_0_0_1  Obj                              3  cs_0_0_0                         4
    xs_0_0_1  cs_0_0_2                         4  cs_0_0_1                        -2
    xs_0_0_2  cs_0_0_0                        -3  Obj                              3
    xs_0_0_2  cs_0_0_1                         6
    xs_0_0_3  cs_0_0_1                        -1  cs_0_0_0                         1
    xs_0_0_3  Obj                              1  cs_0_0_2                         1
RHS
    RHS       cs_0_0_0                         2  cs_0_0_1                         1
    RHS       cs_0_0_2                         4
BOUNDS
 PL Bound     xs_0_0_0                           
 PL Bound     xs_0_0_1                           
 PL Bound     xs_0_0_2                           
 PL Bound     xs_0_0_3                           
ENDATA
