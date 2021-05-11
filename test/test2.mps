* min 24 y1 + 60 y2
* s.t.
* 1/2 y1 +   y2 ≥  6,
*   2 y1 + 2 y2 ≥ 14,
*     y1 + 4 y2 ≥ 13,
*     y1 ≥ 0,
*     y2 ≥ 0
*
NAME         test2
OBJSENSE
  MIN
ROWS
 N  Obj 
 G  cs_0_0_0 
 G  cs_0_0_1 
 G  cs_0_0_2 
COLUMNS
    xs_0_0_0  Obj                             24  cs_0_0_2                         1
    xs_0_0_0  cs_0_0_0                       0.5  cs_0_0_1                         2
    xs_0_0_1  cs_0_0_0                         1  cs_0_0_1                         2
    xs_0_0_1  Obj                             60  cs_0_0_2                         4
RHS
    RHS       cs_0_0_0                         6  cs_0_0_1                        14
    RHS       cs_0_0_2                        13
BOUNDS
 PL Bound     xs_0_0_0
 PL Bound     xs_0_0_1
ENDATA
