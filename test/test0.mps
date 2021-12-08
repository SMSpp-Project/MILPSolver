* this is known in some circles as "the PIntel problem"
* min  - 5 * x_1 - 2 * x_2
* s.t.   2 * x_1 +     x_2 <= 9
*            x_1           <= 4
*                      x_2 <= 7
*            x_1           >= 0
*                      x_2 >= 0
*
NAME          Pintel
OBJSENSE
  MIN
ROWS
 N  Obj 
 L  cs_0_0_0
COLUMNS
    xs_0_0_0  Obj                             -5  cs_0_0_0                         2
    xs_0_0_1  Obj                             -2  cs_0_0_0                         1
RHS
    RHS       cs_0_0_0                         9
BOUNDS
 UP Bound     xs_0_0_0                         4
 UP Bound     xs_0_0_1                         7
ENDATA
