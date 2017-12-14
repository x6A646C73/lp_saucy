NAME    simple
ROWS
 N COST
 G LIM1
 G LIM2
 G LIM3
 G LIM4
COLUMNS
    XONE  COST    1   LIM1    1
    XONE  LIM4    1
    XTWO  COST    1   LIM1    1
    XTWO  LIM2    1
    XTHREE  COST    1   LIM2    2
    XTHREE  LIM3    2
    XFOUR  COST    1   LIM3    1
    XFOUR  LIM4    1
RHS
    RHS1    LIM1    1   LIM2    1
    RHS1    LIM3    1   LIM4    1
BOUNDS
 LO  BND1    XONE  0
 UP  BND1    XONE  1
 LO  BND1    XTWO  0
 UP  BND1    XTWO  1
 LO  BND1    XTHREE  0
 UP  BND1    XTHREE  1
 LO  BND1    XFOUR  0
 UP  BND1    XFOUR  1
ENDATA

xone + xtwo + xthree + xfour
lim1: xone + xtwo >= 1
lim2: xtwo + 2*xthree >= 1
lim3: 2*xthree + xfour >= 1
lim4: xfour + xone >= 1

//TODO: build the graph correctly, doing it wrong right now...
1 1 0 0
0 1 2 0
0 0 2 1
1 0 0 1
above implies the following edges and weights:
x1 c1 1
x2 c1 1
x2 c2 1
x3 c2 2
x3 c3 2
x4 c3 1
x1 c4 1
x4 c4 1
