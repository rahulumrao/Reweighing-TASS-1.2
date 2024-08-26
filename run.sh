#!/bin/bash

file_name='input.tass'

cat > $file_name << EOF
NUMBER OF CV
5
CODE NAME
CPMD
NUMBER OF UMBRELLA
23
UCV COLUMN
1
MTD ENABLED
n
mtd cv column
2
mtd bias factor
1200
SYSTEM TEMP
300
CV TEMP
1000
TMIN
1000
TMAX
14000
REWEIGHTING TOOL
pmf
B-SPLINE INTERPOLATION
PROBABILITY DIMENSION
2
PROBABILITY CV INDEX
1 3 
statistical errors block size
4
CV PRINT FREQUENCY
1
mtd print frequency
10
GRIDS
1.8 6.0  0.01
0.5 5.0  0.05
1.0 8.0  0.01
1.0 6.0  0.01
1.0 8.0  0.01
EOF

tass_analysis.x < $file_name
