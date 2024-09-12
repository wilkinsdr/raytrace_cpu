#!/bin/bash

run_angdist()
{
	cat > ../par/angdist_jetpoint.par <<EOF
../dat/angdist_jetpoint_h${1}_V${2}.dat
0 $1 1E-3 1.507
$2
0.998
0.001 0.001
-1 1
-3.141 3.141
0.01
1 -3
EOF

	./angdist_jetpoint
}

cd ../bin

#run_angdist 5 0
#run_angdist 5 0.01
#run_angdist 5 0.05
#run_angdist 5 0.1
#run_angdist 5 0.2
run_angdist 5 0.3
#run_angdist 5 0.4
#run_angdist 5 0.5
#run_angdist 5 0.6
