#!/bin/bash
export OMP_NUM_THREADS=12
echo run_0064_00064mpc_01
cd run_0064_00064mpc_01; ../../source/bin/solve run.input >& run.log; cd ../

echo run_0512_00512mpc_01
#cd run_0512_00512mpc_01; ../../source/bin/solve run.input >& run.log; cd ../

echo run.sh:done
