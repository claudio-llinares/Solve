#!/bin/bash
export OMP_NUM_THREADS=8
echo
echo Bringing source from ../../../../../source_solve/source/
cp -r ../../../../../source_solve/source/* ../source
echo
echo run_0064_0064mpc_newt
cd run_0064_0064mpc_newt; time ../../source//bin/solve >& run.log; cd ../
a=$( tail -n 2 run_0064_0064mpc_newt/run.log | head -n 1 )
echo run_0064_0064mpc_newt $a >> times.log

echo run_0064_0064mpc_mond_small_a0_without_m
cd run_0064_0064mpc_mond_small_a0_without_m; time ../../source//bin/solve >& run.log; cd ../
a=$( tail -n 2 run_0064_0064mpc_mond_small_a0_without_m/run.log | head -n 1 )
echo run_0064_0064mpc_mond_small_a0_without_m $a >> times.log

echo run_0064_0064mpc_mond_without_m
cd run_0064_0064mpc_mond_without_m; time ../../source//bin/solve >& run.log; cd ../
a=$( tail -n 2 run_0064_0064mpc_mond_without_m/run.log | head -n 1 )
echo run_0064_0064mpc_mond_without_m $a >> times.log

echo run_0064_0064mpc_mond_with_m
cd run_0064_0064mpc_mond_with_m; time ../../source//bin/solve >& run.log; cd ../
a=$( tail -n 2 run_0064_0064mpc_mond_with_m/run.log | head -n 1 )
echo run_0064_0064mpc_mond_with_m $a >> times.log

echo run.sh:done
