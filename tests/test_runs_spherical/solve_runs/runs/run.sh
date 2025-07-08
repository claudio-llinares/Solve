#!/bin/bash
export OMP_NUM_THREADS=8
echo
echo Bringing source from ../../../../source/
cp -r ../../../../source/* ../source
echo
echo run_032_016mpc
cd run_032_016mpc; time ../../source//bin/solve >& run.log; cd ../
a=$( tail -n 2 run_032_016mpc/run.log | head -n 1 )
echo run_032_016mpc $a >> times.log

echo run.sh:done
