#!/bin/bash
echo run_lcdm_small
cd run_lcdm_small; rm lingerg.dat linger.dat; time time ../../source/lingers < run.input >& run.log; cd ../
echo run_lcdm_large
cd run_lcdm_large; rm lingerg.dat linger.dat; time time ../../source/lingers < run.input >& run.log; cd ../
echo run_lcdm_bar
#cd run_lcdm_bar; rm lingerg.dat linger.dat; time time ../../source/lingers < run.input >& run.log; cd ../
echo run.sh:done
