#!/bin/bash

# Generates a set of input files from a bla using one field and 
# a list of output files (without the .input).
#
# You should do something like this 
# ls ../../../amiga_runs/run_15/OCBMond.s8=0.4.z?.??? > field1.txt
# ls ../../../amiga_runs/run_15/OCBMond.s8=0.4.z?.??? > field2.txt
# in order to get field1.txt, inputs.txt
#
# Run this file in ahf_runs/run_## (the output will go there).
# Input file need places field1, field2.
#

executable="time ../../source/lingers"
rm run.sh

# Identify  inputs files of this program
if [ "$1" ]; then
    inp_file=$1.inputs.txt
    mod_file=$1.model.txt
    fie_file=$1.fields.txt
else
    inp_file=inputs.txt
    mod_file=model.txt
    fie_file=fields.txt
fi
exec 4<$inp_file

echo \#!\/bin\/bash >> run.sh
#echo \#!\/bin\/tcsh >> run.sh
#echo setenv OMP_NUM_THREADS 5 >> run.sh
while read -u 4 inputs
  do	
	# Identify directory and input
    directory=$(echo $inputs | awk '{print $1}' )
    input_file=$(echo $inputs | awk '{print $2}' )

	a=$(echo $input_file | awk '{if(substr($0,1,1)=="#") print 0; else print 1;}')
	if [ "$a" = "1" ]; then
	    echo "echo "$directory >> run.sh
		echo cd $directory\; rm lingerg.dat linger.dat\; time $executable \< $input_file.input \>\& $input_file.log\; cd ../ >> run.sh  
	fi
done

echo echo run.sh\:done >> run.sh

chmod u+x run.sh


