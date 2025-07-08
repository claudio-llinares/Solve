#!/bin/bash



#executable="numactl --cpunodebind=0 --membind=0 time ../../source/bin/solve"
#executable="nice +19 time ../../source/bin/solve"
executable="time ../../source//bin/solve"
#executable="time mpirun -n 8 ../../source/godunov"
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
#echo setenv OMP_NUM_THREADS 56 >> run.sh
echo export OMP_NUM_THREADS=8 >> run.sh

echo "echo" >> run.sh
echo "echo Bringing source from ../../../../../source_solve/source/" >> run.sh
echo cp -r '../../../../../source_solve/source/* ../source' >> run.sh
echo "echo" >> run.sh

while read -u 4 inputs
  do	
	# Identify directory and input
    directory=$(echo $inputs | awk '{print $1}' )
    input_file=$(echo $inputs | awk '{print $2}' )

	a=$(echo $input_file | awk '{if(substr($0,1,1)=="#") print 0; else print 1;}')
	if [ "$a" = "1" ]; then
            #for i in `seq 1 5`; do
	        echo "echo "$directory >> run.sh
                echo cd $directory\; $executable \>\& $input_file.log\; cd ../ >> run.sh  #\&
                echo a=\$\( tail -n 2 $directory/run.log \| head -n 1 \) >> run.sh
                echo "echo "$directory \$a \>\> times.log >> run.sh
                echo "" >> run.sh
            #done
	fi
done

echo echo run.sh\:done >> run.sh

chmod u+x run.sh


