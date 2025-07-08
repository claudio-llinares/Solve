sed 's/ run/ p_table/g' inputs.txt > p_table.inputs.txt
#cut -d" " -f 1 inputs.txt > p_table.fields.txt

./gen_inputsdir9.sh
./gen_inputsdir9.sh p_table
#./gen_inputsdir9.sh p
./gen_runsdir.sh
