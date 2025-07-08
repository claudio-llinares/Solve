#!/bin/bash

# Generates a set of input files from a bla.
# Needs a file model.txt and inputs.txt
#
# Arguments:
#            - prefix for model.txt and inputs.txt
#            If it is present, look for prefix.model.txt, is not, 
#            used model.txt
#
# You should do something like this 
# ls ../../../amiga_runs/run_15/OCBMond.s8=0.4.z?.??? > field1.txt
# ls ../../../amiga_runs/run_15/OCBMond.s8=0.4.z?.??? > field2.txt
# in order to get field1.txt, inputs.txt
#
# Run this file in ahf_runs/run_## (the output will go there).
# Input file need places field1, field2.
#

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
exec 5<$fie_file

while read -u 4 inputs
  do
  read -u 5 fields

  # Identify directory and input
  directory=$(echo $inputs | awk '{print $1}' )
  input_file=$(echo $inputs | awk '{print $2}' )

  # Make the directory and go there
  mkdir $directory
  cd $directory
  
  # Identify the fields
  field01=$(echo $fields | awk '{print $1}' )
  field02=$(echo $fields | awk '{print $2}' )
  field03=$(echo $fields | awk '{print $3}' )
  field04=$(echo $fields | awk '{print $4}' )
  field05=$(echo $fields | awk '{print $5}' )
  field06=$(echo $fields | awk '{print $6}' )
  field07=$(echo $fields | awk '{print $7}' )
  field08=$(echo $fields | awk '{print $8}' )
  field09=$(echo $fields | awk '{print $9}' )
  field10=$(echo $fields | awk '{print $10}' )
  field11=$(echo $fields | awk '{print $11}' )
  
  # Solve the problem with the \/
  # Use this if your field is a file
  ffield01=$(echo $field01 | sed 's/\//\\\//g')
  ffield02=$(echo $field02 | sed 's/\//\\\//g')
  ffield03=$(echo $field03 | sed 's/\//\\\//g')
  ffield04=$(echo $field04 | sed 's/\//\\\//g')
  ffield05=$(echo $field05 | sed 's/\//\\\//g')
  ffield06=$(echo $field06 | sed 's/\//\\\//g')
  ffield07=$(echo $field07 | sed 's/\//\\\//g')
  ffield08=$(echo $field08 | sed 's/\//\\\//g')
  ffield09=$(echo $field09 | sed 's/\//\\\//g')
  ffield10=$(echo $field10 | sed 's/\//\\\//g')
  ffield11=$(echo $field11 | sed 's/\//\\\//g')
  #ffield01=$field01 
  #ffield02=$field02 
  #ffield03=$field03 
  #ffield04=$field04 
  #ffield05=$field05 
  #ffield06=$field06 
  #ffield07=$field07 
  #ffield08=$field08 
  #ffield09=$field09 
  #ffield10=$field10
  #ffield11=$field11
 
  # Make the input file
  sed 's/field01/'$ffield01'/' < ../$mod_file |
  sed 's/field02/'$ffield02'/' | 
  sed 's/field03/'$ffield03'/' | 
  sed 's/field04/'$ffield04'/' | 
  sed 's/field05/'$ffield05'/' | 
  sed 's/field06/'$ffield06'/' | 
  sed 's/field07/'$ffield07'/' | 
  sed 's/field08/'$ffield08'/' | 
  sed 's/field09/'$ffield09'/' | 
  sed 's/field10/'$ffield10'/' | 
  sed 's/field11/'$ffield11'/' > $input_file.input

  # Come back to the original directory
  cd ../

done

