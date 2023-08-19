#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 204 1 254);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./1dtfiST $i

done
