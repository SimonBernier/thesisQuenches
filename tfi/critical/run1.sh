#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 16 32 240);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./tfi1dCrit $i

done
