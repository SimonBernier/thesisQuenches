#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 32 32 256);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./tfi1dCrit $i 2.4

done
