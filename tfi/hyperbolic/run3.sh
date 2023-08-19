#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 50 1 74);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./tfi1d_hyper $i

done
