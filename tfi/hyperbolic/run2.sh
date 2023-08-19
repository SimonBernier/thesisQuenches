#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 25 1 49);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./tfi1d_hyper $i

done
