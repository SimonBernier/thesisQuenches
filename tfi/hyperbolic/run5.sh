#!/bin/bash

set=("$@")
str="Starting run number"
runNumber=0

for i in $(seq 100 1 143);
do
let "runNumber = $(( runNumber + 1 ))"
echo " "
echo "$str $runNumber"

./tfi1d_hyper $i

done
