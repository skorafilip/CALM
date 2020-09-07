#!/bin/bash
cp events.ini events.ini.bak
cores=4

#0.setting amount of events to generate
let amount=$2/$cores
sed -i "s,^EventType =.*$,EventType = $1," events.ini
sed -i "s,^NumberOfEvents =.*$,NumberOfEvents = $amount," events.ini

#1.making name and directiories for each process
EventDir=$(awk -F "=" '/EventDir/ {print $2}' events.ini)
mkdir $EventDir


for i in $( seq 1 $cores )
do
declare calm${i}dir="${EventDir}calm${i}"
eval "mkdir \${calm${i}dir}"

done


#2.change directiory and run CALM
for i in $( seq 1 $cores )
do

eval "sed -i \"s,^EventDir =.*$,EventDir = \${calm${i}dir},\" events.ini"

gnome-terminal -- ./calm
sleep 1s

done

sleep 5s
cp events.ini.bak events.ini
rm events.ini.bak