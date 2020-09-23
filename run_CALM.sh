#!/bin/bash
#args: EventType, NumberOfEvents
cp events.ini events.ini.bak
NumberOfCores=4

#0.setting amount of events to generate
let amount=$2/$NumberOfCores
sed -i "s,^EventType =.*$,EventType = $1," events.ini
sed -i "s,^NumberOfEvents =.*$,NumberOfEvents = $amount," events.ini

#1.creating root directory for simulation
EventDir=$(awk -F "=" '/EventDir/ {print $2}' events.ini)
#EventDir="events_EventType$1/"
mkdir $EventDir


function CALM(){
#2.making name and directiories for each process
sleep ${1}s
declare calm${1}dir="${EventDir}calm${1}/"
eval "mkdir \${calm${1}dir}"

#3.change directiory and run CALM
eval "sed -i \"s,^EventDir =.*$,EventDir = \${calm${1}dir},\" events.ini"

./calm
}

declare command="CALM 1"

for i in $( seq 2 $NumberOfCores )
do
command="${command} & CALM $i"

done

eval $command

sleep 5s
cp events.ini.bak events.ini
rm events.ini.bak
