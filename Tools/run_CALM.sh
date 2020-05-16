#!/bin/bash

#0.setting amount of events to generate
let amount=$2/4
sed -i "32 s/.*/NumberOfEvents = $amount/" events.ini
sed -i "83 s/.*/EventType = $1/" events.ini

#1.making name and directiories for each process

#EventDir="eventsReggae_minijets_10_20_newTotalMom"
dirname="events_EventType$1"
EventDir=$dirname
mkdir $EventDir
EventDir=$EventDir"/"
#77min 79max 83eventtype


for i in {0..3}
do
declare calm${i}dir="${EventDir}calm${i}dir"
eval "mkdir \${calm${i}dir}"

done


#2.change directiory and run CALM
n="67"
#line="EventDir = eventsReggae_minijets_10_20_newTotalMom"
line="EventDir = $dirname"

for i in {0..3}
do

sed -i "$n s/.*/$line\/calm${i}dir\//" events.ini
gnome-terminal -- ./calm
sleep 1s

done

sleep 5s
sed -i "$n s/.*/EventDir = eventsReggae_minijets_10_20_newTotalMom//" events.ini


