#!/bin/bash


tpi_path=$HOME"/Downloads/tpi_CALM-master/tpi/"
(find events_EventType${1}/ -type f)>$tpi_path"list.txt"

echo "Tpi runnig"

if [ "$2" -eq "-1" ]
then    
    for i in 0 1 12 16 21
    do

    gnome-terminal -- ./run_tpi.sh $1 $i
    sleep 4s

    done
else
    $tpi_path"tpi" 202 0 0 $2 10 1000000 1 0 1 0 $tpi_path"list.txt"
fi
#cp $tpi_path"outfilecf202a.root" ./

echo "outputfile copied"
