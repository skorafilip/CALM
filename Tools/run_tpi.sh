#!/bin/bash


#./create_list.sh $1


tpi_path=$HOME"/Downloads/tpi_CALM-master/tpi/"
(find events_EventType${1}/ -type f)>$tpi_path"list.txt"

echo "Tpi runnig"

$tpi_path"tpi" 202 0 0 $2 10 1000000 1 0 1 0 $tpi_path"list.txt"

#cp $tpi_path"outfilecf202a.root" ./

echo "outputfile copied"
