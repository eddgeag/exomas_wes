#!/bin/bash

header=$(head $1 -n 1);
id=$(echo $header | head -n 1 | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g');
sm=$(echo $header | head -n 1 | grep -Eo "[ATGCN]+$");
echo "@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA" ;



bwa mem \
-M \
-t 8 \
-v 3 \
-R $(echo "@RG\tID:$id\tSM:$id"_"$sm\tLB:$id"_"$sm\tPL:ILLUMINA") \
$3 \
$1 $2 -o $4 
