#!/bin/bash


grep ">" $1 |grep -v "_"|grep -v "M"|sed 's|[>,]||g'|xargs samtools faidx $1
