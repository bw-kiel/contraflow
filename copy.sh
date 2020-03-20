#!/bin/bash

simulator=../ogs_kb1
wdc_folder=$simulator/Libs/Contraflow

mkdir -p $wdc_folder

cp build/src/libcontraflow.a $wdc_folder
cp ./src/*.h $wdc_folder
cp ./src/matrix_stru3/*.h $wdc_folder

