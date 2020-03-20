#!/bin/bash

simulator_folder=../ogs_kb1
wdc_folder=$simulator_folder/Libs/Contraflow

mkdir -p $wdc_folder

cp build/src/libcontraflow.a $wdc_folder
cp ./src/*.h $wdc_folder
cp ./src/matrix_stru3/*.h $wdc_folder

