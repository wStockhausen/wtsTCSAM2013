#!/bin/sh
echo on
DIR="$( cd "$( dirname "$0" )" && pwd )"
cd ${DIR}
cp ../tcsam2013 ./tcsam2013
./tcsam2013 -rs -nox  -configFile ../data/TCSAM2013e_ModelConfig.txt -nohess
echo off
