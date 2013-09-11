#!/bin/sh
#PBS -l nodes=1:ppn=1,walltime=3:59:00
source /home/xingy/.bashrc
cd /tigress/xingy/nate/code
./NMME-SPI.exe
rm -f NMME-SPI.exe
