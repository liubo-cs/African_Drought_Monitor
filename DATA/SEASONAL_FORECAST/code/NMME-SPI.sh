#!/bin/bash
#source /opt/intel/composer_xe_2013/bin/compilervars.sh intel64
#ifort -g -traceback -fpe0 -CB NMME-SPI.f90 -L/home/xingy/netcdf-4.1.1/lib -lnetcdf -I/home/xingy/netcdf-4.1.1/include -o NMME-SPI.exe
ulimit -s unlimited
#/opt/intel/composer_xe_2013/bin/ifort -g -traceback -fpe0 -CB NMME-SPI.f90 -lnetcdf -lnetcdff -o NMME-SPI.exe
#/opt/intel/composer_xe_2013/bin/ifort -g -traceback -fpe0 -CB NMME-SPI.f90 -L/home/freeze/water_monitor/African_Drought_Monitor/LIBRARIES/netcdf-4.1.1/netcdf-4.1.1/lib -lnetcdf -I/home/freeze/water_monitor/African_Drought_Monitor/LIBRARIES/netcdf-4.1.1/netcdf-4.1.1/include -o NMME-SPI.exe
#ifort -g -traceback -fpe0 -CB NMME-SPI.f90 -L/home/freeze/water_monitor/African_Drought_Monitor/LIBRARIES/NETCDF-4.2/lib -lnetcdff -lnetcdf -I/home/freeze/water_monitor/African_Drought_Monitor/LIBRARIES/NETCDF-4.2/include -o NMME-SPI.exe
#gfortran NMME-SPI_NEW.f90 -I/home/stream1/monitor/African_Drought_Monitor/LIBRARIES/NETCDF-4.2/include -L/home/stream1/monitor/African_Drought_Monitor/LIBRARIES/NETCDF-4.2/lib -lnetcdf -lhdf5_hl -lhdf5 -lz -lm -lcurl -o NMME-SPI.exe
gfortran NMME-SPI_NEW.f90 -I/usr/lib64/gfortran/modules -lnetcdff -o NMME-SPI.exe
#ifort -g -traceback -fpe0 -CB NMME-SPI.f90 -L/usr/local/lib -lnetcdf -lhdf5_hl -lhdf5 -lz -lm -lcurl -I/usr/local/include -o NMME-SPI.exe
#export LD_PRELOAD=/home/stream1/monitor/African_Drought_Monitor/LIBRARIES/netcdf-4.1.1_intel_compiled/lib/libnetcdf.so.7
#ifort -O3 -xHost -ip -no-prec-div -g -traceback -fpe0 -CB -static-intel NMME-SPI.f90 -L/home/stream1/monitor/African_Drought_Monitor/LIBRARIES/netcdf-4.1.1_intel_compiled/lib -lnetcdff -lnetcdf -I/home/stream1/monitor/African_Drought_Monitor/LIBRARIES/netcdf-4.1.1_intel_compiled/include -o NMME-SPI.exe
./NMME-SPI.exe
rm NMME-SPI.exe
