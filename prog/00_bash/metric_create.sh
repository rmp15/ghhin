#!/bin/bash

# this script
# downloads the desired climate variable globally for the desired model, climate scenario in the chosen years
# picks out a location for creating a time series
# processes the netcdf file year-by-year
# binds the time series together

clear

declare -a years=($(seq 1950 1950))
declare -a rcps=('rcp45' 'rcp85')
declare -a vars=('tasmin' 'tasmax')
declare -a models=('ACCESS1-0')

for year in "${years[@]}"; do
for rcp in "${rcps[@]}"; do
for var in "${vars[@]}"; do
for model in "${models[@]}"; do

# create file architecture to store files
dir="//data/climate/net_cdf/" ; dir+="$var/" ; dir+="$rcp/" ; dir+="$model/" ; dir+="raw/" ;
mkdir -p "$dir"

# downloads the desired climate variable globally
wget http://nasanex.s3.amazonaws.com/NEX-GDDP/BCSD/rcp85/day/atmos/tasmin/r1i1p1/v1.0/tasmin_day_BCSD_rcp85_r1i1p1_ACCESS1-0_2059.nc -P ~/Desktop

done; done; done; done;

declare lon=
declare lat=

# process a particular location
