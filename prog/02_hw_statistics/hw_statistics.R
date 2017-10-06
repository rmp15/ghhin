# Extracting NetCDF files

# From:
# https://github.com/rmp15/climate/blob/master/extract_netcdf_data.R
# http://geog.uoregon.edu/bartlein/courses/geog607/Rmd/netCDF_01.htm

rm(list=ls())

library(RColorBrewer)
library(lattice)
library(ncdf4)
library(lubridate)
library(ggplot2)
library(plyr)

# arguments from Rscript
args <- commandArgs(trailingOnly=TRUE)

# year of interest
#year <- as.numeric(args[1])

#print(paste0('running extracting_netcdf_files.R for ',year))

# load for 1961-2000

# names for files
dname <- 'itmax'
#dname <- as.character(args[2])
scenario <- 'A2'
#freq <- as.character(args[3])
model <- 'csiro_mk3'
#num <- as.character(args[4])
rest <- '0_20c3m_72.5714E_23.0225N_n_su_00'
ncname <- paste0(dname,'_',scenario,'_',model,'_',rest,'.nc')

#year <- unlist(strsplit(ncname, "_"))
#year <- year[5]
#year <- unlist(strsplit(year, ".nc"))

# open NetCDF file
ncin <- nc_open(paste0('~/git/wmo/data/knmi/ahmedabad/ghhin/1961_2000/',dname,'/raw/',ncname))

# extract climate variable
varname <- 'tasmax'
temp <- ncvar_get(ncin, varname)
#nlon <- dim(temp)

#lat <- ncvar_get(ncin, "latitude", verbose = F)
#nlat <- dim(lat)

# get time variable and convert to days
t <- ncvar_get(ncin, "time")
tunits <- ncatt_get(ncin, "time", "units")
nt <- dim(t)

# global attributes
title <- ncatt_get(ncin, 0, "title")
institution <- ncatt_get(ncin, 0, "institution")
datasource <- ncatt_get(ncin, 0, "source")
references <- ncatt_get(ncin, 0, "references")
history <- ncatt_get(ncin, 0, "history")
Conventions <- ncatt_get(ncin, 0, "Conventions")

# close NetCDF file
nc_close(ncin)

# split the time units string into fields
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
#t.names <- as.POSIXct(t, origin = "2046-01-01")
t.names <- as.Date(t,origin = "1961-01-01")

# stamp as character names
timeStamp <-  strptime(t.names,"%Y-%m-%d")

# dataframe with time and temp
dat0 <- data.frame(time=timeStamp,temp=temp)
dat0 <- na.omit(dat0)

# test for days in Ahmedabad schema
dat0$color <- ifelse(dat0$temp<41,'white',
ifelse(dat0$temp<43.4,'yellow',
ifelse(dat0$temp<45,'orange','red')))

# load for 2046-2065

# names for files
dname <- 'itmax'
#dname <- as.character(args[2])
scenario <- 'A2'
#freq <- as.character(args[3])
model <- 'csiro_mk3'
#num <- as.character(args[4])
rest <- '0_sresa1b_20_72.5714E_23.0225N_n_su'
ncname <- paste0(dname,'_',scenario,'_',model,'_',rest,'.nc')

#year <- unlist(strsplit(ncname, "_"))
#year <- year[5]
#year <- unlist(strsplit(year, ".nc"))

# open NetCDF file
ncin <- nc_open(paste0('~/data/ghhin/net_cdf/',dname,'/raw/',ncname))

# extract climate variable
varname <- 'tasmax'
temp <- ncvar_get(ncin, varname)
#nlon <- dim(temp)

#lat <- ncvar_get(ncin, "latitude", verbose = F)
#nlat <- dim(lat)

# get time variable and convert to days
t <- ncvar_get(ncin, "time")
tunits <- ncatt_get(ncin, "time", "units")
nt <- dim(t)

# global attributes
title <- ncatt_get(ncin, 0, "title")
institution <- ncatt_get(ncin, 0, "institution")
datasource <- ncatt_get(ncin, 0, "source")
references <- ncatt_get(ncin, 0, "references")
history <- ncatt_get(ncin, 0, "history")
Conventions <- ncatt_get(ncin, 0, "Conventions")

# close NetCDF file
nc_close(ncin)

# split the time units string into fields
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
#t.names <- as.POSIXct(t, origin = "2046-01-01")
t.names <- as.Date(t,origin = "2046-01-01")

# stamp as character names
timeStamp <-  strptime(t.names,"%Y-%m-%d")

# dataframe with time and temp
dat <- data.frame(time=timeStamp,temp=temp)
dat <- na.omit(dat)

# test for days in Ahmedabad schema
dat$color <- ifelse(dat$temp<41,'white',
                ifelse(dat$temp<43.4,'yellow',
                    ifelse(dat$temp<45,'orange','red')))

# load for 2081-2100

# names for files
dname <- 'itmax'
#dname <- as.character(args[2])
scenario <- 'A2'
#freq <- as.character(args[3])
model <- 'csiro_mk3'
#num <- as.character(args[4])
rest <- '0_sresa1b_21_72.5714E_23.0225N_n_su'
ncname <- paste0(dname,'_',scenario,'_',model,'_',rest,'.nc')

#year <- unlist(strsplit(ncname, "_"))
#year <- year[5]
#year <- unlist(strsplit(year, ".nc"))

# open NetCDF file
ncin <- nc_open(paste0('~/data/ghhin/net_cdf/',dname,'/raw/',ncname))

# extract climate variable
varname <- 'tasmax'
temp <- ncvar_get(ncin, varname)
#nlon <- dim(temp)

#lat <- ncvar_get(ncin, "latitude", verbose = F)
#nlat <- dim(lat)

# get time variable and convert to days
t <- ncvar_get(ncin, "time")
tunits <- ncatt_get(ncin, "time", "units")
nt <- dim(t)

# global attributes
title <- ncatt_get(ncin, 0, "title")
institution <- ncatt_get(ncin, 0, "institution")
datasource <- ncatt_get(ncin, 0, "source")
references <- ncatt_get(ncin, 0, "references")
history <- ncatt_get(ncin, 0, "history")
Conventions <- ncatt_get(ncin, 0, "Conventions")

# close NetCDF file
nc_close(ncin)

# split the time units string into fields
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth = as.integer(unlist(tdstr)[2])
tday = as.integer(unlist(tdstr)[3])
tyear = as.integer(unlist(tdstr)[1])
#t.names <- as.POSIXct(t, origin = "2046-01-01")
t.names <- as.Date(t,origin = "2081-01-01")

# stamp as character names
timeStamp <-  strptime(t.names,"%Y-%m-%d")

# dataframe with time and temp
dat2 <- data.frame(time=timeStamp,month=month(timeStamp),temp=temp,count=1)
dat2 <- na.omit(dat2)

# test for days in Ahmedabad schema
dat2$color <- ifelse(dat2$temp<41,'white',
ifelse(dat2$temp<43.4,'yellow',
ifelse(dat2$temp<45,'orange','red')))

# create statistics for two time periods
dat.stats2 <- ddply(dat2,.(month),summarize,days=sum(count))

# limits for plots
ymax = max(dat0$temp,dat$temp,dat2$temp) + 2
ymin = min(dat0$temp,dat$temp,dat2$temp) - 2

# plot timeseries

ggplot() +
geom_line(data=subset(dat0,year %in% c(1981,2000)),aes(x=time,y=temp)) +
geom_point(data=dat0,aes(x=time,y=temp,color=color)) +
ylim(c(ymin,ymax)) +
ggtitle('A2 Ahmedabad 1981-2000') +
scale_color_manual(values=c("white"="black", "yellow"="yellow", "orange"="orange", "red"="red")) +
theme_bw()

ggplot() +
geom_line(data=dat,aes(x=time,y=temp)) +
geom_point(data=dat,aes(x=time,y=temp,color=color)) +
ylim(c(ymin,ymax)) +
ggtitle('A2 Ahmedabad 2046-2065') +
scale_color_manual(values=c("white"="black", "yellow"="yellow", "orange"="orange", "red"="red")) +
theme_bw()

ggplot() +
geom_line(data=dat2,aes(x=time,y=temp)) +
geom_point(data=dat2,aes(x=time,y=temp,color=color)) +
ylim(c(ymin,ymax)) +
ggtitle('A2 Ahmedabad 2081-2100') +
scale_color_manual(values=c("white"="black", "yellow"="yellow", "orange"="orange", "red"="red")) +
theme_bw()

# save timeseries as RDS
saveRDS(



















