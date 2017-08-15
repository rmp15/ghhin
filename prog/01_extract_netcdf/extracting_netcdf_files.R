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

# function to load appropriate climate data
climate.timeseries = function(dname,scenario,model,rest,date) {

    ncname <- paste0(dname,'_',scenario,'_',model,'_',rest,'.nc')
    
    # open NetCDF file
    ncin <- nc_open(paste0('~/data/wmo/net_cdf/',dname,'/raw/',ncname))
    
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
    t.names <- as.Date(t,origin = date)
    
    # stamp as character names
    timeStamp <-  strptime(t.names,"%Y-%m-%d")
    
    # dataframe with time and temp
    dat0 <- data.frame(time=timeStamp,temp=temp)
    dat0 <- na.omit(dat0)
    
    # test for days in Ahmedabad schema
    dat0$color <- ifelse(dat0$temp<41.1,'white',
    ifelse(dat0$temp<43.1,'yellow',
    ifelse(dat0$temp<45,'orange','red')))
    
    # add year
    dat0$year =  as.numeric(substr(dat0$time,1,4))
    
    return(dat0)
}

# 1961-2000
dat0 = climate.timeseries('itmax','A2','csiro_mk3','0_20c3m_72.5714E_23.0225N_n_su_01','1961-01-01')

# 2046-2065
dat = climate.timeseries('itmax','A2','csiro_mk3','0_sresa1b_20_72.5714E_23.0225N_n_su','2046-01-01')

# 2081-2100
dat2 = climate.timeseries('itmax','A2','csiro_mk3','0_sresa1b_21_72.5714E_23.0225N_n_su','2081-01-01')

# create statistics for two time periods
dat.stats2 <- ddply(dat2,.(month),summarize,days=sum(count))

# limits for plots
ymax = max(dat0$temp,dat$temp,dat2$temp) + 2
ymin = min(dat0$temp,dat$temp,dat2$temp) - 2

# plot timeseries

pdf(paste0('../../output/plot1.pdf'),height=0,width=0,paper='a4r')

print(ggplot() +
geom_line(data=subset(dat0,year %in% c(1981:2000)),aes(x=time,y=temp)) +
geom_point(data=subset(dat0,year %in% c(1981:2000)),aes(x=time,y=temp,color=color)) +
ylim(c(ymin,ymax)) +
ggtitle('A2 Ahmedabad 1981-2000') +
scale_color_manual(values=c("white"="black", "yellow"="yellow", "orange"="orange", "red"="red")) +
guides(color=FALSE) +
theme_bw())

dev.off()

pdf(paste0('../../output/plot2.pdf'),height=0,width=0,paper='a4r')

print(ggplot() +
geom_line(data=dat,aes(x=time,y=temp)) +
geom_point(data=dat,aes(x=time,y=temp,color=color)) +
ylim(c(ymin,ymax)) +
ggtitle('A2 Ahmedabad 2046-2065') +
scale_color_manual(values=c("white"="black", "yellow"="yellow", "orange"="orange", "red"="red")) +
guides(color=FALSE) +
theme_bw())

dev.off()

pdf(paste0('../../output/plot3.pdf'),height=0,width=0,paper='a4r')

print(ggplot() +
geom_line(data=dat2,aes(x=time,y=temp)) +
geom_point(data=dat2,aes(x=time,y=temp,color=color)) +
ylim(c(ymin,ymax)) +
ggtitle('A2 Ahmedabad 2081-2100') +
scale_color_manual(values=c("white"="black", "yellow"="yellow", "orange"="orange", "red"="red")) +
guides(color=FALSE) +
theme_bw())

dev.off()


# save timeseries as RDS
#saveRDS()



















