library(stats4)
library(methods)
# creating the class signal
# 
rm(list = ls())


Peaks <- setRefClass("Peaks",fields = list(peak = "numeric",
                                   variance = "numeric",
                                   intensity = "numeric",
                                   lowerlimit = "numeric",
                                   upperlimit = "numeric") )


peakset<-Peaks(peak =c(500,1100,1400,1750,1900,2300,2500),
               variance =c(2,2,1,1,1,1,1),
               intensity = c(1,0.25,0.25,0.25,180,700,400),
               lowerlimit = 0,
               upperlimit = 5000)




setClass("Peaks",
         slots = list(peak = "numeric", 
                      variance = "numeric", 
                      intensity = "numeric",
                      lowerlimit = "numeric",
                      upperlimit = "numeric"))





peakset<-Peaks(peak =c(500,1100,1400,1750,1900,2300,2500),
               variance =c(2,2,1,1,1,1,1),
               intensity = c(1,0.25,0.25,0.25,180,700,400),
               lowerlimit = 0,
               upperlimit = 5000)

  
