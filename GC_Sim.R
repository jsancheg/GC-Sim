library(stats4)
library(methods)
# creating the class signal
# 
rm(list = ls())

movies <- setRefClass("movies", fields = list(name = "character",
                                              leadActor = "character", rating = "numeric"))

#now we can use the generator to create objects
movieList <- movies(name = "Iron Man", 
                    leadActor = "Robert downey Jr", rating = 7)
movieList


Peaks <- setRefClass("Peaks",fields = list(peak = "numeric",
                                   variance = "numeric",
                                   intensity = "numeric",
                                   lowerlimit = "numeric",
                                   upperlimit = "numeric") )

peak1 <- Peaks(peak = 500,
               variance = 2,
               intensity = 20)

peakset<-Peaks(peak =c(500,1100,1400,1750,1900,2300,2500),
               variance =c(2,2,1,1,1,1,1),
               intensity = c(1,0.25,0.25,0.25,180,700,400))




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

  
