library(stats4)
library(methods)
# creating the class signal
# 
rm(list = ls())


Peaks <- setRefClass("Peaks",fields = list(peak = "numeric",
                                   variance = "numeric",
                                   intensity = "numeric",
                                   window = "numeric",
                                   lowerlimit = "numeric",
                                   upperlimit = "numeric"),
                     methods = list(
                       print = function()
                       {
                         if(!any(length(peak) == length(variance) | 
                                 length(variance) == length(intensity) |
                                 length(peak == length(intensity)))) 
                           stop("vector peak, variance, and intensity must be same length")
                          
                         # number of peaks
                         np <- length(peak)
                         # create the vector for all the range of retention time
                          x <- rep(0,length(seq(lowerlimit,upperlimit,1)))
                          
                          # for each peak identify the window to be changed
                          sapply(1:np, function(j){
                            if(peak[j]-window[j]>0)
                            {
                              ini<- peak[j]-window[j]
                              fin<- peak[j]+window[j]
                              # obtain position to be changed
                              positionToChange <- seq(ini,fin,1)
                              
                              sapply(1:length(positionToChange),function(pos)
                              {
                                # generate a value from normal dist. with mean equal to
                                # intensity
                                x[pos]<-x[pos]+abs(rnorm(1,intensity[j],sqrt(variance[j])))
                                # exponential decay when the retention time differ from 
                                # the position of the peak
                                x[pos]<-x[pos]/exp(abs(peak[j]-pos)/2)
                              })                  
                              
                            }
                          })
                        print = x 
                       }
                     ))



peakset<-Peaks(peak =c(500,1100,1400,1750,1900,2300,2500),
               variance =c(2,2,1,1,1,1,1),
               intensity = c(10,35,55,35,180,700,400),
               window = c(3,3,3,5,5,10,10),
               lowerlimit = 0,
               upperlimit = 5000)


GC1<- peakset$print()
plot(GC1,type = "l")

x1 <- rep(0,5000)
for(j in 1:7 )
{
  if(peakset$peak[j]-peakset$window[j]>0)
  {
    ini<- peakset$peak[j]-peakset$window[j]
    fin<- peakset$peak[j]+peakset$window[j]
    # obtain position to be changed
    positionToChange <- seq(ini,fin,1)
    cat("\n",positionToChange,"\n")
    f1 <- abs(peakset$peak[j]-positionToChange)
    cat("\n",f1,"\n")
    flag <- 1
    
    x1[ini:fin]<-x1[ini:fin] + sapply(positionToChange,function(pos)
    {
      # generate a value from normal dist. with mean equal to
      # intensity
      aux1 <- abs(rnorm(1,peakset$intensity[j],sqrt(peakset$variance[j]) ) )/exp( abs(pos-peakset$peak[j]))
    })                  
  }
}

x1
summary(x1)
plot(x1,type ="l")

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

  
