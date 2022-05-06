source("GC_Sim.R")

target_signal =  c(0,1,1,2,5,6,5,5,4,4,3,2,1,1,0,0)
x = 1:length(target_signal)

sample_signal = c(0,0,1,1,1,2,3,5,5,6,6,5,3,1,1,0)
xsignal = 1:length(sample_signal)

P<-list(x = xsignal, f = sample_signal)
P

slack = 1
segments = 3

