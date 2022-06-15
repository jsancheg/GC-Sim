source("GC_Sim.R")


T = c(0,1,1,2,5,6,5,5,4,4,3,2,1,1,0,0)
X = c(0,0,1,1,1,2,3,5,5,6,6,5,3,1,1,0)
  
slack = 1
Seg = 5
Slack = 1
Options = c(0,1,1,0,0)

cow(T,X,Seg,Slack,Options)



