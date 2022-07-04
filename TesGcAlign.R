source("GC_Sim.R")


T1 = c(0,1,1,2,5,6,5,5,4,4,3,2,1,1,0,0)
X1 = c(0,0,1,1,1,2,3,5,5,6,6,5,3,1,1,0)

X2  = c(0,0,1,1,1,2,3,5,5,6,6,5,2,1,1,0)


slack = 1
Seg = 5
Slack = 1
Options = c(0,1,1,0,0)

# align signal X1
align(T1,X1,Seg,Slack)
alignGCMS(T1,X1,Seg,Slack)

# align signal X2
align(T1,X2,Seg,Slack)
alignGCMS(T1,X2,Seg,Slack)

MX <- cbind(X1,X2)
MT <- cbind(T1,T1)

alignGCMS(MT,MX,Seg,Slack)




# align signal X1 and X2


cow(T,X,Seg,Slack,Options)



