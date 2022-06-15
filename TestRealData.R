
ruta<-getwd()
Data <- readRDS(paste0(ruta,"/Data/sample_data.rds"))
str(Data)
Data[1,,]
dim(Data[1,1:2900,])
plot(Data[1,,100],type = "l")
plot(Data[2,,100],type = "l")


T<-Data[1,,100]
X<-Data[2,,100]

Seg = 5
Slack = 1
Options = c(0,1,1,0,0)

cow(T,X,Seg,Slack,Options)

length(T)
length(X)

align(T,X,10,15)
