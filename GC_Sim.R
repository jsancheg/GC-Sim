library(stats4)
library(methods)
library(pracma)
library(tictoc)
# creating the class signal.
# 
rm(list = ls())
# Create a class with three fields.
# peak: A vector that contains the position of peaks in the range of spectra.
# variance: A vector that contains the variances for different peaks.
# intensity: A vector that contains the intensity of the peaks.
# Window: A vector that contains the number of retention time to be affected
#         around the main peak.
# lowerlimit: A variable that set the minimum value for retention time.
# upperlimit: A variable that set the maximum value for retention time.


InterpCoeff <- function(n, nprime, offs)
{
  # Function to calculate coefficients for interpolation
   p = length(nprime)
   q = n-1
   Coeff = matrix(0, nrow = p, ncol = n)
   Index = matrix(0, nrow = p, ncol = n)
   output <- list()
   for(i_p in 1:p)
   {
     pp <- 1:nprime[i_p]
     p1 <- (0:q) * (nprime[i_p]-1)/q + 1
     aux <- histc(p1,pp)
     ignore <- aux$cnt
     k <- aux$bin
     if(any(p1<1))
         k[p1<1] = 1
     
     if(any(p1>=nprime[i_p]))
       k[p1>=nprime[i_p]] = nprime[i_p]-1
     
     Coeff[i_p,] = (p1-pp[k])
     Index[i_p,] = k - offs[i_p]
   }
   output <- list(Coeff = Coeff, Index = Index)
   return(output)
}


cow <- function(T,X,Seg,Slack, Options)
{
  # T (1xnT) target vector
  # X (mP x nP) matrix with data for mP row vectors of length nP to be warped
  # Seg (1x1) segment length; number of segments N= floor(nP/m)
  #     or (2 x N+1) matrix with segment (pre-determined) boundary-points
  # slack (1x1) 'slack' - maximum range or degree of warping in segment length "m"
  # Options (1x5) 1: triggers plot and progress-text (note: only last row/object in
  #                  "xP" is plotted)
  #               2: correlation power (minimum 1th power, maximum is 4th power)
  #               3: force equal segment lengths in "xt" and "xP" will generate
  #                    on error)
  #               4: fix maximum correction to + or - options(4) points from the
  #                  diagonal
  #               5: save in "diagnos" the table with the optimal values of loss
  #                  function and procedure (memory consuming for large problems -
  #                  on how to read tables are in the m-file)
  #         default[0 1 0 0 0](no plot; power 1; no forced equal segment lengths;
  #                 no band constraints; no Table in "diagnos")


# Check input values ------------------------------------------------------
nargin <- nargs()
  



# Initialise --------------------------------------------------------------
if(is.matrix(X))
{
  nX = nrow(X) # number of signals to be aligned
  pX = ncol(X) # number of data points in each signal
} else if(is.vector(X))
{
  nX = 1 # number of signals to be aligned
  pX = length(X) # number of data points in each signal
  
}

pT = length(T) # number of points in the target

# Xwarped intialised matrix of warped signals
Xwarped = matrix(0, nrow = nX, ncol = pT)
Time = 1    #  Time : processing time


# Initialise segments -----------------------------------------------------
Seg = round(Seg)              # Only integers are currently allowed as segment boundaries
Pred_Bound = length(Seg) > 1  # True if segment boundaries are predefined
if (Pred_Bound)
  {
    if(!(setequal(Seg[,1],rep(1,2)) & setequal(Seg[,ncol(Seg)],c(pT,pX)) ) )
        stop("End points must be equal to 1 and to the length of the pattern/target")
  
    LenSeg = t(apply(Seg,1,diff,1)) # LengSeg[1,] length of the segment in the - 1
    if(!all(LenSeg>=2))
      stop("Segments must contain at least two points")
    
    if(is.matrix(Seg)) nSeg = ncol(Seg) # number of segments
    else nSeg = length(Seg)

}else {
  
  if(Seg > min(pX,pT)) 
    stop("Segment length is larger than length of the signal")    
  
  
  if(Options[3])
  {
    nSeg =  floor( (pT-1) / Seg)
    LenSeg <- matrix(0,nrow = 2,ncol = nSeg)
    LenSeg[1,1:nSeg]= floor( (pT-1)/ nSeg)
    LenSeg[2,1:nSeg] = floor( (pX-1)/ nSeg)
    #cat("\n Segment length adjusted to cover the remainders\n")
  } else {
    nSeg = floor( (pT-1)/ (Seg - 1))
    LenSeg[1:2,1:nSeg] = Seg - 1
    if(floor( (pX-1)/ (Seg - 1) )!= nSeg )
      stop("For non-fixed segment lengths the target and the 
           signal do not have the same number of segments (try Options(3))")
  }
  
  temp = (pX-1) %% LenSeg[1,1] # The remainders are attached to the last segment
                               # in the target and in the reference
  
  if(temp > 0)
  {
    LenSeg(1,Seg)= LenSeg(1,nSeg) + temp
    if(Options[1])
      cat("\n Segments: ",LengSeg[1,1]+1, " points x ", nSeg-1, " segments +",
          LenSeg[1,ncol(LenSeg)] + 1)
  } else {
    
    if (Options[1])
      cat("\n Segments: ",LengSeg[2,1]+1,"points x",nSeg, " segments (target)")
    
  }
  temp = (pX-1) %% LenSeg[2,1]  
  if(temp> 0)
  {
    LenSeg[2,nSeg] = LenSeg[2,nSeg] + temp ;
    if (Options[1])
      cat("\n",LenSeg[2,1]+1, " points x ",Seg - 1, " segments + ",
          LenSeg[2,ncol(LenSeg)] + 1, " sigmals \n")
  } else {
    if (Options[1])
      cat("\n ", LenSeg[2,1] + 1, " points x ",nSeg," segments (signals)\n")
  }
  
}

if( any(LenSeg <= Slack + 2)) # Two points are the minimum required for linear interpolation
      stop("The slack cannot be longer than the lengthof the segments")
  
bT = cumsum(c(1,LenSeg[1,]))
bP = cumsum(c(1,LenSeg[2,]))
Warping = matrix(0, nrow = nX,ncol = nSeg + 1)

# Check Slack

if(length(Slack) > 1) # Different slacks for the segment boundaries will be implemented
  {
    if (ncol(Slack)<= nSeg)  
      stop("The number of slack parameter is not equal to the number of optimised segments")
  stop("\n Multiple slacks have not been implemented yet")
  }


Slacks_vec = -Slack:Slack # All possible slacks for a segment boundary


# Set feasible points for boundaries --------------------------------------

Bounds = matrix(1,nrow = 2, ncol = nSeg + 1)

# Slope constrains
offs = Slack * c(-1,1) %*% t(c(0:nSeg))
# all maximum and minimum values for nodes
Bounds_a = matrix(rbind(bP,bP),nrow = 2,ncol=nSeg+1,byrow = T) + offs
Bounds_b = matrix(rbind(bP,bP),nrow = 2,ncol=nSeg+1,byrow = T) + offs[,seq((nSeg+1),1,-1)]


Bounds[1,] = sapply(1:(nSeg+1),function(i) max(Bounds_a[1,i],Bounds_b[1,i]))
Bounds[2,] = sapply(1:(nSeg+1),function(i) min(Bounds_a[2,i],Bounds_b[2,i]))

# Band Constrains
if (Options[4])
  if(abs(pT-pX) > Options[4])
    stop("The band is too narrow and proper correction is not possible")

# Calculate first derivatives for interpolation
if(is.matrix(X))
{
  Xdiff = t(apply(X,1,diff,1))
  
}else if(is.vector(X))  Xdiff =diff(X)

  
# Calculate coefficients and indexes for interpolation
Int_Coeff = list()
Int_Index = list()

if(!Pred_Bound)
{
  aux1 <- list()
  aux2 <- list()
  aux3 <- list()
  aux1 = InterpCoeff(LenSeg[1,1]+1, LenSeg[2,1] + Slacks_vec + 1, Slacks_vec)
  for(i in 1:(nSeg-1))
  {
    Int_Coeff[[i]] = aux1$Coeff
    Int_Index[[i]] = aux1$Index
  }
  aux2 = InterpCoeff(LenSeg[1,nSeg]+1, LenSeg[2,nSeg] + Slacks_vec + 1, Slacks_vec)
  Int_Coeff[[nSeg]] = aux2$Coeff
  Int_Index[[nSeg]] = aux2$Index
} else {
  for (i_seg in 1:nSeg)
  {
      aux1 = InterpCoeff(LenSeg[1,i_seg] + 1, LenSeg(2,i_seg)+Slacks_vec+1,Slacks_vec) 
      Int_Coeff[[i_seg]] = aux1$Coeff
      Int_Index[[i_seg]] = aux1$Index
  }
}
  
# Dynamic Programming section ---------------------------------------------

# Indexes for the first node (boundary point) of each segment in Table
Table_Index = cumsum(c(0,diff(Bounds)+1))
# each column refer to a node
# (1,i) position of the boundary point in the signal
# (2,i) optimal
# value of the loss function up to node (i)
# (3,i) pointer to optimal preceding node (in Table)

Table = array(0, dim = c(3,Table_Index[nSeg+2],nX) )

# All loss function values apart from mode(1) are set to -Inf
Table[2,2:Table_Index[nSeg+2],] <- - Inf
v=list()
# Initialise Table
for (i_seg in 1:(nSeg+1))
{
    v  = (Bounds[1,i_seg]:Bounds[2,i_seg])
   Table[1,(Table_Index[i_seg] + 1):Table_Index[i_seg + 1],] = v

}
tic("Time dinamyc programming part")
# Forward phase
# Loop over segments
Int_Index_Seg = list()
Int_Coeff_Seg = list()
for (i_seg in 1:nSeg)
{
  # a,b,c auxiliary variables that depends on segment number and not node
  a = Slacks_vec + LenSeg[2,i_seg]
  b = Table_Index[i_seg] + 1 - Bounds[1,i_seg]
  c = LenSeg[1,i_seg] + 1
  # Counter for local table for segment i_seg
  Count = 1
  # Last node for segment i_seg
  Node_Z = Table_Index[i_seg+2]
  # First node for sement i_seg
  Node_A = Table_Index[i_seg+1]+1
  # Initialize local table for boundary
  Bound_k_Table = array(0, dim= c(2, Node_Z - Node_A + 1, nX) )
  
  # Indexes for interpolation of segment i_seg
  Int_Index_Seg = t(Int_Index[[i_seg]]) - (LenSeg[2,i_seg]+1) 
  # Coefficients for interpolation of segment i_seg
  Int_Coeff_Seg = t(Int_Coeff[[i_seg]])
  
  
  # Segment i_seg of target T
  TSeg = T[bT[i_seg]:bT[i_seg + 1]]
  TSeg_centred = TSeg - sum(TSeg)/length(TSeg)
  # norm of the vector
  Norm_TSeg_cen = sqrt(sum(TSeg_centred^2))
  
  # Loop over nodes (possible boundary positions) for 
  # segment i_seg
  for (i_node in Node_A:Node_Z)
  {
    # Possible predecessors given the allowed segment lengths
    Prec_Nodes = Table[1,i_node,]-a
    # Arcs allowed by local and global constraints
    Allowed_Arcs = Prec_Nodes>=Bounds[1,i_seg] & Prec_Nodes <= Bounds[2,i_seg]
    # Pointer to predecesors in Table
    Nodes_TablePointer = b + Prec_Nodes[Allowed_Arcs]
    # Number of allowed arcs
    N_AA = sum(Allowed_Arcs)
    # Sometimes boundaries are ineffective and few nodes are allowed that cannot
    # be reached. It has to be further investigated
    
    if (N_AA)
    {
      # Interpolation signal indexes for all allowed arcs for node i_node
      Index_Node = Table[1,i_node,] + Int_Index_Seg[,Allowed_Arcs]
      # Interpolation coefficients for all the allowed arcs for node i_node
      Coeff_b = Int_Coeff_Seg[,Allowed_Arcs]
      Coeff_b = as.vector(Coeff_b) # Check when there are more than query signal
      Coeff_b0 = matrix(rep(Coeff_b,nX))
      Coeff_b1 = matrix(Coeff_b0,nrow = nX,byrow = TRUE)
      Xi_Seg = X[Index_Node]
      Xi_diff = Xdiff[Index_Node]
      # Interpolate for all allowed predecessors
      Xi_Seg1 = array(Xi_Seg + (Coeff_b1 * Xi_diff), dim = c(c,N_AA*nX))
      # Mean of the interpolate segments
      if(N_AA == 1){
        Xi_Seg_mean = sum(Xi_Seg1)/length(Xi_Seg1)
      }else{
        Xi_Seg_mean = apply(Xi_Seg1,2,sum)/nrow(Xi_Seg1)
      }

      # Fast method for calculating the covariance of T and X
      # no centering of X is needed
      if(is.vector(Xi_Seg1))
      {
        Norm_Xi_Seg_cen = sqrt(sum(Xi_Seg1^2)-length(Xi_Seg1)*Xi_Seg_mean^2 )
      }else {
        Norm_Xi_Seg_cen = sqrt(apply(Xi_Seg1^2,2,sum)-nrow(Xi_Seg1)*Xi_Seg_mean^2 )
      }
      # Correlation coefficients relative to all probable predecessors
      CCs_Node = (t(TSeg_centred) %*% Xi_Seg1)/(Norm_TSeg_cen * Norm_Xi_Seg_cen)
      # If standard deviation is zero, update is not chosen
      CCs_Node[!is.finite(CCs_Node)] = 0
      CCs_Node1 = array(CCs_Node,dim=c(N_AA,nX)) 
      
      if (Options[2]==1)
      {
        Cost_Fun = array(Table[2,Nodes_TablePointer,],dim = c(N_AA,nX)) + CCs_Node1
      } else {
        Cost_Fun = array(Table[2,Nodes_TablePointer,],dim = c(N_AA,nX)) + CCs_Node1^Options[2]
      }
      ind = apply(Cost_Fun,2,max)
      pos = which(Cost_Fun == ind)
      Bound_k_Table[1,Count,] = ind
      Bound_k_Table[2,Count,] = Nodes_TablePointer[pos]
      Count = Count + 1
    }
  }# end-for i_node
  
  
  # Update general table (it turned out to be faster than using
  # Table directly in the loop over nodes)
  Table[2:3,Node_A:Node_Z,] = Bound_k_Table
  
} # end-for i_seg

Time = toc()

# Loop over samples/signals
for (i_sam in 1:nX)
{
  # Backward phase
  # Backtrace optimal boundaries using pointers in Table
  Pointer = ncol(Table)
  Warping[i_sam, (nSeg+1)] = pX
  for(i_bound in nSeg:1)
  {
    Pointer = Table[3,Pointer,i_sam]
    Warping[i_sam,i_bound] = Table[1,Pointer,i_sam]
  }
}

# Reconstruct aligned signals
for (i_seg in 1:nSeg)
{
  indT = bT[i_seg]:bT[i_seg + 1]
  indT = bT(i_seg + 1)-bT(i_seg)
  for (i_sam in 1:nX)
  {
    indX = Warping[i_sam,i_seg]:Warping[i_sam,i_Seg+1]
    ienX = Warping[i_sam,i_seg + 1]-Warping[i_sam.i_seg]
#    Xwarped[i_sam, indT] = approx(t(indX) - Warping[i_sam,i_seg]+1, 
#                                  t(X[i_sam,indX])) 
  }
}

 return(list(Warping = Warping)) 
}  
  























gc <- setRefClass("gc",fields = list(peak = "numeric",
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
                      for(j in 1:np)
                      {
                        if(peak[j]-window[j]>0)
                        {
                          ini<- peak[j]-window[j]
                          fin<- peak[j]+window[j]
                          # obtain position to be changed
                          positionToChange <- seq(ini,fin,1)
                          
                          x[ini:fin]<-x[ini:fin] + sapply(positionToChange,function(pos)
                          {
                            # generate a value from normal dist. with mean equal to
                            # intensity
                            abs(rnorm(1,intensity[j],sqrt(variance[j])))/exp(abs(peak[j]-pos)/2)
                            # exponential decay when the retention time differ from 
                            # the position of the peak
                            
                          })                  
                          
                        }
                      }
                      print = x 
                    },
                    range = function()
                    {
                      range = seq(lowerlimit,upperlimit,1)
                    }
                  ))



set.seed(1)
v1 <- c(50,110,140,175,190,230,250)
u <- round(rnorm(7,5,2),0)
u
v2 <- v1 + u
v2

gc1<-gc(peak =v1,
        variance =c(2,2,1,1,1,1,1),
        intensity = c(10,35,55,35,180,700,400),
        window = c(3,3,3,5,5,10,10),
        lowerlimit = 0,
        upperlimit = 500)

# signal shifted by adding a random normal vector to the peak vectors
gc2<-gc(peak =v2,
        variance =c(2,2,1,1,1,1,1),
        intensity = c(10,35,55,35,180,700,400),
        window = c(3,3,3,5,5,10,10),
        lowerlimit = 0,
        upperlimit = 500)

B1 <- gc1$print()
# signal shifted
T1 <- gc2$print()

plot(gc1$range(), B1, type = "l", col = "blue")
lines(gc1$range(),T1, col = "green")





Pprime <- function(P,T1,xsw,xew)
{
  output<-list()
  
  if(length(which(P$x == xsw)) > 0 )
      xi <-which(P$x == xsw)
  else if(length(which(P$x < xsw) > 0) )
      xpi<- which(P$x < xsw)
  else if(length(which(P$x > xsw) > 0 ))
      xpi<-which(P$x > xsw)
  
  if(length(which(P$x == xew)) > 0 )
      xi1<- which(P$x == xew)
  else if(length(which(P$x < xew)) > 0)
      xpi1<-which(P$x < xew)
  else if(length(which(P$x > xew)) > 0 )
      xpi1<-which(P$x > xew) 
    
  nT <- length(T1$x[xi:xi1])

 
  pj <- rep(o,nT)
  pj[1] <- P$x[xpi]
  pj[nT] <- P$x[xpi1]

  fpj<- rep(0,nT)
  fpj[1] <- P$f[xpi]
  fpj[nT] <- P$f[xpi1]
  
  
  
  for (j in 2:(nT-1))
    pj[j] <-P$x[xpi] + j * (P5x[xpi1] - P$x[xpi])/(xew - xsw)

  for (j in 2:(nT-1))
  {
      lb = which(P$x < pj[j])
      ub = which(P$x > pj[j])
      xs = P$x[lb]
      xe = P$x[ub]
      fxs = P$f[lb]
      fxe = P$f[ub]
      m = (fxe-fxs)/(xe-xs)
      c = fxe-m*xe      
      fpj[j] = m*pj[j] + c
      
  }
  
  output <- list(x = pj, f = fpj)
  return(output)
  
}


# benefit function 
f <- function(P,T1,xsw,xew)
{
  # xsw : start warping point
  # xew : end warping point
  # xi : warping start position
  # xi1 : warping end position
  
  xi <- which(P$x == xsw)
  xi1 <- which(T1$x == xew)
  
  temp <- list()
  nTi <- length(T1$x[xi:xi1])
  PPf <- list(x = rep(0.0,nTi), f = rep(0.0,nTi) )
  
  # Use the warping points to find the unwarped points and obtain pj's
  # calculate the Pprime using pj's in order to compare with the target
  # Calculate correlation between Pprime and T in a segment.
  # if Pw and T1 are the same length
  #ini<- which(T$x == xi)
  #fin <-which(T$x == xi1)
  
  if (all(T1$x[xi:xi1] == P$x[xi:xi1]) & length(T1$x[xi:xi1]) == length( P$x[xi:xi1] ) )
  {
      b <- P$f[xi:xi1]
  } else
  {
  
      PPf <-Pprime(P,T1,xsw,xse)
      b <- PPf$f
  }
    a <- T1$f[xi:xi1]
    output <- cov(a,b)
    return(output)
  
}

align <- function(P,Target,m,t)
{
  # P:  class GC signal that contains the retention time and intensity
  # Target: class Gc signal that contains the target retention and intensity
  # m:      Number of segments
  # t:       
  # F1: matrix containing the cumulated benefit function
  
  # Pre-aligming length of chromatogram
  # Intervals of 1 unit in P
  Lp<- max(P$x)-min(P$x)
  
  # Post-aligning length of chromatogram and length of 
  # target chromatogram
  Lt <- max(Target$x)-min(Target$x)
  
  # Calculate number of sections for P
  N <- floor(Lp/m)
  
  # calculate difference in mean section length between P and T
  d <- floor(Lt/N) - m
  
  F1 <- matrix(0.0, nrow = (N+1), ncol = (Lt+1))
  U <- matrix(0.0, nrow = (N+1), ncol = (Lt+1))
  
  for (i in 0:N)
  {
    for(x in 0:Lt)
    {
      F1[i,x] <- -Inf  
    }
  }
  
  F1[N,0] = 0
  
  for (i in (N-1): 0)
  {
    xstart <- max(i*(m+d-t),Lt-(N-i)*(m+d+t))
    xend <- min(i*(m+d + t),Lt-(N-i)*(m+d-t))
    cat("\n \\(", xstart, "-" ,xend, "\\) \n")
#    for (x in xstart:xend)
#    {
#      for (u in (d-t):(d+t))
#      {
#        if(F1[i+1,x+m+u] == -Inf)
#            fsum <- f(x,x+m+u)
#        else  fsum <- F1[i+1,x+m+u] + f(x,x+m+u)
#        if (fsum > F1[i,x]) then
#        F1[i,x] <- fsum
#        U[i,x] <- u
#      }
#    }
  }
  
  x(0) = 0
  for(i in 0:N-1)
  {
    u[i] <- U[i,x[i]]
    x[i+1] <- x[i]+m+u(i)
  }
  return(x = x, u=u,F1=F1,U=U)
}
