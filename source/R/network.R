
#c11

#' Function to generate a ring lattice
#' @param N the number of nodes
#' @param K the number of neighbors to which each node is connected so degree = 2xK
#' @return An object of class CM (contact matrix)
#' @examples
#' cm=ringlattice(N=20,K=4)
#' @export
#' @importFrom stats toeplitz
ringlattice=function(N,K){
CM=toeplitz(c(0,rep(1,K),rep(0,N-2*K-1),rep(1,K)) )
    class(CM)="cm"
    return(CM)
}

#' Function to plot an object of class CM
#' @param x an object of class cm
#' @param ... other arguments 
#' @return A plot of the contract matrix
#' @examples
#' cm=ringlattice(N=20,K=4)
#' \dontrun{plot(cm)}
#' @export
#' @importFrom graphics symbols
#' @importFrom graphics segments
plot.cm=function(x, ...){
N=dim(x)[1]
theta=seq(0,2*pi,length=N+1)
x2=cos(theta[1:N])
y2=sin(theta[1:N])
symbols(x2,y2, fg=0, circles=rep(1, N), inches=0.1, bg=1, xlab="", ylab="")
segx1=as.vector(matrix(x2, ncol=length(x2), nrow=length(x2), byrow=TRUE))
segx2=as.vector(matrix(x2, ncol=length(x2), nrow=length(x2), byrow=FALSE))
segy1=as.vector(matrix(y2, ncol=length(x2), nrow=length(x2), byrow=TRUE))
segy2=as.vector(matrix(y2, ncol=length(x2), nrow=length(x2), byrow=FALSE))
segments(segx1,segy1, segx2, segy2, lty=as.vector(x))
}

#' Function to generate a Watts-Strogatz network
#' @param N the number of nodes
#' @param K the number of neighbors to which each node is connected so degree = 2*K
#' @param Prw the rewiring probability
#' @return An object of class CM (contact matrix)
#' @examples
#' cm2=wattsStrogatz(N=20, K=4, Prw=.3)
#' @export
#' @importFrom stats runif
wattsStrogatz=function(N, K, Prw){
CM=ringlattice(N=N, K=K)
CMWS=CM
tri=CM[upper.tri(CM)]
Br=rbinom(length(tri),1,Prw)  # specify which edges to break
a=0
for(i in 1:(N-1)){                  
   for(j in (i+1):N){
  a=a+1               
  if(Br[a]==1 & CMWS[i,j]==1){ # if "break" is specified in Br matrix
    CMWS[i,j]=CMWS[j,i]=0 # break edge
    tmp=i
    tmp2=c(i, which(CMWS[i,]==1))             
    while(any(tmp2==tmp)){tmp=ceiling(N*runif(1))} # search new edge  
    CMWS[i,tmp]=CMWS[tmp,i]=1 # make new edge
    }
  }
}
class(CMWS)="cm"
return(CMWS)
}

#' Function to calculate the degree distribution for an object of class CM
#' @param object an object of class cm
#' @param plot if TRUE a bar plot of the degree distribution is produced 
#' @param ... other arguments 
#' @return A plot of the contract matrix
#' @examples
#' cm=wattsStrogatz(N=20, K=4, Prw=.3)
#' summary(cm)
#' @export
#' @importFrom graphics barplot
summary.cm=function(object, plot=FALSE, ...){
  x=table(apply(object, 2, sum))
  res=data.frame(n=x)
  names(res)=c("degree", "freq")
  if(plot) barplot(x, xlab="degree")
    return(res)
}

#' Function to generate a Barabasi-Albert network
#' @param N the number of nodes
#' @param K the number of neighbors to which each node is connected so degree = 2*K
#' @return An object of class CM (contact matrix)
#' @examples
#' cm3=barabasiAlbert(200, 4)
#' @export
barabasiAlbert=function(N, K){
#https://en.wikipedia.org/wiki/Barabasi-Albert_model#Algorithm
CM=matrix(0, ncol=N, nrow=N)
CM[1,2]=1
CM[2,1]=1

for(i in 3:N){                  
  probs=apply(CM, 1, sum)
  link=unique(sample(c(1:N)[-i], size=min(c(K, i-1)), prob=probs[-i]))
    CM[i, link]=CM[link, i]=1
}
class(CM)="cm"
return(CM)
}

#' Function to simulate an epidemic on a network
#'
#' Function to simulate a stochastic (discrete time) Reed-Frost SIR model on a social network
#' 
#' @param CM a contact matrix
#' @param tau the transmission probability
#' @param gamma the recovery probability
#' @return An object of class netSir with infectious status for each node through time
#' @examples
#' cm1=barabasiAlbert(N=200,K=2)
#' sim1=sirNetmod(cm1,.3,0.1)
#' summary(sim1)
#' \dontrun{plot(sim1)}
#' @export
sirNetmod=function(CM,tau,gamma){
#generate SIR epidemic on a network specified by the contact matrix 
#CM = contact matrix
#tau = probability of infection across an edge
#gamma = probability of removal per time step 
N=dim(CM)[1]
I=matrix(rep(0,N),nrow=N,ncol=1)   # initialize infecteds 
S=matrix(rep(1,N),nrow=N,ncol=1)  # initialize susceptibles 
R=matrix(rep(0,N),nrow=N,ncol=1)  # initialize removed 
I1=sample(1:N, size=1)
I[I1,1]=1 # initialize 1 random infected
S[I1,1]=0     

t=1
while(sum(I[,t-1])>0 | t==1){
    t=t+1
    infneigh=CM%*%I[,t-1]
    pinf=1-(1-tau)^infneigh
    newI=rbinom(N, S[,t-1], pinf)
    newR=rbinom(N, I[,t-1], gamma)

    nextS=S[,t-1]-newI
    nextI=I[,t-1]+newI-newR
    nextR=R[,t-1]+newR

    I=cbind(I, nextI)
    S=cbind(S, nextS)
    R=cbind(R, nextR)
    }
  
res=list(I=I,S=S,R=R)
class(res)="netSir"
return(res)
}

#' Function to summarize a netSir object
#' @param object an object of class netSir
#' @return A data-frame with the time series of susceptible, infected and recovered individuals
#' @param ... other arguments 
#' @seealso
#' \code{\link{sirNetmod}}
#' @export
summary.netSir=function(object, ...){
t=dim(object$S)[2]
S=apply(object$S,2,sum)
I=apply(object$I,2,sum)
R=apply(object$R,2,sum)
res=data.frame(S=S,I=I,R=R)
return(res)
}

#' Function to plot a netSir object
#' @param x an object of class netSir
#' @param ... other arguments 
#' @seealso
#' \code{\link{sirNetmod}}
#' @export
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom graphics lines
plot.netSir = function(x, ...){
   y = summary(x)	
   plot(y$S, type = "b", xlab = "time", ylab = "", ylim=range(y))
   lines(y$I, type = "b", col = "red")	
   lines(y$R, type = "b", col = "blue")	
   legend("left", legend = c("S", "I", "R"),
     lty = c(1,1,1), pch = c(1,1,1),
     col = c("black", "red", "blue"))
}

#' Function to calculate R0 from a contact matrix
#' @param CM an object of class CM
#' @param tau = probability of infection across an edge
#' @param gamma = probability of removal per time step 
#' @return the R0
#' @examples
#' cm1=barabasiAlbert(N=200,K=2)
#' r0fun(cm1, 0.3, 0.1)
#' @export
r0fun=function(CM, tau, gamma){
x=apply(CM, 2, sum)
(tau/(tau+gamma))*(mean(x^2)-(mean(x)))/mean(x)
}
