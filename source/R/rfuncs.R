#c2

#' Negative log-likelihood function for the chain-binomial model
#' @param S0 a scalar with value for S0
#' @param beta a scalar with value for beta
#' @param I a vector incidence aggregated at serial interval
#' @return the negative log-likelihood for the model 
#' @examples
#' twoweek=rep(1:15, each=2)
#' niamey_cases1=sapply(split(niamey$cases_1[1:30], twoweek), sum)
#' llik.cb(S0=6500, beta=23, I=niamey_cases1)
#' @export
#' @importFrom stats dbinom
llik.cb = function(S0,beta,I){
    n = length(I)
    S = floor(S0-cumsum(I[-n]))
    p = 1-exp(-beta*(I[-n])/S0)  
    L = -sum(dbinom(I[-1],S,p,log=TRUE))  
    return(L)
}


#' Function to simulate the chain-binomial model
#' @param S0 a scalar with value for S0
#' @param beta a scalar with value for beta
#' @return A data-frame with time series of susceptibles and infecteds
#' @examples
#' sim=sim.cb(S0=6500, beta=23)
#' @export
#' @importFrom stats rbinom
sim.cb=function(S0, beta){
I=1
S=S0
i=1
while(!any(I==0)){
i=i+1
I[i]=rbinom(1, size=S[i-1], prob=1-exp(-beta*I[i-1]/S0))
S[i]=S[i-1]-I[i]
}
out=data.frame(S=S, I=I)
return(out)
}

#c3

#' Auxilliary function used by llik.pc 
#' @param a a vector with the ages 
#' @param up a vector with upper age-bracket cut-offs
#' @param foi a vector with FoI
#' @return A vector with FoIs matched to data
#' @seealso{llik.pc}
#' @export
integrandpc=function(a, up, foi){
wh=findInterval(a, sort(c(0,up)))
dur=diff(sort(c(0,up)))
inte=ifelse(wh==1, foi[1]*a, sum(foi[1:(wh-1)]*dur[1:(wh-1)])+foi[wh]*(a-up[wh-1]))
return(inte)
 }


#' Function to estimate parameters for the picewise-constant catalytic model
#'
#' This function uses binomial likelihoods to estimate the picewise-constant FoI model from age-incidence data
#'
#' @param par a vector with initial guesses 
#' @param age a vector with the ages 
#' @param num a vector with number infected by age
#' @param denom a vector with number tested by age
#' @param up a vector with upper age-bracket cut-offs
#' @return The negative log-likelihood for a candidate piecewise constant catalytic model
#' @examples
#' x=c(1,4,8,12,18,24)
#' para=rep(.1,length(x))
#' \dontrun{optim(par=log(para),fn=loglikpc, age=rabbit$a, num=rabbit$inf, denom=rabbit$n, up=x)}
#' @export
#' @importFrom stats integrate
llik.pc = function(par, age, num, denom, up) {
ll = 0
for (i in 1:length(age)) {
p = 1 - exp(-integrandpc(a=age[i], up = up, foi = exp(par)))
ll = ll + dbinom(num[i], denom[i], p, log = T)
}
return(-ll)
}

#logspline not documented

#c4

#c6


#c7 

#' Gillespie exact algorithm
#'
#' Function simulating a dynamical system using the Gillespie exact algorithm
#'
#' @param rateqs a list with rate equations 
#' @param eventmatrix a matrix of changes in state variables associated with each event
#' @param parameters a vector of parameter values
#' @param initialvals a vector of initial values for the states
#' @param numevents number of events to be simulated
#' @return A data frame with simulated time series
#' @examples
#' rlist=c(quote(mu * (S+I+R)), quote(mu * S), quote(beta * S * I /(S+I+R)), 
#'  quote(mu * I), quote(gamma * I), quote(mu*R))
#' emat=matrix(c(1,0,0,-1,0,0,-1,1,0,0,-1,0,0,-1,1,0,0,-1),ncol=3, byrow=TRUE)
#' paras  = c(mu = 1, beta =  1000, gamma = 365/20)
#' inits = c(S=100, I=2, R=0)
#' sim=gillespie(rlist, emat, paras, inits, 100)
#' @export
#' @importFrom stats rexp
gillespie=function(rateqs, eventmatrix, parameters, initialvals, numevents){
res=data.frame(matrix(NA, ncol=length(initialvals)+1, nrow=numevents+1))
names(res)=c("time", names(initialvals))
res[1,]=c(0, initialvals)
for(i in 1:numevents){
rat=sapply(rateqs, eval, as.list(c(parameters, res[i,])))
res[i+1,1]=res[i,1]+rexp(1, sum(rat))
whichevent=sample(1:nrow(eventmatrix), 1, prob=rat)
res[i+1,-1]=res[i,-1]+eventmatrix[whichevent,]
}
return(res)
}


#' Gillespie tau-leap algorithm
#'
#' Function simulating a dynamical system using the Gillespie tau-leap approximation
#'
#' @param rateqs a list with rate equations 
#' @param eventmatrix a matrix of changes in state variables associated with each event
#' @param parameters a vector of parameter values
#' @param initialvals a vector of initial values for the states
#' @param deltaT the tau-leap time interval
#' @param endT the time length of simulation
#' @return A data frame with simulated time series
#' @examples
#' rlist2=c(quote(mu * (S+E+I+R)), quote(mu * S), quote(beta * S * I/(S+E+I+R)), 
#'  quote(mu*E), quote(sigma * E), quote(mu * I), quote(gamma * I), quote(mu*R))
#' emat2=matrix(c(1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1,0,0,0,-1,1,0,0,0,-1),
#' ncol=4, byrow=TRUE)
#' paras  = c(mu = 1, beta =  1000, sigma = 365/8, gamma = 365/5)
#' inits = c(S=999, E=0, I=1, R = 0)
#' sim2=tau(rlist2, emat2, paras, inits, 1/365, 1)
#' @export
tau=function(rateqs, eventmatrix, parameters, initialvals, deltaT, endT){
time=seq(0, endT, by=deltaT)
res=data.frame(matrix(NA, ncol=length(initialvals)+1, nrow=length(time)))
res[,1]=time
names(res)=c("time", names(initialvals))
res[1,]=c(0, initialvals)
for(i in 1:(length(time)-1)){
rat=sapply(rateqs, eval, as.list(c(parameters, res[i,])))
evts=rpois(1,  sum(rat)*deltaT)
if(evts>0){
whichevent=sample(1:nrow(eventmatrix), evts, prob=rat, replace=TRUE)
mt=rbind(eventmatrix[whichevent,], t(matrix(res[i,-1])))
#strange fix:
mt=matrix(as.numeric(mt), ncol=ncol(mt))
res[i+1,-1]=apply(mt,2,sum)
res[i+1, ][res[i+1,]<0]=0
}
else{
res[i+1,-1]=res[i,-1]
}}
return(res)
}

#' Function to predict efficacy of outbreak-response vaccination campaign
#' @param R reproductive ratio
#' @param day first day of ORV campaign
#' @param vaccine_efficacy Vaccine efficacy
#' @param target_vaccination fraction of population vaccinated during ORV campaign
#' @param intervention_length duration of ORV campaign
#' @param mtime length of simulation
#' @param LP length of latent period
#' @param IP length of infectious period
#' @param N initial susceptible population size
#' @return A list of gradients 
#' @examples
#' red1=retrospec(R=1.8, 161, vaccine_efficacy=0.85, target_vaccination=0.5, 
#'  intervention_length=10, mtime=250, LP=8, IP=5, N=16000)
#' 1-red1$redn
#' @export
retrospec<-function(R,day, vaccine_efficacy,target_vaccination,intervention_length, mtime, LP=7, IP=7, N=10000){
  steps<-1:mtime
        out<-matrix(NA,nrow=mtime, ncol=3)
  xstrt<-c(S=1-1/N,E=0,I=1/N,R=0,K=0)   #starting values
  beta<- R/IP         #transmission rate
  par<-c(B=beta, r=1/LP, g = 1/IP, q = vaccine_efficacy,
                P = 0, Dt = 0, T = Inf)
  outv<-as.data.frame(ode(xstrt,steps,sivmod,par))
  fsv<-max(outv$K)

  par<-c(B=beta, r=1/LP, g = 1/IP, q = vaccine_efficacy,
             P = target_vaccination, Dt = intervention_length, T = day)
outi<-as.data.frame(ode(xstrt,steps,sivmod,par))
  fsi<-max(outi$K)

         res<-list(redn=fsi/fsv)
  return(res)
}


#' Gradient-function for the SIR model with outbreak-response vaccination
#' @param t Implicit argument for time
#' @param x  A vector with values for the states
#' @param parms A vector with parameter values for the SIR system
#' @return A list of gradients 
#' @seealso \code{\link{retrospec}}
#' @export
sivmod<-function(t,x,parms){
    S<-x[1]
    E<-x[2]
    I<-x[3]
    R<-x[4]
    K<-x[5]
    with(as.list(parms),{
      Q<- ifelse(t<T | t>T+Dt,0,(-log(1-P)/Dt))
      dS<- -B*S*I-q*Q*S
      dE<- B*S*I-r*E
      dI<- r*E - g*I
      dR<- g*I+q*Q*S
      dK<-r*E
      res<-c(dS,dE,dI,dR,dK)
      list(res)
    })
  }



#c10


#' Gradient-function for the SIRWS model
#' @param t Implicit argument for time
#' @param logy  A vector with values for the log(states)
#' @param parms A vector with parameter values for the SIRWS system
#' @return A list of gradients (in log-coordinates)
#' @examples
#' require(deSolve)
#' times  = seq(0, 26, by=1/10)
#' paras  = c(mu = 1/70, p=0.2, N = 1, beta = 200, omega = 1/10, gamma = 17, kappa=30)
#' start = log(c(S=0.06, I=0.01, R=0.92, W = 0.01))
#' out = as.data.frame(ode(start, times, sirwmod, paras))
#' @export
sirwmod=function(t, logy, parms){
  y=exp(logy)
   S=y[1]
   I=y[2]
   R=y[3]
   W=y[4]

   with(as.list(parms),{
   dS = mu * (1-p) * N  - mu * S  - beta * S * I / N + 2*omega * W
   dI = beta * S * I / N - (mu + gamma) * I
   dR = gamma * I - mu * R - 2*omega * R +  kappa * beta * W * I / N + mu*p*N
   dW = 2*omega * R - kappa * beta * W * I / N - (2*omega +mu)* W 
   res=c(dS/S, dI/I, dR/R, dW/W)
   list(res)
 })
 }


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
#' sim1=networkSir(cm1,.3,0.1)
#' summary(sim1)
#' \dontrun{plot(sim1)}
#' @aliases netSir
#' @export
networkSir=function(CM,tau,gamma){
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
#' \code{\link{netSir}}
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
#' \code{\link{netSir}}
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

#c14

#' The Nicholson-Bailey model
#'
#' Function to simulate the Nicholson-Bailey Parasitoid-host model
#'
#' @param R the host reproductive rate
#' @param a the parasitoid search efficiency
#' @param T the length of simulation (number of time-steps)
#' @param H0 initial host numbers
#' @param P0 initial parasitoid numbers
#' @return A list of simulated Host and Parasitoid numbers
#' @examples
#' sim= nbmod(R=1.1,a=0.1)
#' @export
 nbmod = function(R, a, T = 100, H0 = 10, P0 = 1){
   #T is length of simulation (number of time-steps)
   #H0 and P0 are initial numbers
   #we provide default parameters except for R and a
   H=rep(NA, T) #host series
   P=rep(NA, T) #parasitoid series

   H[1] = H0 #Initiating the host series
   P[1] = P0 #Initiating the parasitoid series

   for(t in 2:T){
     H[t] = R * H[t-1] * exp(- a * P[t-1])
     P[t] = R * H[t-1] * (1-exp(- a * P[t-1]))
   } #end of loop

   #the two vectors of results are stored in a "list"
   res= list(H = H, P = P) 
   return(res)
} 

