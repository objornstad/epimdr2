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
  outv<-as.data.frame(ode(xstrt,steps,sirvmod,par))
  fsv<-max(outv$K)

  par<-c(B=beta, r=1/LP, g = 1/IP, q = vaccine_efficacy,
             P = target_vaccination, Dt = intervention_length, T = day)
outi<-as.data.frame(ode(xstrt,steps,sirvmod,par))
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
sirvmod<-function(t,x,parms){
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

