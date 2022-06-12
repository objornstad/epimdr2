#' Function to simulate the stochastic TSIR
#'
#' Function to simulate the stochastic TSIR assuming stochasticity in transmission and a Poisson birth-death process
#'
#' @param alpha the exponent on I
#' @param B the birth rate
#' @param beta the transmission rate
#' @param sdbeta the standard deviation on beta
#' @param S0 the initial susceptible fraction
#' @param I0 the initial number of infecteds
#' @param IT the length of simulation
#' @param N the population size
#' @return A list with time series of simulated infected and susceptible hosts
#' @examples
#' out = tsirSim()
#' @export
#' @importFrom stats rnorm
#' @importFrom stats rpois
tsirSim=function(alpha=0.97, B=2300, beta=25, sdbeta=0,
    S0 = 0.06, I0=180, IT=520, N=3.3E6){
    lambda = rep(NA, IT)
    I = rep(NA, IT)
    S = rep(NA, IT)
    I[1] = I0
    lambda[1] = I0
    S[1] = S0*N
    for(i in 2:IT) {
        lambda[i] = rnorm(1, mean=beta, sd=sdbeta) * I[i - 1]^alpha * S[i - 1] /N
        if(lambda[i]<0) {lambda[i]=0}
        I[i] = rpois(1, lambda[i])
        S[i] = S[i - 1] + B - I[i]
    }
    list(I = I, S = S)
}

#' Function to simulate the seasonally-forced TSIR
#'
#' Function to simulate the stochastic TSIR assuming stochasticity in transmission and a Poisson birth-death process
#'
#' @param beta the seasonal transmission coefficients
#' @param alpha the exponent on I
#' @param B a vector of Births (the length of which determines the length of the simulation)
#' @param N the population size
#' @param inits a list containing initial S and I
#' @param type an argument "det" or "stoc" that determines whether a deterministic or stochastic simulation is done
#' @return A list with time series of simulated infected and susceptible hosts
#' @examples
#' \dontrun{see chapter 8 in book}
#' @export
tsirSim2=function(beta, alpha, B, N,  inits = list(Snull = 0, Inull = 0), type = "det"){
    type = charmatch(type, c("det", "stoc"), nomatch = NA)
    if(is.na(type))
        stop("method should be \"det\", \"stoc\"")
        
    IT = length(B)
    s = length(beta)
    lambda = rep(NA, IT)  
    I = rep(NA, IT)
    S = rep(NA, IT)
    
    I[1] = inits$Inull
    lambda[1] = inits$Inull
    S[1] = inits$Snull
    
    for(i in 2:IT) {
        lambda[i] = beta[((i - 2) %% s) + 1] * S[i - 1] * (I[i - 1]^alpha)/N
        if(type == 2) {
                I[i] = rpois(1, lambda[i])
            }
        if(type == 1) {
            I[i] = lambda[i]
        }
        S[i] =S[i - 1] + B[i] - I[i]
    }
    return(list(I = I, S = S))
}
#' Function to do  Lyapunov exponent calculations from a TSIR simulation
#'
#' Function to do  Lyapunov exponent calculations from a TSIR simulation
#'
#' @param I a vector containing the time series of Is
#' @param S vector containing the time series of Ss
#' @param bt the seasonal transmission coefficients
#' @param alpha the exponent on I
#' @param N the population size
#' @return An object of class lyap with the Lyapunov exponent, values for the Jacobians, parameters and data
#' @examples
#' \dontrun{see chapter 10 in book}
#' @export
tsirLyap=function(I, S, alpha, bt, N){
  IT <- length(I)
  s <- length(bt)
  j11=rep(NA, IT)
  j12=rep(NA, IT)
  j21=rep(NA, IT)
  j22=rep(NA, IT)
  #initial unit vector
  J=matrix(c(1,0),ncol=1)
  #loop over the attractor
  for(i in 1:IT) {
  j11=1 -  bt[((i - 1) %% s) + 1] * I^alpha/N
  j12=-( bt[((i - 1) %% s) + 1] * S * (I^(alpha - 1) * alpha)/N)
  j21= bt[((i - 1) %% s) + 1] * I^alpha/N
  j22= bt[((i - 1) %% s) + 1] * S * (I^(alpha - 1) * alpha)/N
    J<-matrix(c(j11[i],j12[i],j21[i],j22[i]), ncol=2, byrow=TRUE)%*%J
      }
res=list(lyap=log(norm(J))/IT, j11=j11, j12=j12, j21=j21, j22=j22, I=I, S=S, alpha=alpha, bt=bt, N=N)
class(res)="lyap"
return(res)
}

#' Function to calculate the local Lyapunov exponents for the TSIR 
#'
#' Function to calculate the local Lyapunov exponents from an object of class \code{lyap}.
#'
#' @param x an object of class \code{lyap} (normally from a call to \code{tsirLyap})
#' @param m number of forward iterations on the attractor
#' @return An object of class llyap with the local Lyapunov exponent and S-I data
#' @examples
#' \dontrun{see chapter 10 in book}
#' @export
tsirLlyap=function(x, m=1){
llyap=rep(NA, length(x$I))
for(i in 1:(length(x$I)-m)){
J=matrix(c(1,0,0,1), ncol=2)
for(k in 0:(m-1)){J = matrix(c(x$j11[(i+k)], x$j12[(i+k)], x$j21[(i+k)], x$j22[(i+k)]), ncol = 2, byrow=TRUE)%*%J}
llyap[i]=log(max(abs(eigen(J)$values)))/m
}
res=list(llyap=llyap, I=x$I, S=x$S)
class(res)="llyap"
return(res)
}

 #' Function to simulate the spatially-extended seasonally-forced TSIR
#'
#' Function to simulate the spatially-extended seasonally-forced TSIR for a patchily
#' distributed host population. Coupling is assumed to be global and according to a 
#' commuter model so with p patches and a coupling of c, local transmission
#' is reduced by a fraction (1-c*p).
#'
#' @param beta the seasonal transmission coefficients
#' @param alpha the exponent on I
#' @param B a vector of Births (the length of which determines the length of the simulation)
#' @param N the population size
#' @param p the number of patches
#' @param c the the spatial coupling
#' @param inits a list containing a vector of initial S and a vector of initial I for each patch
#' @param type an argument "det" or "stoc" that determines whether a deterministic or stochastic simulation is done
#' @return A list with time series of simulated infected and susceptible hosts
#' @examples
#' \dontrun{see chapter XX in book}
#' @export
tsirSpat=function(beta, alpha, B, N, p, c, inits, type = "det"){
  type = charmatch(type, c("det", "stoc"), nomatch = NA)
  if(is.na(type))
    stop("method should be \"det\", \"stoc\"")
  IT = dim(B)[1]
  s = length(beta)
  lambda = matrix(NA, nrow=IT, ncol=p)  
  I = matrix(NA, nrow=IT, ncol=p)
  S = matrix(NA, nrow=IT, ncol=p)
  
  I[1,] = inits$Inull
  lambda[1,] = inits$Inull
  S[1,] = inits$Snull
  cmat=matrix(c, ncol=p, nrow=p)
  #diag(cmat)=1  
  diag(cmat)=1 - c*p 
  for(i in 2:IT) {
    lambda = beta[((i - 2) %% s) + 1]*cmat%*%(I[i - 1,]^alpha)*S[i - 1,]/N 
    if(type == 2) {
      I[i,] = rpois(p, lambda)
    }
    if(type == 1) {
      I[i,] = lambda
    }
    S[i,] =S[i - 1,] + B[i,] - I[i,]
  }
  return(list(I = I, S = S))
}

