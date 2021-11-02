#c1

#' Gradient-function for the SIR model
#' @param t Implicit argument for time
#' @param y  A vector with initial values for the states
#' @param parameters A vector with parameter values for the SIR system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 26, by=1/10)
#' paras  = c(mu = 0, N = 1, beta =  2, gamma = 1/2)
#' start = c(S=0.999, I=0.001, R = 0)
#' out=ode(start, times, sirmod, paras)
#' @export
 sirmod=function(t, y, parameters){
   #Pull state variables from y vector
   S=y[1]
   I=y[2]
   R=y[3]
   #Pull parameter values from the input vector
   beta=parameters["beta"]
   mu=parameters["mu"]
   gamma=parameters["gamma"]
   N=parameters["N"]
   #Define equations
   dS = mu * (N  - S)  - beta * S * I / N
   dI = beta * S * I / N - (mu + gamma) * I
   dR = gamma * I - mu * R
   res=c(dS, dI, dR)
   #Return list of gradients
   list(res)
 }

#' Gradient-function for the chain-SIR model
#' @param t Implicit argument for time
#' @param logx  A vector with values for the log-states
#' @param parameters A vector with parameter values for the chain-SIR system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 10, by=1/52)
#' paras  = c(mu = 1/75, N = 1, beta =  625, gamma = 365/14, u=5)
#' xstart2 = log(c(S=.06, I=c(0.001, rep(0.0001, paras["u"]-1)), R = 0.0001))
#' out = as.data.frame(ode(xstart2, times, chainSir, paras))
#' @export
chainSir=function(t, logx, parameters) {
    x = exp(logx)
    u = parameters["u"]
    S = x[1]
    I = x[2:(u + 1)]
    R = x[u + 2]
    
    with(as.list(parameters), {
        dS = mu * (N - S) - sum(beta * S * I)/N
        dI = rep(0, u)
        dI[1] = sum(beta * S * I)/N - (mu + u * gamma) * I[1]
        if (u > 1) {
            for (i in 2:u) {
                dI[i] = u * gamma * I[i - 1] - (mu + u * gamma) * 
                  I[i]
            }
        }
        dR = u * gamma * I[u] - mu * R
        res = c(dS/S, dI/I, dR/R)
        list(res)
    })
}


#' Gradient-function for the SEIR model
#' @param t Implicit argument for time
#' @param y  A vector with initial values for the states
#' @param parameters A vector with parameter values for the SEIR system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 10, by=1/120)
#' paras  = c(mu = 1/50, N = 1, beta =  1000, sigma = 365/8, gamma = 365/5)
#' start = c(S=0.06, E=0, I=0.001, R = 0.939)
#' out=ode(start, times, seirmod, paras)
#' @export
seirmod=function(t, y, parameters){
  S=y[1]
  E=y[2]
  I=y[3]
  R=y[4]

  mu=parameters["mu"]
  N=parameters["N"]
  beta=parameters["beta"]
  sigma=parameters["sigma"]
  gamma=parameters["gamma"]

  dS = mu * (N  - S)  - beta * S * I / N
  dE = beta * S * I / N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R 
  res=c(dS, dE, dI, dR)
  list(res)
}

#c4
#' Gradient-function for the forced SEIR model
#' @param t Implicit argument for time
#' @param y  A vector with initial values for the states
#' @param parameters A vector with parameter values for the SIR system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 10, by=1/120)
#' paras  = c(mu = 1/50, N = 1, beta0 = 1000, beta1 = 0.2, sigma = 365/8, gamma = 365/5)
#' start = c(S=0.06, E=0, I=0.001, R = 0.939)
#' out=ode(start, times, seirmod2, paras)
#' @export
seirmod2=function(t, y, parameters){
  S=y[1]
  E=y[2]
  I=y[3]
  R=y[4]

 with(as.list(parameters),{
  dS = mu * (N  - S)  - beta0 * (1+beta1*cos(2*pi*t))* S * I / N
  dE = beta0 * (1+beta1*cos(2*pi * t))* S * I / N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R
  res=c(dS, dE, dI, dR)
  list(res)
})
} 

#' Gradient-function for the age-structured SIR model with possibly heterogeneous mixing
#' @param t Implicit argument for time
#' @param logx  A vector with initial values for the log-states
#' @param parameters A named list with parameter values for the age-structured SIR system. N is population size, gamma is recovery rate, mu is birth/death rate, beta is transmission rate, W is the normalized contact matrix, v is vector of age-class specific vaccination rates and r is class-specific aging rates (since age brackets may differ in width).
#' @return A list of gradients 
#' @examples
#' ra=rep(1,4)
#' n=length(ra)
#' W=matrix(1, ncol=4, nrow=4)
#' paras =list(N=1, gamma=365/14, mu=0.02, beta=500, W=W,v=rep(0,4), r=ra)
#' xstart=log(c(S=rep(0.099/n,n), I=rep(0.001/n,n), R=rep(0.9/n,n)))
#' times=seq(0,10,by=14/365)
#' out=as.data.frame(ode(xstart, times, siragemod, paras))
#' @export
siragemod = function(t, logx,  parameters){
	n=length(parameters$r)
	xx = exp(logx)
        S = xx[1:n]
        I = xx[(n+1):(2*n)]
        R = xx[(2*n+1):(3*n)]
        
	with(as.list(parameters), {
                phi = (beta*W%*%I)/N
                dS = c(mu,rep(0,n-1))*N - (phi+r)*S + c(0,r[1:(n-1)]*S[1:(n-1)]) - mu*S - v*S
                dI = phi*S + c(0,r[1:(n-1)]*I[1:(n-1)]) -(gamma+r)*I - mu*I
                dR =  v*S + c(0,r[1:(n-1)]*R[1:(n-1)]) + gamma*I - r*R - mu*R
                res = c(dS/S,dI/I,dR/R)
		list((res))
	})
}

