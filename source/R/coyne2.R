#c9
#' Edition 2 Gradient-function for Coyne et al's rabies model
#' @param t Implicit argument for time
#' @param logx A vector with values for the log-states
#' @param parms A vector with parameter values for the dynamical system
#' @return A list of gradients for the log system
#' @examples
#' require(deSolve)
#' times  = seq(0, 50, by=1/520)
#' paras  = c(gamma = 0.0397, b = 0.836, a = 1.34, sigma = 7.5, 
#' alpha = 66.36, beta = 33.25, c = 0, rho = 0.8)
#' start = log(c(S=12.69/2, E1=0.1, E2=0.1, I = 0.1, R = 0.1))
#' out = as.data.frame(ode(start, times, coyne, paras))
#' @export
coyne2 = function(t, logx, parms){
  x = exp(logx)
  S = x[1]
  E1 = x[2]
  E2 = x[3]
  I = x[4]
  R = x[5]
  N = sum(x)
  with(as.list(parms),{
  dS = a * (S + R) - beta * S * I - 
       d * N * S  - (b + c) * S
  dE1= lambda * beta * S * I  - d * N * E1  -
       (b + sigma + c) * E1
  dE2= (1-lambda) * beta * S * I  - d * N * E2  -
       (b + sigma + c) * E2
  dI = sigma * E1  - d * N * I  -
       (b + alpha + c) * I
  dR = sigma * E2 - d * N * R  - (b + c) * R
  res = c(dS/S, dE1/E1, dE2/E2, dI/I, dR/R)
  list(res)
})
}



#' Edition 1 Gradient-function for Coyne et al's rabies model
#' @param t Implicit argument for time
#' @param logx A vector with values for the log-states
#' @param parms A vector with parameter values for the dynamical system
#' @return A list of gradients for the log system
#' @examples
#' require(deSolve)
#' times  = seq(0, 50, by=1/520)
#' paras  = c(gamma = 0.0397, b = 0.836, a = 1.34, sigma = 7.5, 
#' alpha = 66.36, beta = 33.25, c = 0, rho = 0.8)
#' start = log(c(X=12.69/2, H1=0.1, H2=0.1, Y = 0.1, I = 0.1))
#' out = as.data.frame(ode(start, times, coyne, paras))
#' @export
coyne=function(t, logx, parms){
  x=exp(logx)
  X=x[1]
  H1=x[2]
  H2=x[3]
  Y=x[4]
  I=x[5]
  N = sum(x)
  
  with(as.list(parms),{
  dX = a * (X + I) - beta * X * Y - gamma * N * X  - (b + c) * X
  dH1= rho * beta * X * Y  - gamma * N * H1  - (b + sigma + c) * H1
  dH2= (1-rho) * beta * X * Y  - gamma * N * H2  - (b + sigma + c) * H2
  dY = sigma * H1  - gamma * N * Y  - (b + alpha + c) * Y
  dI = sigma * H2 - gamma * N * I  - (b + c) * I
  res=c(dX/X, dH1/H1, dH2/H2, dY/Y, dI/I)
  list(res)
})
}

