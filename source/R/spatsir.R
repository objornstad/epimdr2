#' A function to caculate the matrix of gravity coupling based on distance and population size
#' @param tau1 recipient exponent
#' @param tau2 donor exponent
#' @param rho distance exponent
#' @param pop a vector of population sizes
#' @param distance a matrix of distances
#' @return A matrix of gravity coupling 
#' @examples
#' require(ncf)
#' data(usflu)
#' usdist = gcdist(usflu$Longitude, usflu$Latitude)
#' G = gravity(0.3, 0.6, 3, usflu$Pop, usdist)
#' @seealso
#' \code{\link{sirSpatmod}}
#' @export
gravity = function(tau1, tau2, rho, pop, distance){
   gravity = outer(pop^tau1, pop^tau2)/distance^rho
   diag(gravity) = 0
   gravity}


#' Gradient function for a spatially-extended SIR model
#'
#' Gradient function for a spatially-extended SIR model given some spatial topology
#' @param t Implicit argument for time
#' @param y  A vector of length L*3 with initial values for the states. The first 1:L represents intial S's, (L+1):2*L are initial I's and the last (2*L+1):3*L are initial R's
#' @param parameters A vector with parameter values for the spatial SIR system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' require(ncf)
#' data(usflu)
#' usdist = gcdist(usflu$Longitude, usflu$Latitude)
#' G = gravity(0.3, 0.6, 3, usflu$Pop, usdist)
#' gamma = 1/3.5
#' R0 = 1.8
#' beta = R0 * gamma/usflu$Pop
#' m = 1 / 1000 / sum(usflu$Pop)
#' parms = list(beta = beta, m = m, gamma  =  gamma, G = G) 
#' S = usflu$Pop
#' R = I = rep(0, length(usflu$Pop))
#' I[31] = 1 
#' inits = c(S = S, I = I, R = R)  
#' times = 0:200
#' out = ode(inits, times, sirSpatmod, parms)
#' L=length(usflu$Pop)
#' \dontrun{matplot(out[, 50+(1:L)], type = "l", ylab = "Prevalence", xlab = "Day")}
#' @export
 sirSpatmod = function(t, y, parameters){
  L=length(y)/3
  i = c(1:L)
  S = y[i]
  I = y[L+i]
  R = y[2*L+i]
  with(parameters,{
  beta = beta[i]
  dS = -(beta*I + m*G%*%I)*S
  dI = (beta*I + m*G%*%I)*S - gamma*I
  dR = gamma*I 
  list(c(dS, dI, dR)) 
})
}

#