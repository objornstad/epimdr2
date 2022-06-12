#' Gradient-function for the two-strain SIR model
#' @param t Implicit argument for time
#' @param y  A vector with initial values for the states
#' @param parameters A vector with parameter values for the two-strain SIR  system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 30, by=1/200)
#' paras  = c(mu = 0.02, N = 1, beta1=500, beta2=750, gamma = 365/5, Theta=0.15, Xi=0.15, Pi=0.8)
#' start = c(S = 0.999, I1 = 0.001, I2 = 0.00, R1=0, R2=0, J1=0, J2=0, R = 0)
#' out = as.data.frame(ode(start, times, twostrain, paras))
#' @export
twostrain=function(t, y, parameters){
  S=ifelse(y[1]<0, 0, y[1])
  I1=ifelse(y[2]<0, 0, y[2])
  I2=ifelse(y[3]<0, 0, y[3])
  R1=ifelse(y[4]<0, 0, y[4])
  R2=ifelse(y[5]<0, 0, y[5])
  J1=ifelse(y[6]<0, 0, y[6])
  J2=ifelse(y[7]<0, 0, y[7])
  R=ifelse(y[8]<0, 0, y[8])

 with(as.list(parameters),{
  phi=(beta1*I1+beta2*I2+Theta*(beta1*J1+beta2*J2))/N
  dS = mu*N - phi *S -mu * S
  dI1= (beta1*I1+Theta*beta1*J1)*S/N - (gamma+mu)*I1
  dI2= (beta2*I2+Theta*beta2*J2)*S/N - (gamma+mu)*I2
  dR1= Pi*gamma*I1 - (beta2*I2+Theta*beta2*J2)*Xi*R1/N - mu*R1
  dR2= Pi*gamma*I2 - (beta1*I1+Theta*beta1*J1)*Xi*R2/N - mu*R2
  dJ1= (beta1*I1+Theta*beta1*J1)*Xi*R2/N - gamma*J1 - mu*J1
  dJ2= (beta2*I2+Theta*beta2*J2)*Xi*R1/N - gamma*J2 - mu*J2
  dR=(1-Pi)*gamma*(I1+I2)+gamma*(J1+J2)-mu*R

  res=c(dS, dI1, dI2, dR1,dR2,dJ1,dJ2,dR)
  list(res)
})
} 
