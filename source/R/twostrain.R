#' Gradient-function for the twostrain model
#' @param t Implicit argument for time
#' @param y  A vector with initial values for the states
#' @param parameters A vector with parameter values for the system
#' @return A list of gradients 
#' @examples
#' require(deSolve)
#' times  = seq(0, 5, by=1/100)
#' paras  = c(mu = 5, N = 1, R01=2, R02=1.5, gamma = 365/14, Theta=0.9, Xi=0.9, Pi=0.9)
#' start = c(S=0.998, I1=0.001, I2=0.001,R1 = 0, R2=0, J1=0, J2=0, R=0)
#' out=as.data.frame(ode(start, times, twostrain, paras))
#' @export
 
twostrain=function(t, y, parameters){
  S=y[1]
  I1=y[2]
  I2=y[3]
  R1=y[4]
  R2=y[5]
  J1=y[6]
  J2=y[7]
  R=y[8]
beta1 = with(as.list(parameters), R01*(365*gamma+mu))
beta2 = with(as.list(parameters), R02*(365*gamma+mu))

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
})}
