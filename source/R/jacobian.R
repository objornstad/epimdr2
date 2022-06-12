#' A Jacobian matrix calculator
#'
#' A general-purpose function to construct and evaluate Jacobian matrices
#'
#' @param states a vector naming all state variables
#' @param elist a list that contains equations (as quotes) for all state variables
#' @param parameters a labeled vector of parameters
#' @param pts a labeled vector of the point in the phase plane in which to evaluate the Jacobian 
#' (often the endemic or disease-free equilibrium if working in mathematical epidemiology)
#' @return The Jacobian matrix
#' @examples
#' #The SEIR model
#' states=c("S", "E", "I", "R")
#' elist=c(dS = quote(mu * (N  - S)  - beta * S * I / N),
#'  dE = quote(beta * S * I / N - (mu + sigma) * E),
#' dI = quote(sigma * E - (mu + gamma+alpha) * I),
#' dR = quote(gamma * I - mu * R))
#' paras  = c(mu = 1/50, N = 1, beta =  1000, 
#' sigma = 365/8, gamma = 365/5, alpha=0)
#' deq=list(S = 1, E = 0, I = 0, R = 0)
#' jacobian(states=states, elist=elist, parameters=paras, pts=deq)
#' @export
jacobian=function(states, elist, parameters, pts){
paras = as.list(c(pts, parameters)) 

k=0
jl=list(NULL)
for(i in 1:length(states)){
assign(paste("jj", i, sep = "."), lapply(lapply(elist, deriv, states[i]), eval, paras))
for(j in 1:length(states)){
k=k+1
jl[[k]]=attr(eval(as.name(paste("jj", i, sep=".")))[[j]], "gradient")[1,]
}
}

J=matrix(as.numeric(as.matrix(jl)[,1]), ncol=length(states))
return(J)
}
