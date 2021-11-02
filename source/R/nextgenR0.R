#' Next generation matrix R0 calculator
#'
#' Calculates R0 for arbitrarily complex compartmental flows using the method of Diekmann et al. (1990).
#'
#' @param Istates a vector naming all Infected classes
#' @param Flist a list that contains equations (as quotes) for completely new infections entering each infected 
#' compartment for each class
#' @param Vlist a list that contains the equations (as quotes) for losses out of each infected compartment minus the 
#' equations (as quotes) for all gains into each infected compartment that does not
#' represent new infections but transfers among infectious classes
#' @param parameters a labeled vector of parameters
#' @param dfe a labeled vector of all states at the disease-free equilibrium
#' @return The next generation matrix estimate of R0 
#' @source Diekmann, O., Heesterbeek, J. A. P. and Metz, J. A. J. 1990.
#' On the Definition and the Computation of the Basic Reproduction Ratio R0 in Models for Infectious-Diseases in 
#' Heterogeneous Populations. Journal of Mathematical Biology 28: 365-382.
#' @examples
#' #The SEIR model
#' #Infected classes are $E$ and $I$
#' istates=c("E", "I")
#' flist=c(dEdt=quote(beta * S * I / N), dIdt=quote(0))
#' Vm1=quote(mu * E + sigma * E)
#' Vm2=quote(mu * I + alpha * I + gamma * I)
#' Vp1=0
#' Vp2=quote(sigma * E)
#' V1=substitute(a-b, list(a=Vm1, b=Vp1))
#' V2=substitute(a-b, list(a=Vm2, b=Vp2))
#' vlist = c(V1,V2)
#' para = list(mu = 0, alpha = 0, beta = 5, gamma = .8, sigma = 1.2, N = 1)
#' df = list(S = 1, E = 0, I = 0, R = 0)
#' nextgenR0(Istates=istates, Flist=flist, Vlist=vlist, parameters=para, dfe=df)
#' @export
#' @importFrom stats deriv
nextgenR0=function(Istates, Flist, Vlist, parameters, dfe){
paras = as.list(c(dfe, parameters)) 

k=0
vl=fl=list(NULL)
for(i in 1:length(Istates)){
assign(paste("f", i, sep = "."), lapply(lapply(Flist, deriv, Istates[i]), eval, paras))
assign(paste("v", i, sep = "."), lapply(lapply(Vlist, deriv, Istates[i]), eval, paras))
for(j in 1:length(Istates)){
k=k+1
fl[[k]]=attr(eval(as.name(paste("f", i, sep=".")))[[j]], "gradient")[1,]
vl[[k]]=attr(eval(as.name(paste("v", i, sep=".")))[[j]], "gradient")[1,]
}
}

f=matrix(as.numeric(as.matrix(fl)[,1]), ncol=length(Istates))
v=matrix(as.numeric(as.matrix(vl)[,1]), ncol=length(Istates))
R0=max(eigen(f%*%solve(v))$values)
return(R0)
}
