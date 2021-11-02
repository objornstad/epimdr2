globalVariables("x")

#' Launch a shiny-app simulating May's Parasitoid-host Model model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{may.app}
#' @export
#' @importFrom graphics abline curve legend lines plot title
#' @importFrom shiny renderPlot
#' @importFrom deSolve ode
#' @importFrom polspline lspec
may.app=shinyApp(
# This creates the User Interface (UI)
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("May's Parasitoid-host Model"),
sidebarPanel(
sliderInput("R", "Growth rate (R):", 1.1,
              min = 1, max = 2, step=.01),
sliderInput("a", "Search efficiency (a):", 0.1,
              min = 0, max = .5),
sliderInput("k", "aggregation (k):", 1.5,
              min = 0.1, max = 3, step=0.1),
numericInput("P0", "Initial parasitoid:", 10,
              min = 1, max = 100),
numericInput("H0", "Initial host:", 20,
              min = 1, max = 100),
numericInput("Tmax", "Tmax:", 100,
              min = 1, max = 500)
),
mainPanel(tabsetPanel(
  tabPanel("Simulation", plotOutput("plot1", height = 500)),
  tabPanel("Phase plane", plotOutput("plot2", height = 500)),
  tabPanel("Details", 
    withMathJax(
         helpText("MODEL:"),
             helpText("Host $$H_t = R H_{t-1} (1 + a P_{t-1}/k)^{-k}$$"),
          helpText("Parasitoid $$P_t = R H_{t-1} (1-(1 + a P_{t-1}/k)^{-k})$$"),
          helpText("REFERENCE: May RM (1978) Host-parasitoid systems in patchy 
            environments: a phenomenological model. J Anim Ecol 47: 833-843"))
)
)
)
),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
 nbmod = function(R, a, k, T = 100, H0 = 10, P0 = 1){
   #T is length of simulation (number of time-steps)
   #H0 and P0 are initial numbers
   #we provide default parameters except for R and a

   H=rep(NA, T) #host series
   P=rep(NA, T) #parasitoid series

   H[1] = H0 #Initiating the host series
   P[1] = P0 #Initiating the host series

   for(t in 2:T){
     H[t] = R * H[t-1] * (1+ a * P[t-1]/k)^(-k)
     P[t] = R * H[t-1] * (1-(1+ a * P[t-1]/k)^(-k))
     if(P[t-1]==0) break
   } #end of loop

   #the two vectors of results are stored in a "list"
   res= data.frame(H = H, P = P)
 
   #the list is passed out of this function
   return(res)
} #end of function



  output$plot1 <- renderPlot({

    sim= nbmod(R=input$R, a=input$a, k=input$k, H0=input$H0, P0=input$P0, T=input$Tmax)
    time = 1:input$Tmax

    plot(time, sim$H, type= "b",xlab = "Generations", ylab = "Abundance", 
      ylim = range(sim, na.rm=TRUE))
    points(time, sim$P, type = "b", col = "red")
     legend("topleft",
        legend=c("H", "P"),
        lty=c(1,1),
        pch=c(1,1),
        col=c("black", "red"))
   })

output$plot2 <- renderPlot({

    sim= nbmod(R=input$R, a=input$a, k=input$k, H0=input$H0, P0=input$P0, T=input$Tmax)
    time = 1:input$Tmax
   
    plot(sim$H, sim$P, type= "b",xlab = "Host", ylab = "Parasitoid")
    #Pstar=input$k*(input$R^(1/input$k)-1)/input$a
    #Hstar=Pstar*input$R/(input$R-1)
    #points(Hstar, Pstar, col=2, pch=19)
   })
  }
)



#' Launch a shiny-app simulating TSIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{tsir.app}
#' @export
tsir.app=shinyApp(
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("Simulating with TSIR"),
sidebarPanel(
sliderInput("alpha", "alpha:", 0.97,
              min = 0.8, max = 1),
sliderInput("beta", "beta:", 25,
              min = 0, max = 100),
sliderInput("B", "Births (B):", 2300,
              min = 0, max = 5000),
sliderInput("S0", "Fraction S:", 0.06,
              min = 0, max = 1),
numericInput("sdbeta", "sdbeta:", 3,
              min = 0, max = 10),
numericInput("N", "Popsize:", 3000000,
              min = 0, max = 100),
numericInput("I0", "initial I:", 100,
              min = 1, max = 100),
numericInput("IT", "Iterations:", 520,
              min = 1, max = 1000)
),
mainPanel(
  tabsetPanel(
      tabPanel("Simulation", plotOutput("plot1")), 
      tabPanel("Transfer function", plotOutput("plot2")),
      tabPanel("Details", 
           withMathJax(
            helpText("MODEL:"),
            helpText("Susceptible $$S_{t+1} = S_t - I_{t+1} + B$$"),
            helpText("Expected Inf $$\\lambda_{t+1} = \\frac{\\beta_t I_t^\\alpha S_t}{N}$$"),
           helpText("Infected $$I_{t+1} \\sim \\mbox{Poisson}(\\lambda_{t+1})$$"),
           helpText("Transmission $$\\beta_t \\sim \\mbox{Norm}(\\beta, \\mbox{sd}\\beta^2)$$"),             
            helpText("The transfer function is $$T(\\omega) = (\\vec{I}-e^{-\\imath \\omega} \\vec{J})^{-1} \\cdot \\vec{A}$$")
          )))   
  )),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
  simTsir=function(alpha, B, beta, sdbeta,
    S0, I0, IT, N){
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

  output$plot1 <- renderPlot({
  out = simTsir(alpha=input$alpha, B=input$B, beta=input$beta, sdbeta=input$sdbeta, 
  S0=input$S0, I0=input$I0, IT=input$IT, N=input$N)
par(mfrow=c(1,2))  #This puts two plots side by side each other
plot(out$I, ylab="infected", xlab="time", type="b")
plot(out$S, out$I, ylab="infected", xlab="susceptible", type="b")
})

  output$plot2 <- renderPlot({
  out = simTsir(alpha=input$alpha, B=input$B, beta=input$beta, sdbeta=input$sdbeta, 
  S0=input$S0, I0=input$I0, IT=input$IT, N=input$N)
Seq=expression(S-beta*S*I^alpha/N+B)
Ieq=expression(beta*S*I^alpha/N)
j11=D(Seq, "S")
j12=D(Seq, "I")
j21=D(Ieq, "S")
j22=D(Ieq, "I")
jj=c(j11, j12, j21,j22)

a1=D(Seq, "beta")
a2=D(Ieq, "beta")
aa=c(a1, a2)

paras  = c(B = input$B, beta =  input$beta, alpha = input$alpha, N=input$N)
eqs=sapply(c(quote(B^(1-alpha)*N/beta), quote(B)), eval, as.list(paras))

J=matrix(sapply(jj, eval, as.list(c(paras, c(S=eqs[1], I=eqs[2])))), 2, byrow=TRUE)
evs=eigen(J)$values
rp=2*pi/atan2(Im(evs[1]), Re(evs[1]))

A=matrix(sapply(aa, eval, as.list(c(paras, c(S=eqs[1], I=eqs[2])))), 2, byrow=TRUE)
Id=matrix(c(1,0,0,1),ncol=2)
wseq=seq(0,pi,length=500)

Fr=vector("list",500)  #set up empty list of matrices
# loop to fill those matrices with fourier trasform for the 500 values of w
for(i in 1:500){ 
Fr[[i]]=matrix(solve(Id-exp(1i*wseq[i])*J)%*%A,ncol=1)  #solve gives inverse
}

PS=matrix(NA,ncol=2,nrow=500,dimnames=list(1:500, c("S","I"))) 
#power spectra from real and imaginary parts of Fourier transform
for(i in 1:500){
PS[i,]=sqrt(Re(Fr[[i]])^2+Im(Fr[[i]])^2)
}

sfit=spectrum(out$I[-c(1:104)])

sfit2=lspec(out$I[-c(1:104)])
plot(wseq, PS[,2], type="l", 
  xlab="frequency (in radians)", ylab="amplitude", xlim=c(0,0.6))
title("Simulated spectrum (periodogram and log-spline) \n and T-fn prediction")
lines(pi*sfit$freq/0.5, max(PS[,2])*sfit$spec/max(sfit$spec), col=2)
par(new=TRUE)
plot(sfit2, col=3, xlim=c(0,0.6), axes=FALSE)
legend("topright", c("T-function", "periodogram", "log-spline"), lty=c(1,1,1), col=c(1,2,3))
   })
  }
)

