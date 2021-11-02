

#' Launch a shiny-app simulating the SIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{sir.app}
#' @export
sir.app=shinyApp(
# This creates the User Interface (UI)
ui <- pageWithSidebar(
headerPanel("The SIR model"),
#The sidebar for parameter input
sidebarPanel(
#Sliders:
sliderInput("R0", "R0:", 2,
              min = 0.5, max = 20),
sliderInput("infper", "Infectious period (days)", 5,
              min = 1, max = 100),
sliderInput("mu", "birth rate (yr^-1):", 5,
              min = 0, max = 100),
sliderInput("T", "Time range:",
                  min = 0, max = 1, value = c(0,1))
),
#Main panel for figures and equations
mainPanel(
  #Multiple tabs in main panel
  tabsetPanel(
      #Tab 1: Time plot (plot1 from server)
      tabPanel("Time", plotOutput("plot1")), 
      #Tab 2: Phase plot (plot2 from server)
      tabPanel("Phase plane", plotOutput("plot2", height = 500)),
      #Tab 3: MathJax typeset equations 
      tabPanel("Equations", 
           withMathJax(
             helpText("Susceptible $$\\frac{dS}{dt} = \\underbrace{\\mu N}_{birth} - \\underbrace{\\mu S}_{death} - \\underbrace{\\frac{\\beta(t) I S}{N}}_{infection}$$"),
               helpText("Infe ctious $$\\frac{dI}{dt} = \\underbrace{\\frac{\\beta(t) I S}{N}}_{infection} - \\underbrace{\\mu I}_{death}-\\underbrace{\\gamma I}_{recovery}$$"),
            helpText("Removed $$\\frac{dR}{dt} = \\underbrace{\\gamma I}_{recovery}  - \\underbrace{(\\mu R}_{death}$$"),
           helpText("Transmission rate $$\\beta = R_0  (\\gamma + \\mu)$$"),
           helpText("Reproduction number $$R_0 =  \\frac{1}{\\gamma + \\mu} \\frac{\\beta N}{N} = \\frac{\\beta}{\\gamma + \\mu}$$")             
           ))
  ))), #End of ui()


# This creates the 'behind the scenes' code (Server)
server <- function(input, output) {
  #Gradient function for SIR model
  sirmod=function(t, x, parameters){
    S=x[1]
    I=x[2]
    R=x[3]
    R0=parameters["R0"]
    mu=parameters["mu"]
    gamma=parameters["gamma"]
    N=parameters["N"]
    beta=R0*(gamma+mu)
    dS = mu * (N  - S)  - beta * S * I / N
    dI = beta * S * I / N - (mu + gamma) * I
    dR = gamma * I - mu * R
    res=c(dS, dI, dR)
    list(res)
  }

 #Plot1: renderPlot to be passed to UI tab 1
  output$plot1 <- renderPlot({
  #input\$xx's are pulled from UI
  times  = seq(0, input$T[2], by=1/1000)
  paras  = c(mu = input$mu, N = 1, R0 =  input$R0, gamma =
    365/input$infper)
  start = c(S=0.999, I=0.001, R = 0)
  paras["beta"] = with(as.list(paras), R0*(gamma+mu))
  #Resonant period
  AA=with(as.list(paras), 1/(mu*(R0-1)))
  GG=with(as.list(paras), 1/(mu+gamma))
  rp=round(2*pi*sqrt(AA*GG),2)

  #Integrate ode with parameters pulled from UI
  out=ode(start,  times, sirmod, paras)
  out=as.data.frame(out)

  #Plot1
  sel=out$time>input$T[1]&out$time<input$T[2]
  plot(x=out$time[sel], y=out$S[sel], ylab="fraction", xlab="time", type="l",
  ylim=range(out[sel,-c(1,4)]))
  title(paste("R0=", paras["R0"], "Period=", rp))
  lines(x=out$time[sel], y=out$I[sel], col="red")
  lines(x=out$time[sel], y=out$R[sel], col="green")
  legend("right",
        legend=c("S", "I", "R"),
        lty=c(1,1,1),
         col=c("black", "red", "green"))
   })

 #Plot2: renderPlot to be passed to UI tab 2
  output$plot2 <- renderPlot({
  times  = seq(0, input$T[2], by=1/1000)
  paras  = c(mu = input$mu, N = 1, R0 =  input$R0, gamma =
    365/input$infper)
  paras["beta"] = with(as.list(paras), R0*(gamma+mu))

  start = c(S=0.999, I=0.001, R = 0)
 
  #Gradient function used for phaseR phase-plot
  simod=function(t, y, parameters){
   S=y[1]
   I=y[2]
   beta=parameters["beta"]
   mu=parameters["mu"]
   gamma=parameters["gamma"]
   N=parameters["N"]   
   dS = mu * (N  - S)  - beta * S * I / N
   dI = beta * S * I / N - (mu + gamma) * I
   res=c(dS, dI)
   list(res)
  }

  #Integrate simod
  out=ode(start[-3], times, simod, paras)
  out=as.data.frame(out)

  AA=with(as.list(paras), 1/(mu*(R0-1)))
  GG=with(as.list(paras), 1/(mu+gamma))
  rp=round(2*pi*sqrt(AA*GG),2)
  
  plot(x=out$S, y=out$I, xlab="Fraction suceptible", ylab="Fraction infected", type="l")
   title(paste("R0=", paras["R0"], "Period=", rp))
 #Add vector field
  fld=flowField(simod, xlim=range(out$S), ylim=range(out$I), 
  parameters=paras, system="two.dim", add=TRUE,
  ylab="I", xlab="S")
  #Add isoclines
  abline(v=1/paras["R0"], col="green")
  curve(paras["mu"]*(1-x)/(paras["beta"]*x), min(out$S), max(out$S), add=TRUE, col="red")
    legend("topright",
        legend=c("S-isocline", "I-isocline"),
        lty=c(1,1),
         col=c("red", "green"))
   })
  } #End of server()
)
