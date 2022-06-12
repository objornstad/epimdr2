
#' Launch a shiny-app simulating the seasonal SEIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{seir.app}
#' @export
seir.app=shinyApp(
# This creates the User Interface (UI)
ui = pageWithSidebar(
headerPanel("Seasonally forced SEIR"),
sidebarPanel(
sliderInput("R0", "R0:", 15,
              min = 0.5, max = 20),
sliderInput("beta1", "Seasonality (\\%):", 10,
              min = 0, max = 100),
sliderInput("Ip", "Infectious period (days)", 5,
              min = 1, max = 100),
sliderInput("oneoversigma", "Latent period (days):", 8,
              min = 1, max = 100),
sliderInput("mu", "birth rate (yr^-1):", 0.02,
              min = 0, max = .1),
sliderInput("T", "Time range:",
                  min = 0, max = 100, value = c(50,70)),
checkboxInput("lg", "un-Log", TRUE)
),
mainPanel(
  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("Phase plane", plotOutput("plot2")),
       tabPanel("Details", 
           withMathJax(
       helpText("MODEL:"),
            helpText("Susceptible $$\\frac{dS}{dt} = \\underbrace{\\mu N}_{birth} - \\underbrace{\\mu S}_{death} - \\underbrace{\\frac{\\beta(t) I S}{N}}_{infection}$$"),
            helpText("Exposed $$\\frac{dE}{dt} = \\underbrace{\\frac{\\beta(t) I S}{N}}_{infection} - \\underbrace{\\mu E}_{death}-\\underbrace{\\sigma E}_{latency}$$"),
            helpText("Infectious $$\\frac{dI}{dt} = \\underbrace{\\sigma E}_{latency} -\\underbrace{(\\mu I}_{death} - \\underbrace{\\gamma I}_{recovery}$$"),
           helpText("Removed $$\\frac{dR}{dt} = \\underbrace{\\gamma I}_{recovery}  - \\underbrace{(\\mu R}_{death}$$"),
           helpText("Mean transmission $$\\beta_0 = R_0 (\\sigma+\\mu) (\\gamma+\\mu) / \\sigma$$"),
           helpText("Seasonality $$\\beta(t) =  \\beta_0 (1 + \\beta_1 cos(2 \\pi t))$$"),
           helpText("Reproduction number $$R_0 =  \\frac{\\sigma}{\\sigma +\\mu} \\frac{1}{\\gamma+\\mu} \\frac{\\beta N}{N} = \\frac{\\sigma}{\\sigma +\\mu} \\frac{\\beta}{\\gamma+\\mu}$$"),             
            helpText("REFERENCE: Earn, D.J.D., Rohani, P., Bolker, B.M. and Grenfell, B.T. (2000) A simple model for complex dynamical transitions in epidemics.
             Science 287: 667-670"))
           ))
  
  )
),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
  seirmod2=function(t, x, params){
  S=x[1]
  E=x[2]
  I=x[3]
  R=x[4]

  mu=params["mu"]
  N=params["N"]
  beta0=params["R0"]*(params["sigma"]+params["mu"])*(params["gamma"]+params["mu"])/params["sigma"]
  beta1=params["beta1"]/100
  sigma=params["sigma"]
  gamma=params["gamma"]

  dS = mu * (N  - S)  - beta0 * (1+beta1*cos(2*pi*t))* S * I / N
  dE = beta0 * (1+beta1*cos(2*pi * t))* S * I / N - (mu + sigma) * E
  dI = sigma * E - (mu + gamma) * I
  dR = gamma * I - mu * R
  res=c(dS, dE, dI, dR)
  list(res)
} 



  output$plot1 <- renderPlot({

  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = input$mu, N = 1, R0 = input$R0, beta1 = input$beta1, sigma = 365/input$oneoversigma, gamma = 365/input$Ip)
  xstart = c(S=0.06, E=0, I=0.001, R = 0.939)
  R0 = round(input$R0, 1)
 
out=ode(y=xstart,
  times=times,
  func=seirmod2,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

par(mar = c(5,5,2,5))
#lg=ifelse(input$lg==TRUE, "y", "")
plot(x=out$time[sel], y=out$I[sel], ylab="fraction", xlab="time", type="l",
ylim=range(out[sel,-c(1,2, 5)]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""), col="red")
 lines(x=out$time, y=out$E, col="blue")
title(paste("R0=", R0))
# lines(x=out$time, y=out$S, col="green")
par(new=T)
plot(x=out$time, y=out$S, type="l", col="green", axes=FALSE, xlab=NA, ylab=NA, 
    ylim=range(out[sel,2]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""))
axis(side = 4, col="green")
mtext(side = 4, line = 4, "S", col="green")
  legend("right",
        legend=c("I", "E", "S"),
        lty=c(1,1,1),
         col=c("red", "blue", "green"))
   })
  
output$plot2 <- renderPlot({
  times  = seq(0, input$T[2], by=1/100)
  paras  = c(mu = input$mu, N = 1, R0 = input$R0, beta1 = input$beta1, sigma = 365/input$oneoversigma, gamma = 365/input$Ip)
  xstart = c(S=0.06, E=0, I=0.001, R = 0.939)
  R0 = round(input$R0, 1)
 
  out=ode(y=xstart,
  times=times,
  func=seirmod2,
  parms=paras)

  out=as.data.frame(out)

  sel=out$time>input$T[1]&out$time<input$T[2]

  plot(out$S[sel], out$I[sel], log=ifelse(input$lg==TRUE, "xy", ""), type="l", xlab="fraction susceptible", ylab="fraction infected")
  abline(v=1/R0, col="green")
  curve(paras["mu"]*(1-x)/(paras["beta0"]*x), min(out$S), max(out$S), add=TRUE, col="red")
    legend("topright",
        legend=c("S-isocline", "I-isocline"),
        lty=c(1,1),
         col=c("red", "green"))
 
  })

  }
)
