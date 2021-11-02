#' Launch a shiny-app simulating a two-strain SIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{twostrain.app}
#' @export
# This creates the User Interface (UI)
twostrain.app=shinyApp(
ui = pageWithSidebar(
headerPanel("A two-strain SIR model"),
sidebarPanel(
sliderInput("R01", "R01:", 2,
              min = 0.5, max = 5),
sliderInput("R02", "R02:", 1.5,
              min = 0.5, max = 5),
sliderInput("Ip", "Infectious period (days)", 5,
              min = 1, max = 100),
sliderInput("mu", "birth rate (per year):", 0.02,
              min = 0, max = .1),
sliderInput("theta", "Theta:", 0.9,
              min = 0, max = 1),
sliderInput("xi", "Xi:", 0.9,
              min = 0, max = 1),
sliderInput("pi", "Pi:", 0.9,
              min = 0, max = 1),
sliderInput("T", "Time range:",
                  min = 0, max = 10, value = c(0,5)),
checkboxInput("lg", "un-Log", TRUE),
width=3
),
mainPanel(
  tabsetPanel(
      tabPanel("Time", plotOutput("plot1")), 
      tabPanel("Phase plane", plotOutput("plot2")),
       tabPanel("Details", 
           withMathJax(
       helpText("MODEL:"),
            helpText("Susceptible $$\\frac{dS}{dt} = \\mu (N - S) - \\phi$$"),
            helpText("$$\\phi = \\frac{\\beta_1 * I1 +\\beta_2 * I2+\\Theta * (\\beta_1 * J1+\\beta_2* J_2)}{N}$$"),
            helpText("Strain 1 I $$\\frac{dI_1}{dt} = \\frac{(\\beta_1 * I_1 + \\Theta * \\beta_1 * J_1) * S}{N} - (\\gamma + \\mu) * I_1$$"),
            helpText("Strain 2 I $$\\frac{dI_2}{dt} = \\frac{(\\beta_2 * I_2 + \\Theta * \\beta_2 * J_2) * S}{N} - (\\gamma + \\mu) * I_2$$"),
            helpText("Immune to 1 only $$\\frac{dR_1}{dt} = \\Pi * \\gamma * I_1 - (\\beta_2 * I_2 + \\Theta * \\beta_2 * J_2) * \\Xi * R_1 / N -  \\mu * R_1$$"),
            helpText("Immune to 2 only $$\\frac{dR_2}{dt} = \\Pi * \\gamma * I_2 - (\\beta_1 * I_1 + \\Theta * \\beta_1 * J_1) * \\Xi * R_2 / N -  \\mu * R_2$$"),
            helpText("Strain 1 J $$\\frac{dJ_1}{dt} = (\\beta_1 * I_1 + \\Theta * \\beta_1 * J_1) * \\Xi * R_2 / N - \\gamma * J_1 - \\mu * J_1$$"),
            helpText("Strain 2 J $$\\frac{dJ_2}{dt} = (\\beta_2 * I_2 + \\Theta * \\beta_2 * J_2) * \\Xi * R_1 / N - \\gamma * J_2 - \\mu * J_2$$"),
            helpText("Removed $$\\frac{dR}{dt} = (1-\\Pi) * \\gamma * (I_1 + I_2) + \\gamma * (J_1 + J_2) - \\mu * R$$"),  
            helpText("REFERENCE: ")
           ))
  
  )
)
),

# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
twostrain=function(t, y, parameters){
  S=y[1]
  I1=y[2]
  I2=y[3]
  R1=y[4]
  R2=y[5]
  J1=y[6]
  J2=y[7]
  R=y[8]

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



  output$plot1 <- renderPlot({
      times  = seq(0, input$T[2], by=1/100)
paras  = c(mu = input$mu, N = 1, R01=input$R01, R02=input$R02, 
     gamma = 365/input$Ip, Theta=input$theta, Xi=input$xi, Pi=input$pi)
paras["beta1"] = with(as.list(paras), R01*(365*gamma+mu))
paras["beta2"] = with(as.list(paras), R02*(365*gamma+mu))
start = c(S = 0.06, I1 = 0.001, I2 = 0.001, R1=0, R2=0, J1=0, J2=0, R = 0.938)

out = as.data.frame(ode(start, times, twostrain, paras))

sel=out$time>input$T[1]&out$time<input$T[2]

par(mar = c(5,5,2,5))
#lg=ifelse(input$lg==TRUE, "y", "")
plot(x=out$time[sel], y=out$I1[sel], ylab="fraction", xlab="time", type="l",
ylim=range(out[sel, 3:8]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""), col="black")
 lines(x=out$time, y=out$I2, col="red")
 lines(x=out$time, y=out$J1, lty=2)
 lines(x=out$time, y=out$J2, col="red", lty=2)

 par(new=T)
#plot(x=out$time, y=out$S, type="l", col="green", axes=FALSE, xlab=NA, ylab=NA 
#    ylim=range(out[sel,2]), xlim=c(input$T[1], input$T[2]), #log=ifelse(input$lg==TRUE, "y", ""))
#axis(side = 4, col="green")
#mtext(side = 4, line = 4, "S", col="green")
  legend("right",
        legend=c("I1", "I2", "J1", "J2"),
        lty=c(1,1,2,2),
         col=c("black", "red", "black", "red"))
   })
  
output$plot2 <- renderPlot({
      times  = seq(0, input$T[2], by=1/100)
paras  = c(mu = input$mu, N = 1, R01=input$R01, R02=input$R02, 
     gamma = 365/input$Ip, Theta=input$theta, Xi=input$xi, Pi=input$pi)
paras["beta1"] = with(as.list(paras), R01*(365*gamma+mu))
paras["beta2"] = with(as.list(paras), R02*(365*gamma+mu))
start = c(S = 0.06, I1 = 0.001, I2 = 0.001, R1=0, R2=0, J1=0, J2=0, R = 0.938)

out = as.data.frame(ode(start, times, twostrain, paras))

#  plot(out$I1) 
sel=out$time>input$T[1]&out$time<input$T[2]

  plot(out$I1[sel], out$I2[sel], log=ifelse(input$lg==TRUE, "xy", ""), type="l", 
       xlab="I2", ylab="fraction infected")
  })

  })

