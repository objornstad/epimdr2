#' Launch a shiny-app simulating a two-strain SIR model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{twostrain.app}
#' @export

twostrain.app=shinyApp(
ui = pageWithSidebar(
  headerPanel("A two-strain SIR model"),
  sidebarPanel(
    sliderInput("beta1", "Transmission 1 (yr^-1):", 500,
                min = 0, max = 3000),
    sliderInput("beta2", "Transmission 2 (yr^-1):", 750,
                min = 0, max = 3000),
    sliderInput("Ip", "Infectious period (days)", 5,
                min = 1, max = 100),
    sliderInput("mu", "birth rate (per year):", 0.02,
                min = 0, max = .1),
    sliderInput("theta", "Theta:", 0.5,
                min = 0, max = 1),
    sliderInput("xi", "Xi:", 0.5,
                min = 0, max = 1),
    sliderInput("pi", "Pi:", 0.5,
                min = 0, max = 1),
    sliderInput("T", "Time range:",
                min = 0, max = 30, value = c(10,30)),
    #checkboxInput("lg", "Log", FALSE), width=3),
    checkboxInput("lg", "4th root", FALSE), width=3),
  mainPanel(
    tabsetPanel(
      tabPanel("Parameters", dataTableOutput("table1")), 
      tabPanel("Time", plotOutput("plot1")),
      tabPanel("S* plot", plotOutput("plot2")),
      tabPanel("Phase plane", plotOutput("plot3")),
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


#R_{0,1}=&\frac{\beta}{\gamma+\mu} \frac{N}{N}\\
#S_1^* =& 1/R_0\\
#I_1^* =& \mu (R_0 - 1)/\beta\\
#R_1^* =& \frac{\gamma I_1^*}{(1-\Pi) + mu}\\
#R^* =& (1 - \Pi) * R_1^* / \mu\\
#Q_0=\frac{\beta_2}{\gamma+\mu} \frac{S^*}{N} + \frac{\Xi \beta_2 R_1^* I_2}{N (\gamma+\mu)}


# This creates the 'behind the scenes' code (Server)
server = function(input, output) {
  output$table1 = renderDataTable({
    R01 = (input$beta1/365)/(1/input$Ip + input$mu/365)
    S1star=1/R01
    I1star=input$mu/365 *(R01-1)/(input$beta1/365)
    R1star=1-(S1star+I1star)
    #Rstar= (1-input$pi)*R1star/input$mu/365
    R02 = (input$beta2/365)/(1/input$Ip + input$mu/365)
    Q02 = ((input$beta2/365)*S1star)/(input$mu/365 + 1/input$Ip) + (input$beta2/365)*input$xi*R1star/(input$mu/365 + 1/input$Ip) 
    S2star=1/R02
    I2star=input$mu/365 *(R02-1)/(input$beta1/365)
    R2star=1-(S2star+I2star)
    #Rstar= 1-(S1star+I1star+S1star+R1star)
    r=log(Q02)/input$Ip
    Td=log(2)/r
    data.frame(R01=round(R01,2), S1star=round(S1star,3), I1star=round(I1star,4),R02=round(R02,2), Q02=round(Q02,2),
               S2star=round(S2star,3), I2star=round(I2star,4), Td=round(Td,1))
    #  data.frame(R01=round(R01,2), S1star=round(S1star,3), I1star=round(I1star,4),R1star=round(R1star,4), R02=round(R02,2), Q02=round(Q02,2),
    #             S2star=round(S2star,3), I2star=round(I2star,4),R2star=round(R2star,4), Td=round(Td,1))
  })
  output$plot1 <- renderPlot({
    times  = seq(0, input$T[2], by=1/100)
    paras  = c(mu = input$mu, N = 1, beta1=input$beta1, beta2=input$beta2, 
               gamma = 365/input$Ip, Theta=input$theta, Xi=input$xi, Pi=input$pi)
    #paras["beta1"] = with(as.list(paras), R01*(gamma+mu))
    #paras["beta2"] = with(as.list(paras), R02*(gamma+mu))
    start = c(S = 0.06, I1 = 0.001, I2 = 0.001, R1=0, R2=0, J1=0, J2=0, R = 0.938)
    
    out = as.data.frame(ode(start, times, twostrain, paras))
    
    sel=out$time>input$T[1]&out$time<input$T[2]
    
    par(mar = c(5,5,2,5))
    if(input$lg==FALSE){
      plot(x=out$time[sel], y=out$I2[sel], ylab="fraction", xlab="time", type="l",
           ylim=range(out[sel, c(3,4,7,8)]), xlim=c(input$T[1], input$T[2]), log=ifelse(input$lg==TRUE, "y", ""), col="red")
      lines(x=out$time, y=out$I1, col="black")
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
    }
    
    if(input$lg==TRUE){
      plot(x=out$time[sel], y=out$I2[sel]^(1/4), ylab="fourth root fraction", xlab="time", type="l",
           ylim=range(out[sel, c(3,4,7,8)]^(1/4)), xlim=c(input$T[1], input$T[2]), col="red")
      lines(x=out$time, y=out$I1^(1/4), col="black")
      lines(x=out$time, y=out$J1^(1/4), lty=2)
      lines(x=out$time, y=out$J2^(1/4), col="red", lty=2)
      
      par(new=T)
      #plot(x=out$time, y=out$S, type="l", col="green", axes=FALSE, xlab=NA, ylab=NA 
      #    ylim=range(out[sel,2]), xlim=c(input$T[1], input$T[2]), #log=ifelse(input$lg==TRUE, "y", ""))
      #axis(side = 4, col="green")
      #mtext(side = 4, line = 4, "S", col="green")
      legend("right",
             legend=c("I1", "I2", "J1", "J2"),
             lty=c(1,1,2,2),
             col=c("black", "red", "black", "red"))
    }
  })
  
  # output$plot2b <- renderPlot({
  #   times  = seq(0, input$T[2], by=1/200)
  #   paras  = c(mu = input$mu, N = 1, beta1=input$beta1, beta2=input$beta2, 
  #              gamma = 365/input$Ip, Theta=input$theta, Xi=input$xi, Pi=input$pi)
  #   #paras["beta1"] = with(as.list(paras), R01*(gamma+mu))
  #   #paras["beta2"] = with(as.list(paras), R02*(gamma+mu))
  #   start = c(S = 0.999, I1 = 0.001, I2 = 0.00, R1=0, R2=0, J1=0, J2=0, R = 0)
  #   
  #   out1 = as.data.frame(ode(start, times, twostrain, paras))
  #   ta=out1[out1[,1]<input$T[1],]
  #   #ta=tail(out1, 1)
  #   start2 = c(S = ta[1,2], I1 = ta[1,3], I2 = 0.001, R1=ta[1,5], R2=ta[1,6], J1=ta[1,7], J2=ta[1,8], R = ta[1,9])
  #   out2 = as.data.frame(ode(start2, times, twostrain, paras))
  #   
  #   sel=out1$time>0
  #   #sel=TRUE
  #   R01 = (input$beta1/365)/(1/input$Ip + input$mu/365)
  #   R02 = (input$beta2/365)/(1/input$Ip + input$mu/365)
  #   plot(out1$S[sel], R01*out1$S[sel], ylim=c(0, R01), type="l", xlab="S", ylab="Re")
  #   lines(out2$S[sel], R01*out2$S[sel], lty=2)
  #   lines(out2$S[sel], R02*out2$S[sel], col=2)
  #   abline(h=1)
  #   
  #   legend("right",
  #          legend=c("Re1", "Re2"),
  #          col=c("black", "red"))
  # })
  # 
  output$plot2 <- renderPlot({
    times  = seq(0, input$T[2], by=1/100)
    paras  = c(mu = input$mu, N = 1, beta1=input$beta1, beta2=input$beta2, 
               gamma = 365/input$Ip, Theta=input$theta, Xi=input$xi, Pi=input$pi)
    #paras["beta1"] = with(as.list(paras), R01*(gamma+mu))
    #paras["beta2"] = with(as.list(paras), R02*(gamma+mu))
    start = c(S = 0.999, I1 = 0.001, I2 = 0.00, R1=0, R2=0, J1=0, J2=0, R = 0)
    
    out1 = as.data.frame(ode(start, times, twostrain, paras))
    ta=out1[out1[,1]<input$T[1],]
    #ta=tail(out1, 1)
    start2 = c(S = ta[1,2], I1 = ta[1,3], I2 = 0.001, R1=ta[1,5], R2=ta[1,6], J1=ta[1,7], J2=ta[1,8], R = ta[1,9])
    out2 = as.data.frame(ode(start2, times, twostrain, paras))
    
    sel=out1$time>0
    #sel=TRUE
    R01 = (input$beta1/365)/(1/input$Ip + input$mu/365)
    R02 = (input$beta2/365)/(1/input$Ip + input$mu/365)
    plot(out1$S[sel], R01*out1$S[sel], ylim=c(0, R01), type="l", xlab="S", ylab="Re")
    lines(out2$S[sel], R01*out2$S[sel], lty=2)
    lines(out2$S[sel], R02*out2$S[sel], col=2)
    abline(h=1)
    points(1/R01, 1, pch="X")
    points(1/R02, 1, pch="O")
    
    legend("right",
           legend=c("Re1", "Re2"),
           lty=c(1,1),
           col=c("black", "red"))
  })
  
  output$plot3 <- renderPlot({
    times  = seq(0, input$T[2], by=1/200)
    paras  = c(mu = input$mu, N = 1, beta1=input$beta1, beta2=input$beta2, 
               gamma = 365/input$Ip, Theta=input$theta, Xi=input$xi, Pi=input$pi)
    #paras["beta1"] = with(as.list(paras), R01*(gamma+mu))
    #paras["beta2"] = with(as.list(paras), R02*(gamma+mu))
    start = c(S = 0.999, I1 = 0.001, I2 = 0.00, R1=0, R2=0, J1=0, J2=0, R = 0)
    
    out1 = as.data.frame(ode(start, times, twostrain, paras))
    ta=out1[out1[,1]<input$T[1],]
    #ta=tail(out1, 1)
    start2 = c(S = ta[1,2], I1 = ta[1,3], I2 = 0.001, R1=ta[1,5], R2=ta[1,6], J1=ta[1,7], J2=ta[1,8], R = ta[1,9])
    out2 = as.data.frame(ode(start2, times, twostrain, paras))
    
    sel=out1$time>0
    #sel=TRUE
    R01 = (input$beta1/365)/(1/input$Ip + input$mu/365)
    R02 = (input$beta2/365)/(1/input$Ip + input$mu/365)
    plot(out2$S[sel], out1$I1[sel], log="", xlim=c(0, 1), type="l", xlab="S", ylab="I")
    lines(out2$S[sel], out2$I1[sel], lty=2)
    lines(out2$S[sel], out2$I2[sel], col=2)
    # abline(h=1)
    
    legend("right",
           legend=c("I1only", "I2", "I1tostrain"),
           lty=c(1,1,2),
           col=c("black", "red", "black"))
  })
}
)
