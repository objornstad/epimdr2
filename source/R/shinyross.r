#' Launch a shiny-app simulating a Ross-Macdonald model
#' @details
#' Launch app for details.
#' @examples
#' \dontrun{ross.app}
#' @export
# This creates the User Interface (UI)
ross.app=shinyApp(
ui <- pageWithSidebar(
headerPanel(""),
sidebarPanel(
sliderInput("ip",  HTML("Human infectious period (1/&gamma;)"), 7,
              min = 0, max = 30),
sliderInput("g", "Gonotrophic cycle length (1/a)", 4,
              min = 1, max = 20),
sliderInput("b", "V-to-H transmission probabiliy (b)", 0.1,
              min = 0, max = 1),
sliderInput("c", "H-to-V transmission probability (c)", 1,
              min = 0, max = 1),
sliderInput("vl", HTML("Vector longevity (1/&mu;)"), 7,
              min = 0, max = 21),
sliderInput("m", "Numbers of vectors per host (m)", 50,
              min = 1, max = 200)
),
mainPanel(
  #Multiple tabs in main panel
  tabsetPanel(
      #Tab 1: Time plot (plot1 from server)
      tabPanel("Phase plane", plotOutput("plot1")), 
      #Tab 2: Phase plot (plot2 from server)
      tabPanel("Equations", 
           withMathJax(
             helpText("Host infected fraction $$\\frac{dx}{dt} = \\underbrace{(a b m) y (1-x)}_{\\mbox{V-to-H transmission}} - \\underbrace{\\gamma x}_{\\mbox{H recovery}}$$"),
            helpText("Vector infected fraction $$\\frac{dy}{dt} = \\underbrace{a c x (1-y)}_{\\mbox{H-to-V transmission}}-  \\underbrace{\\mu y}_{\\mbox{V death}}$$"),
            helpText("Reproduction number $$R_0 =\\frac{a c}{\\gamma} \\frac{a b m}{\\mu} = \\frac{a^2 c b m}{\\gamma \\mu}$$"),
            helpText("Endemic equilibrium $$x^*=(R_0-1)/(R_0 + a c / \\mu)$$, $$y^*=\\frac{R_0 - 1}{R_0}\\frac{a c/\\mu}{1+ac/\\mu}$$"),
            helpText("Host isocline $$y =\\frac{\\gamma x}{(a b m)(1-x)}$$"),
            helpText("Vector isocline $$y = \\frac{c x}{a c x + \\mu}$$"),
            helpText("a: biting rate; b: host-to-vector transmission probability; c: vector-to-host transmission probability; 1/gamma: host infectious period; 1/mu: vector longevity"),
            helpText("REFERENCE: Aron, J. L., & May, R. M. (1982). The population dynamics of malaria. In The population dynamics of infectious diseases: theory and applications (pp. 139-179). Springer, Boston, MA.")

            ))
  ))), #End of ui()

# This creates the 'behind the scenes' code (Server)
server <- function(input, output) {
  require("phaseR", "deSolve")
grfn=function(t, y, parameters){
  X=y[1]
  Y=y[2]
  with(as.list(parameters),{
  dx=a*b*m*Y*(1-X)-gamma*X
  dy=a*c*X*(1-Y)-mu*Y
  gr=c(dx, dy)
  list(gr)
  })
}


  output$plot1 <- renderPlot({
times=seq(0, 365*2, by=.1)

parameters  = c(gamma = 1/input$ip, a =  1/input$g, b=input$b, c=input$c, mu=1/input$vl, m=input$m)
start=c(0.01, 0.01)

out=ode(y=start,
  times=times,
  func=grfn,
  parms=parameters)

  out=as.data.frame(out)

with(as.list(parameters),{
curve(gamma*x/((a*b*m)*(1-x)), 0,1, ylim=c(0,1), xlab="Human prevalence (x)", ylab="Mosquito prevalence (y)")
R0=m*a^2*b*c/(mu*gamma)
title(paste ("R0=",round(R0,2)))
curve(a*c*x/(a*c*x+mu), 0,1, add=TRUE, col="red")
fld=flowField(grfn, xlim=c(0,1), ylim=c(0,1), 
parameters=parameters, system="two.dim", add=TRUE,
ylab="H", xlab="M")

})
points(out[,2], out[, 3])
legend("topleft", c("H isocline", "M isocline", "Trajectory"), lty=c(1,1,0), col=c(1,2, 1), pch=c(NA,NA, 1))
   })
  })

