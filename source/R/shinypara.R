#' Launch a shiny-app simulating the spatially-extended Nicholson-Bailey parasitoid model
#' @details
#' Launch app for details
#' @examples
#' \dontrun{nbspat.app}
#' @export
#' @importFrom phaseR flowField
#' @importFrom ggplot2 ggplot
#' @importFrom plotly ggplotly
nbspat.app=shinyApp(
# This creates the User Interface (UI)
ui <- fluidPage(

  # App title ----
  titlePanel("Spatial Nicholson-Bailey parasitoid model"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar to demonstrate various slider options ----
    sidebarPanel(
      sliderInput("Dh", "Host dispersal:", 0.5,
              min = 0, max = 1),
      sliderInput("Dp", "Parasitoid dispersal:", 0.5,
              min = 0, max = 1),
      sliderInput("R", "Host growth:", 2,
              min = 1, max = 5, step=0.1),
      sliderInput("a", "Parasitoid efficiency:", 1,
              min = 0, max = 5),
      numericInput("IT", "Generations:", 200,
               min = 0, max = 1000)
    ),

    # Main panel for displaying outputs ----
 mainPanel(tabsetPanel(
  tabPanel("Simulation", plotlyOutput("plot1", height = 500)),
     tabPanel("Details", 
     withMathJax(
         helpText("Local dynamics:"),
             helpText("Host $$H_t = R H_{t-1} exp(-a P_{t-1})$$"),
          helpText("Parasitoid $$P_t = R H_{t-1} (1- exp(a P_{t-1})$$"),
          helpText("Redistribution"),
          helpText("$$\\vec{H_t} = \\vec{H_t} \\vec{k_h}$$"),
          helpText("$$\\vec{P_t} = \\vec{P_t} \\vec{k_p}$$"),          
        helpText("where k_h and k_p are redistribution matrices."),          
          helpText("REFERENCE: Bjornstad, O. N. and Bascompte, J. (2001). Synchrony and second-order spatial correlation
in host-parasitoid systems. Journal of Animal Ecology, 70(6), 924-933.

Hassell, M. P., Comins, H. N. and May, R. M. (1991). Spatial structure and chaos in
insect population dynamics. Nature, 353(6341), 255-258.

Nicholson, A. J. and Bailey, V. A. (1935). The balance of animal populations. Part I.
Proceedings of the Zoological Society of London, 105(3), 551-598.")))
           )
)
)),

server = function(input, output) {

#xlen is width of the lattice (E-W)
#ylen is height of the lattice (N-S)
xlen = 30
ylen = 30

hp.dyn = function(h, p, R, a){ 
   #hnew is the post-interaction host density
   hnew = R * h * exp(- a * p)
   #pnew is the post-interaction parasitoid density
   pnew = R * h * (1 - exp(- a * p))
   #the two vectors of results are stored in a "list"
   res = list(h = hnew, p = pnew)
   return(res)
} 

xy = expand.grid(1:xlen, 1:ylen)
dmat = as.matrix(dist(xy))
  # Show the values in an HTML table ----
output$plot1 <- renderPlotly({
kh = ifelse(dmat<1.5,input$Dh/8,0)
kp = ifelse(dmat<1.5,input$Dp/8,0)
diag(kh) = 1-input$Dh
diag(kp) = 1-input$Dp

IT = input$IT
hmat = matrix(NA, nrow=xlen*ylen, ncol = IT)
pmat = matrix(NA, nrow=xlen*ylen, ncol = IT)
hmat[,1] = 0
pmat[,1] = 0
hmat[23,1] = 4
pmat[23,1] = 1

for(i in 2:IT){
   #growth
   tmp = hp.dyn(h = hmat[,(i-1)], p = pmat[,(i-1)], 
      R = input$R, a = input$a)
   #redistribution
   hmat[,i] = tmp$h%*%kh;
   pmat[,i] = tmp$p%*%kp;
}

hmat2=as.data.frame(hmat)
hmat2$x=xy[,1]
hmat2$y=xy[,2]
longH=reshape(hmat2, direction="long", varying=1:IT, v.names="H")
anim=ggplot(longH, aes(x=x, y=y, frame=longH$time))+geom_point(size=longH$H)
ggplotly(anim, width = 600, height = 600) %>%
    animation_opts(200)
})

}
)