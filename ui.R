
library(shiny)

# Define UI for application that draws a histogram
fixedPage(

    # Application title
    titlePanel("Planning cancer screening trials for the late-stage
cancer incidence reduction"),

    # Sidebar with a slider input for number of bins
    tabsetPanel(
      tabPanel('Introduction',
               includeHTML("intro.html")
               ),
      tabPanel("Cancer natural history inputs",
               withMathJax(),
               sidebarLayout(
                 sidebarPanel(
               #numericInput("w", "overall cancer incident rate", 0.0067, min = 0,max=1),
               numericInput("I.e", "early-stage cancer incident rate \\((I_{E})\\)", 0.0017, min = 0,max=1),
               numericInput("I.l", "late-stage cancer incident rate \\((I_{L})\\)", 0.005, min = 0,max=1),
               numericInput("P.E", "prevalence of early-stage preclinical cancer \\((P^{EP})\\)", 0.018, min = 0,max=1),
               numericInput("P.L", "prevalence of late-stage preclinical cancer \\((P^{LP})\\)", 0.002, min = 0,max=1),
                 ),
               mainPanel(includeHTML('stable_disease.html'),
                         uiOutput('natural.history'),)
               )
      ),
      tabPanel("Screening program inputs & power output", 
               
               sidebarLayout(
                 sidebarPanel(
                   numericInput("se.e.exp", "exp. arm: Sensitivity of detecting early-stage preclinical cancer", 0.51, min = 0,max=1,step=0.05),
                   numericInput("se.l.exp", "exp. arm: Sensitivity of detecting late-stage preclinical cancer", 0.85, min = 0,max=1,step=0.05),
                   numericInput("se.e.ctr", "control arm: Sensitivity of detecting early-stage preclinical cancer", 0.15, min = 0,max=1,step=0.05),
                   numericInput("se.l.ctr", "control arm: Sensitivity of detecting late-stage preclinical cancer", 0.25, min = 0,max=1,step=0.05),
                   numericInput('exp.ctam','exp. arm: contamination rate',value=0.15,min=0,max=1,step=0.1),
                   numericInput('ct.ctam','control arm: contamination rate',value=0.2,min=0,max=1,step=0.1),
                   numericInput('non.scr','never screening rate',value=0.1,min=0,max=1,step=0.1),
                   numericInput("n.test",'number of tests',value=3,min=1,step=1),
                   numericInput("t.interval",'intervals between tests (year)',value=1,min=0.5,step=0.5),
                   numericInput("a.yr",'uniform accural period (year)',value=2,min=0,step=0.5),
                   numericInput("fu",'planned analysis time (years after the 1st enrollment)',value=5,min=0.5,step=0.5),
                   numericInput("n.ctr",'# individuals in control arm',value=15000,min=0,step=5000),
                   numericInput("n.exp",'# individuals in exp arm',value=15000,min=0,step=5000),
                   numericInput('alpha','one-sided Type-I error',value=0.05,min=0,max=1,step=0.01)
                 ),
                 mainPanel(
                           uiOutput('power.rst'),
                           plotOutput("cum.prob.plot"),
                           plotOutput("rel.red.plot"),
                           plotOutput("power.plot")
                           
                           )
               )
               

      )
    )
    
    
)
