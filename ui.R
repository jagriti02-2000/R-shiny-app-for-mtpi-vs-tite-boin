library(shiny)
library(highcharter)
library(tidyverse)

ui <- fluidPage(tags$head(tags$style(
  HTML('
         #sidebar {
            background-color: #34495E;
         }
         
          .form-group {
            margin-bottom: 1% !important;
          }
        
          label { font-size:80%; font-family:Times New Roman; margin-bottom: 
    02px; color : #F7F9F9}

        body, label, input, button, select { 
          font-family: "Arial";
        }
       
       #sim {height : 20px;}
       #cohort {height : 20px;}
       #dtl_window {height : 20px}
        #no_dose {height : 20px}
       #dall {height : 20px}
       #max_sample {height : 20px}
       #max_sample_dose {height : 20px}
       #dall {height : 20px}
       #ac {height : 20px}'
       
  )
)),

titlePanel("Compartive Review of Phase 1 Oncology trials"),
tags$style(HTML("
    body {
            background-color: F2F3F4;
            color: 34495E;
            }")),

tabsetPanel(
  
  tabPanel("Simulation",
           sidebarLayout(
             sidebarPanel(id="sidebar",width=3,
                          sliderInput("pT", "Target Toxicity Probability", value = 0.3, min = 0, max = 1),
                          numericInput("sim", "Number of simulation", value = 10,min=10),
                          numericInput("cohort", "Cohort Size", value = 3,min=3),
                          numericInput("dtl_window", "DLT window(In weeks)", value = 8,min = 4),
                          numericInput("no_dose", "No of doses", value = 3,min = 3),
                          textInput("dall","True Toxicity Probabilty",value = NULL,placeholder = "0.1,0.3,0.5"),
                          
                          numericInput("max_sample","Maximum sample size",value = 50,min=20),
                          numericInput("max_sample_dose","Maximum sample size at a particular dose",value = 9,min=6),
                          
                          numericInput("ac","Accrual Rate (n patients/m weeks)",value = 3,min=1),
                          actionButton("simulate", "Simulate!")
             ),
             mainPanel(
               
               tabsetPanel(
                 
                 tabPanel("Table",
                          tableOutput("table")),
                 tabPanel("Plots",
                          highchartOutput("plot")),
                 
               )
             )
             
           )),
  
  tabPanel("Analysis",sidebarLayout(
    sidebarPanel(id="sidebar",widh=3,
                 
                 selectInput("design_type", "Select Design",c("TITE-BOIN","MTPI")),
                 
                 
                 #Only show this when design MTPI
                 conditionalPanel(condition = "input.design_type=='MTPI'" ,
                                  numericInput("no_pat", "Number of patients at the current dose",value=0,min = 0),
                                  numericInput("no_dlt", "Number of patients with DLT at the current dose",value=0,min = 0),
                                  numericInput("e1", "Value of Epsilon 1",value=0.05,min = 0),
                                  numericInput("e2", "Value of Epsilon 2",value=0.05,min = 0),
                                  numericInput("alpha", "Prior Alpha", value = 1,min=1),
                                  numericInput("beta", "Prior Beta", value = 1,min=1),
                                  actionButton("Analyze", "Analyze!")
                 ),
                 conditionalPanel(condition = "input.design_type=='TITE-BOIN'" ,
                                  sliderInput("pT", "Target Toxicity Probability", value = 0.3, min = 0, max = 1),
                                  numericInput("no_pat1", "Number of patients at the current dose",value=3,min = 0),
                                  numericInput("no_dlt1", "Number of patients with DLT at the current dose",value=0,min = 0),
                                  numericInput("comp", "Number of patients with complete data at the current dose",value=0,min = 0),
                                  numericInput("pending", "Number of pending patients at the current dose",value=0,min = 0),
                                  numericInput("dlt_win", " DTL Assesment time",value=8,min = 4),
                                  textInput("time","Follow up time for pending patients at the current dose (in weeks)",placeholder = "1,2,3"),
                                  actionButton("Run", "Run!")
                                  
                 ),
                 
                 
                 
    ),
    
    mainPanel(
      uiOutput("text"),
      uiOutput("image_dose")
      
    )
    
  )
  
  )
  
)
)

