ui = fluidPage(
   theme = bslib::bs_theme(version = 4, bootswatch = "simplex"),
   fluidRow(
   		column(1, style = "background-color:#023020;height:1250px;color: white",
   		radioButtons("radio1", label = h3("Scenarios"),
   		      choices = list("NED x LCR" = 'A', "NED x DM" = 'B'), 
   		      selected = 'A')),
   		      
    mainPanel(
      tabsetPanel(
        id = "tabset",
        
        tabPanel("Summary", plotOutput(outputId = "freqplot", click = "plot_click"),
							htmlOutput('x_value'),
        					tableOutput("selected_rows"),
							tableOutput("tab2"),
							selectInput("biotypeTable", label="Select Biotype (Top10 Avg)",	choices=LCR$Biotypes),
                			tableOutput("tab1")),

        tabPanel("BIOTYPE BOXPLOT", 
        				    selectInput("biotype", label="Select Biotype (Global Variation)", 
        					choices=LCR$Biotypes),
        				    plotOutput(outputId = "box", height=950)),
        
        tabPanel("SPECIFIC GENE SEARCH", 
        					textInput("search", "", 
        					placeholder = "GENE SYMBOL [e.g: ORM1]"),
        				    tableOutput(outputId = "pubmed")),  
      )
    )
  )
)
