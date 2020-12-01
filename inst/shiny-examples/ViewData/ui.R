# ViewData UI.

fluidPage(shinyjs::useShinyjs(),  # Set up shinyjs,
           theme = shinythemes::shinytheme("lumen"),
          titlePanel("View Cell Data"),
          sidebarLayout(
            sidebarPanel(

              # Button to go back to welcome screen and cancel.
              shinyjs::disabled( actionButton("toggle", "Toggle help screen", width = '250px')),
              # splitLayout(
              #   actionButton("toggle", "Toggle help screen", width = '250px'),
              #   actionButton("show", "Input Data", width = '100px')
              #    ),

              tags$hr(),
              # Text and button used for inputting data.
              helpText("Click the button below to load data"),
              actionButton("show", "Input Data", width = '250px'),

              tags$hr(),

              shinyjs::disabled(

                # Select the cell data you wish to view
                numericInput(inputId = "cell", label = strong("Cell number"), 1, min = 1),

                # Horizontal line ----
                tags$hr(),
                helpText(HTML("Change plot axes?")),
                # Select whether yoy want to change dimensions of the plot
                checkboxInput(inputId = "yaxis", label = strong("Default y axis?"), value = TRUE),

                conditionalPanel(condition = "input.yaxis == false",
                                 splitLayout(
                                   numericInput(inputId = "yMin", label = HTML("Min Ca<sup>2+</sup>  conc"), 0, step = 0.1),
                                   numericInput(inputId = "yMax", label = HTML("Max Ca<sup>2+</sup>  conc"), 1.6, step = 0.1)
                                 ),
                                 conditionalPanel(condition = "input.yMin >= input.yMax",
                                                  helpText(HTML("<font color = 'red'> Min must be smaller then max </font>"))
                                 )

                ),

                checkboxInput(inputId = "xaxis", label = strong("Default x axis?"), value = TRUE),

                conditionalPanel(condition = "input.xaxis == false",
                                 splitLayout(
                                   numericInput(inputId = "xMin", label = HTML("Min time"), 0, step = 2),
                                   numericInput(inputId = "xMax", label = HTML("Max time"), 1000, step = 2)
                                 ),
                                 conditionalPanel(condition = "input.xMin >= input.xMax",
                                                  helpText(HTML("<font color = 'red'> Min must be smaller then max </font>"))
                                 )
                ),
                tags$hr(),

                helpText(HTML("<font size = '+1'> Thresholding </font>")),

                # Add horizontal line at?
                textInput(inputId = "hline", label = strong("Threshold at? "), NULL),

                # Do we want transient periods
                checkboxInput(inputId = "trans1", label = strong("Transient period"), value = FALSE),
                conditionalPanel(condition = "input.trans1 == true",
                                 helpText("Transient area 1:"),
                                 splitLayout(
                                   helpText(HTML("<input type='button' value='0                                   ' disabled>")),
                                   numericInput(inputId = "trans1_high", label = NULL, value = NULL, step = 5)
                                 ),
                                 helpText("Transient area 2:"),
                                 splitLayout(
                                   numericInput(inputId= "trans2_low", label = NULL, value = NULL,step = 5),
                                   numericInput(inputId = "trans2_high", label = NULL, value = NULL, step = 5)
                                 ),
                                 helpText("Transient area 3:"),
                                 splitLayout(
                                   numericInput(inputId = "trans3", label = NULL, value = NULL, step = 5),
                                   helpText(HTML("<input type='button' value='End of Experient time                                  ' disabled>"))
                                 )

                                 # conditionalPanel(condition = "!isNaN(input.trans1) || input.trans1 >= input.trans2_high || input.trans1 >= input.trans2_low ",
                                 #                  helpText(HTML("<font color = 'red'> Transient regions must not overlap! </font>"))
                                 # )

                ),

                checkboxInput(inputId = "Lin.trend", label = strong("Remove linear trends?"), value = FALSE),


                # If we wnt to be ale to remove trends from only 1 not both sections

                checkboxInput(inputId = "doclose", label = strong("Limit closeness of spikes?"), value = F),

                conditionalPanel(condition = "input.doclose == true",
                                   numericInput(inputId = "close", label = HTML("closeness"), 1, step = 2)

                ),

                # conditionalPanel(condition = "Lin.trend == true",
                # checkboxInput(inputId = "Lin.trend1", label = strong("Remove linear trends in 1?"), value = FALSE),
                # checkboxInput(inputId = "Lin.trend2", label = strong("Remove linear trends in 2?"), value = FALSE)
                # ),

                # Button for downloading
                downloadButton("download", "Download", width = '250px')
                # hidden(actionButton("clearBrush", "Clear brush"))
              )),

            # Output: Description, lineplot, and reference
            mainPanel(
              htmlOutput(outputId = "desc"),
              plotOutput(outputId = "lineplot", height = "300px",  brush = "plot_brush"),
              plotOutput(outputId = "ISIplot", height = "300px"),
              htmlOutput(outputId = "Spikes")

            )

          )
)

