# SimSpikes UI.
# library(shiny)
# library(shinythemes)
# library(shinyjs)


# Define UI
fluidPage(shinyjs::useShinyjs(),  # Set up shinyjs,
          theme = shinythemes::shinytheme("lumen"),
                  titlePanel("Simulate spike sequences"),
                  sidebarLayout(
                    sidebarPanel(

                      # Button to go back to welcome screen and cancel.
                      shinyjs::disabled( actionButton("toggle", "Toggle help screen", width = '250px')),

                      # Select type of trend to plot
                      selectInput(inputId = "ISI", label = strong("ISI distribution:"),
                                  choices = c("Exponential","Gamma", "Inverse Gaussian", "Log Normal", "Weibull"),
                                  selected = "Gamma"),

                      conditionalPanel(condition = "input.ISI != 'Poisson'",
                                       numericInput(inputId = "hyper", label = strong("ISI parameter value:"), 10)
                      ),

                      # Select the end time of the simulation
                      numericInput(inputId = "end.time", label = strong("End time:"), 20),

                      # Select the refreactory period
                      numericInput(inputId = "t.min", label = strong("Refractory period:"), 0, step = 0.1, min = 0),

                      # Select the intensity function
                      textInput(inputId = "int.fn", label = strong("Intensity function:"), "cos(t)+1.5"),
                      conditionalPanel(
                        condition = "output.returnedFromFunction",
                        helpText(HTML("<font color = red> The intensity function must be positive! </font>"))
                      ),


                      # Select weather to save multiple spike sequences.
                      numericInput(inputId = "multi", label = strong("How many sequences?"), value = 1, min = 1, step =1),


                      # Add submit button.
                      hr(),
                      actionButton(inputId = "button", label = "Get Spikes"),
                      # Button for downloading
                      shinyjs::disabled( downloadButton("download", "Download", width = '250px'))
                    ),

                    # Output: Description, lineplot, and reference
                    mainPanel(
                      htmlOutput(outputId = "deets"),
                      plotOutput(outputId = "lineplot", height = "300px"),
                      htmlOutput(outputId = "desc")
                    )
                  )
  )
