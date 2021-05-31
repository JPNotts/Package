## Server for SimSpikes.

rm(list = ls())
# library(shiny)
# library(shinyjs)
# library(shinythemes)


# source the code required to simulate spike sequences.
# source("SimSpikes.R")

# Define server function
server <- function(input, output,session) {
  # Global toggle to switch between start screen and spikes.
  toggle <- FALSE

  # observe toggle and swap view screen
  observeEvent(input$toggle, {
    if(toggle == TRUE){
      toggle <<- FALSE
      shinyjs::enable("int.fn"); shinyjs::enable("ISI") ; shinyjs::enable("hyper"); shinyjs::enable("end.time"); shinyjs::enable("t.min"); shinyjs::enable("multi")
      shinyjs::enable("button"); shinyjs::enable("download")

    }
    else{
      toggle <<- TRUE
      shinyjs::disable("int.fn"); shinyjs::disable("ISI") ; shinyjs::disable("hyper"); shinyjs::disable("end.time"); shinyjs::disable("t.min"); shinyjs::disable("multi")
      shinyjs::disable("button"); shinyjs::disable("download")
    }
  })

  # Check if intensity will be positive.
  output$returnedFromFunction <- reactive({
    t <- seq(0,input$end.time,input$end.time/8000)
    ans <- eval(parse(text = input$int.fn))
    if(any(ans<0)){ return(TRUE)}
    else{ return(FALSE)}
  })

   outputOptions(output, "returnedFromFunction", suspendWhenHidden = FALSE)

  int.fn <- eventReactive(input$button, {
    t <- seq(0,input$end.time,input$end.time/8000)
    ans <- eval(parse(text = input$int.fn))
    if(any(ans<0))(stop("intensity function must be positive!"))

    return(ans)
  })

  ISI <- eventReactive(input$button, {
    if(input$ISI == "Exponential"){ISI.type = "Exponential"}
    if(input$ISI == "Gamma"){ISI.type = "Gamma"}
    if(input$ISI == "Inverse Gaussian"){ISI.type = "InverseGaussian"}
    if(input$ISI == "Log Normal"){ISI.type = "LogNormal"}
    if(input$ISI == "Weibull"){ISI.type = "Weibull"}
    return(ISI.type)
  })

  get.spikes <- eventReactive(input$button, {
    shinyjs::enable("download")
    # Get the inputs
    end.time <- input$end.time ; t.min <- input$t.min ; ISI <- input$ISI ; multi <- input$multi ; int.fn <- input$int.fn ; hyper <- input$hyper
    d.spikes <- SimSpikes::simulate_spikes(input$end.time, int.fn(), input$hyper, steps =2000, T.min = input$t.min, ISI.type  = ISI(), multi = input$multi, add.end = TRUE, do.log = T)
    return(list(spikes = d.spikes,  end.time =end.time, t.min=t.min, ISI=ISI, multi=multi, int.fn=int.fn, hyper = hyper))
  })

  check <- eventReactive(input$button, {
    message <- NULL
    shinyjs::enable("toggle")
    # Firstly check that the intensity function won't create too many spikes.
    int <- sum(int.fn()*input$end.time/8000)
    if(int/2000 > 0.3){
      message <- c(message,paste("Warning! The expected number of spikes is", round(int,digits = 2)," this may result in numerical errors as the discritisation isn't fine enough."))
    }
    return(message)
  })

  output$lineplot <- renderPlot({
    input$toggle
    if(input$button == 0 || toggle)(return(NULL))

    # Get the parameters we need but isolate to stop autochanging.
    multi = isolate(input$multi)
    end.time <- isolate(input$end.time)

    t <- seq(0,end.time,end.time/8000)

    y.change <- max(0,min(int.fn())-0.1)
    y.low <- - (min(multi,20)/20)*(max(int.fn()+0.1))
    color = "#434343"
    par(mar = c(4, 4, 1, 1))

    #  now for y axis create label points
    plot(x =  t, y =  int.fn(), type = "l",
         xlab = "Time", ylab = "Intensity", xlim = c(0,end.time), ylim = c(y.low,max(int.fn()+0.1)),yaxt = "n", col = color, fg = color, col.lab = color, col.axis = color)
    axis(2,at = seq(0,round(max(int.fn()+0.1), digit = 2),round(max(int.fn()+0.1))/5), labels = seq(0,round(max(int.fn()+0.1), digit = 2),round(max(int.fn()+0.1))/5) )

    lines(c(0,end.time),c(0,0), col = "grey")
    for(i in 1:min(20,multi)){
      d.spikes <- get.spikes()$spikes
      spikes.cur <- d.spikes[!is.na(d.spikes[,i]),i]
      spikes.cur <- spikes.cur[-length(spikes.cur)]  # Remove the end time.
      y.height <-  - (i/20)*(max(int.fn()+0.1))
      y.axis <- rep(y.height,length(spikes.cur))
      points(spikes.cur, y.axis, pch = 20)
    }


  })

  output$desc <- renderUI({
    # Initial message when someone opens the app.
    input$toggle
    if(input$button == 0 || toggle){
      return(NULL)
    }

    multi = isolate(input$multi)
    mess <- check()
    d.spikes <- get.spikes()$spikes
    spikes.cur <- d.spikes[!is.na(d.spikes[,1]),1]
    spikes.cur <- spikes.cur[-length(spikes.cur)]  # Remove the end time.
    spikes <- round(spikes.cur,digits = 2)
    if(multi == 1){
      if(length(spikes) == 0){
        message <- "There are no spikes in the sequence."
      }
      if(length(spikes) == 1){
        message <- paste("The spike time is: ",spikes[1],".")
      }
      if(length(spikes) == 2){

        message <- (paste("The spike times are: ",spikes[1],",",spikes[2],"."))
      }
      if(length(spikes)>2){
        message <- "The spike times are: <font size ='-1'> "
        for(i in 2:(length(spikes)-1)){message <- c(message, paste(spikes[i],", "))}
        message <- c(message, paste(spikes[length(spikes)],". </font>"))
      }
      all <- paste(message, collapse= "")
      all <- c("<hr>",check(), all)

    }
    if(multi > 1){
      message <- "Below are the spike times for first spike sequence only. <br>"
      if(length(spikes) == 0){
        message <- c(message,"There are no spikes in the sequence.")
      }
      if(length(spikes) == 1){
        message <- c(message,paste("The spike time is: ",spikes[1],"."))
      }
      if(length(spikes) == 2){

        message <- c(message,(paste("The spike times are: ",spikes[1],",",spikes[2],".")))
      }
      if(length(spikes)>2){
        message <- c(message,"The spike times are: <font size ='-1'> ")
        for(i in 2:(length(spikes)-1)){message <- c(message, paste(spikes[i],", "))}
        message <- c(message, paste(spikes[length(spikes)],". </font>"))
      }
      all <- paste(message, collapse= "")
      all <- c("<hr>",check(), all)
    }
    HTML(paste(all, collapse = "<br/>"))
  })

  output$deets <- renderUI({
    input$toggle
    # Initial message when someone opens the app.
    # cat("toggle is",toggle,"\n")
    # cat("input button is is",input$button,"\n")
    if(input$button == 0 || toggle){
      message <- "<p> <font size = '+3'> Welcome to Sim Spikes </font> </p>
<p> This app is for simulating spike sequences, from an inhomogeneous point process, which is described by its interspike interval (ISI) distribution and intensity function.  </p>
      <p> <font size = '+1'> Using the app </font> </p>
      <p> Use the sidebar to choose the model parameters desired, furthermore you can choose to simulate multiple spike sequences. Once you click run, the main panel will plot the intensity function and will show the spike sequences (max 20) below the x-axis. One simulated spike sequence will be wrote in the panel.</p>
      <p> <font size = '+1'> Inputting the intensity function </font> </p>
      <p> Write the intensity function in the given box in the sidebar. You <font color = 'red'> must </font> write this as a function of t and explicitly write all operations.  <br> TIP: You can write completed functions in the box, the method only parses the text into R then applies the vector t. E.g for a step function TRY: c(t[t<10]*0 +2.4, t[t>=10]*0 +0.5) <p>
      <p> <font size = '+1'> Downloading the results </font> </p>
      <p> You can downlad the results by clicking the 'DOWNLOAD button'. This will save a file called SimSpikes.zip to your computer, which will contain the file details.txt containing the model parameters used to simulate the spike sequences, and the file spikes.csv which contains the spike sequences - columns in the spreadsheet.  <p>
      <p> If you experience any issues you can email: pmxjp8@nottingham.ac.uk </p>"
      return(HTML(message))
    }

    return(HTML("<p> <font size = '+2'>Plot of the intensity function and spikes: </font> </p>"))
  })


  # Download the current files.
  output$download <- downloadHandler(
    filename = function(){
      text <- paste0(c("SimSpikes"))
      paste0(text,".zip")

    },
    content = function(file){
      #go to a temp dir to avoid permission issues
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      files <- NULL;

      # Next create the files to put into the zipped folder.

      # Create file details.txt to store information of what parameter values were simulated from.
      fileConn<-file(paste0("Details.txt"), open ='wt')
      writeLines(c(paste0("The inputs for the simulated spike sequences were"),
                   paste0("End time = ",get.spikes()$end.time),
                   paste0("intensity function  = ",get.spikes()$int.fn),
                   paste0("ISI distribution is ",get.spikes()$ISI, " with parameter ", get.spikes()$hyper),
                   paste0("Refractory period was set to ",get.spikes()$t.min),
                   paste0("Number of sequences is ", get.spikes()$multi,".")), fileConn)
      close(fileConn)


      # Save the spikes in csv file.
      d.spikes <- get.spikes()$spikes
      write.csv(d.spikes, file = "spikes.csv")

      #create the zip file
      files <- c("Details.txt","spikes.csv")
      zip(file,files)
    }
  )
}
