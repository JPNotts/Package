## Server for ViewData

shinyjs::useShinyjs()  # Set up shinyjs

function(input, output, session) {

  # Create global objects:
  #     - store: to store the values of threshold, transients and rm linear trend for each cell.
  #     - df: data frame that stores the spikes for each cell after thresholding.
  #     - folder: stores the system time used to create the folder of the outputs in Output/ViewData/
  store <- 1
  df <- 1
  folder <- NULL
  err <- FALSE
  toggle <- FALSE

  # Create time stores to check whether new input for data and load.
  d.old <- Sys.time()
  d.new <- d.old
  l.old <- d.old
  l.new <- d.old



  # reactiveValues object for storing current data set.
  vals <- reactiveValues(data = NULL)

  # The next section deals with creating the modal and errors with inputting the data/loading transients.
  {
    # Return the UI for a modal dialog with data selection input. If 'failed' is
    # TRUE, then display a message that the previous value was invalid.
    dataModal <- function(failed = FALSE, mg = "Something has gone wrong!") {
      modalDialog(
        span(HTML("<p>Please select the file containing the Ca <sup>2+</sup> concentration data.  <br> <font size = '-2'>Required to be .csv type. </font>  </p>")),
        fileInput(inputId = "data", label = "", placeholder = "",  accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv")),
        span(HTML("<p> If you want to load transients details please load below. <br> <font size = '-2'>This must be a .csv file with the correct dimensions to match the data. </font>  </p>")),

        fileInput(inputId = "load", label = "", placeholder = ""),
        if (failed)
          div(tags$b(HTML(mg), style = "color: red;")),

        footer = tagList(
          actionButton("cancel","Cancel"),
          actionButton("ok", "OK")
        )
      )
    }

    # Show modal when button is clicked.
    observeEvent(input$show, {
      showModal(dataModal())
    })

    observeEvent(input$data, {
      d.new <<- Sys.time()
    })

    observeEvent(input$load, {
      l.new <<- Sys.time()
    })

    observeEvent(input$cancel, {
      removeModal()
      # print("cancel")
      err <<- TRUE
      shinyjs::disable("loadcheck") ; shinyjs::disable("load") ; shinyjs::disable("cell") ; shinyjs::disable("yaxis") ; shinyjs::disable("xaxis")
      shinyjs::disable("hline"); shinyjs::disable("yMax") ; shinyjs::disable("yMin") ; shinyjs::disable("xMax") ; shinyjs::disable("xMin")
      shinyjs::disable("Lin.trend") ; shinyjs::disable("trans1_high") ; shinyjs::disable("trans1") ; shinyjs::disable("trans2_high"); shinyjs::disable("trans2_low"); shinyjs::disable('trans3')
      shinyjs::disable("toggle") ; shinyjs::disable("download")
    })

    # When OK button is pressed, attempt to load the data from the inputted file. If successful,
    # remove the modal. If not show another modal, but this time with a failure
    # message.
    observeEvent(input$ok, {
      flag = TRUE

      # Check whether there is a new input for the data.
      # If there isn't show a new modal with error.
      if(d.old == d.new){
        showModal(dataModal(failed = TRUE, mg = "Please select calcium data."))
        flag = FALSE
        err <<- TRUE
      }
      d.old <<- d.new

      #############
      ## Section for checking data.
      #############
      # Check that data object exists and if load is entered that it matches the data.
      if (!is.null(input$data) && flag){
        # Create the data frame that contains
        vals$data <- read.csv((input$data)$datapath)
        okay <- check(vals$data)

        # First case if everything is okay and no load then happy.
        if(!okay$flag && (is.null(input$load) || (l.new == l.old)) ){
          # print("Happy")
          # Enable all the inputs.
          shinyjs::enable("loadcheck") ; shinyjs::enable("load") ; shinyjs::enable("cell") ; shinyjs::enable("yaxis") ; shinyjs::enable("xaxis")
          shinyjs::enable("hline"); shinyjs::enable("yMax") ; shinyjs::enable("yMin") ; shinyjs::enable("xMax") ; shinyjs::enable("xMin")
          shinyjs::enable("Lin.trend") ; shinyjs::enable("trans1_high") ; shinyjs::enable("trans1") ; shinyjs::enable("trans2_high"); shinyjs::enable("trans2_low"); shinyjs::enable('trans3')
          shinyjs::enable("toggle") ; shinyjs::enable("download"); shinyjs::enable("close") ; shinyjs::enable("doclose")

          # Work out how many rows/columns we need for store/df.
          no.cell <- ncol(vals$data)-1
          spike.names <- colnames(vals$data)[-1]

          # Change store to a data frame of the correct size.
          A<- matrix(NA, nrow = no.cell, ncol = 8)
          colnames(A) <- c("threshold", "T1_high","T1_low","T2_high", "T2_low","rm.trend", "T3","close")
          row.names(A) <- spike.names
          store <<- data.frame(A)

          # Change df to allow for enough columns for each of the spikes.
          A <- matrix(NA, nrow = 4, ncol = no.cell)
          colnames(A)<- spike.names
          df <<- data.frame(A)

          # print(store)

          # Update cell input so max is the number of cells in the data.
          updateNumericInput(session, "cell",
                             value = 1, min = 1, max = no.cell, step = 1)
          updateTextInput(session, "hline", value = store[input$cell,1] )
          updateNumericInput(session, "trans1_high",value = store[input$cell,2])
          # updateNumericInput(session, "trans1_low",value = store[input$cell,3])
          updateNumericInput(session, "trans2_high",value = store[input$cell,4])
          updateNumericInput(session, "trans2_low",value = store[input$cell,5])
          updateCheckboxInput(session, "Lin.trend", value = store[input$cell,6])
          updateNumericInput(session, "trans3",value = store[input$cell,7])
          updateNumericInput(session, "close",value = store[input$cell,8])

          if(!is.na(input$trans1_high) || !is.na(input$trans2_high) || !is.na(input$trans2_low) )
            updateCheckboxInput(session, "trans1", value = TRUE)

          err <<- FALSE
          removeModal()
          err <<- FALSE
        }

        # Next case is something is wrong with the inputted data.
        if(okay$flag){
          showModal(dataModal(failed = TRUE, mg = okay$message))
          flag = FALSE
          err<<- TRUE
        }
      }

      #############
      ## Section for checking load.
      #############
      # Next consider case where we also have to load transeient data.
      if(!is.null(input$load) && flag && !(l.new == l.old) ){
        # print("load")
        # Load the file, hopefully containing the saved store.
        store.out <- read.csv((input$load)$datapath, row.names=1,stringsAsFactors = FALSE)

        # Next check whether store.out has the right dimensions.
        if((dim(vals$data)[2]-1) != dim(store.out)[1]){
          showModal(dataModal(failed = TRUE, mg = "load does not correspond to the inputted data!"))
          flag = FALSE
          err <<- TRUE

        }

        if(flag && ncol(store.out) != 8){
          showModal(dataModal(failed = TRUE, mg = "load has incorrect number of columns!"))
          flag = FALSE
          err <<- TRUE

        }


        if(flag){
          # print("Into load good section!")
          # print(head(store.out))
          # If not everything is good and update inputs with the save and enable all boxes.
          colnames(store.out) <- c("threshold", "T1_high","T1_low","T2_high", "T2_low","rm.trend","T3")
          store.out <- data.frame(store.out,stringsAsFactors = FALSE)
          store <<- store.out

          # Update all the current transients according to the new store.
          updateTextInput(session, "hline", value = store[input$cell,1] )
          updateNumericInput(session, "trans1_high",value = store[input$cell,2])
          # updateNumericInput(session, "trans1_low",value = store[input$cell,3])
          updateNumericInput(session, "trans2_high",value = store[input$cell,4])
          updateNumericInput(session, "trans2_low",value = store[input$cell,5])
          updateCheckboxInput(session, "Lin.trend", value = store[input$cell,6])
          updateNumericInput(session, "trans3",value = store[input$cell,7])
          updateNumericInput(session, "close",value = store[input$cell,8])


          if(!is.na(input$trans1_high) || !is.na(input$trans2_high) || !is.na(input$trans2_low) )
            updateCheckboxInput(session, "trans1", value = TRUE)

          # Enable all the inputs.
          shinyjs::enable("loadcheck") ; shinyjs::enable("load") ; shinyjs::enable("cell") ; shinyjs::enable("yaxis") ; shinyjs::enable("xaxis")
          shinyjs::enable("hline"); shinyjs::enable("yMax") ; shinyjs::enable("yMin") ; shinyjs::enable("xMax") ; shinyjs::enable("xMin")
          shinyjs::enable("Lin.trend") ; shinyjs::enable("trans1_high") ; shinyjs::enable("trans1") ; shinyjs::enable("trans2_high"); shinyjs::enable("trans2_low"); shinyjs::enable('trans3')
          shinyjs::enable("toggle") ; shinyjs::enable("download") ; shinyjs::enable("close"); shinyjs::enable("doclose")

          # Work out how many rows/columns we need for store/df.
          no.cell <- ncol(vals$data)-1

          # Change df to allow for enough columns for each of the spikes.
          A <- matrix(NA, nrow = 4, ncol = no.cell)
          colnames(A)<- 1:no.cell
          df <<- data.frame(A)

          # Update cell input so max is the number of cells in the data.
          updateNumericInput(session, "cell",
                             value = 1, min = 1, max = no.cell, step = 1)


          err <<- FALSE
          removeModal()
        }
      }

      l.new <<- l.old
    })

    # Function that checks the inputted data is okay.
    check <- function(data){
      flag <- FALSE
      message <- NULL

      # Initially check that the experiment time is complete
      if(any(is.na(data[,1]))){
        flag <- TRUE
        value <- which(is.na(data[,1]))
        message <- c(message,"Error! The time indexing of your data has a missing value!  <br>")
      }

      # Check that the data doesn't contain any missing values.
      for(i in 2:ncol(data)){
        A <- which(is.na(suppressWarnings(as.numeric(as.character(data[,i])))))
        B <- which(is.na(data[,i]))
        if(!setequal(A,B)){
          flag = TRUE
          # at <- colnames(data)[i]  # use this if you wanted the actual column name.
          message <- c(message, paste("Cell",i, " contains a non-numeric concentration value! <br/>"))
        }
      }

      message <- paste0(message, collapse = "")
      return(list(flag = flag, message = message))

    }
  }

  # This observe event is used to save the values for the cell we just used and reset the argumements for the new cell.
  observeEvent(input$cell, {
    # print("Changed Cell!!")
    # print(input$cell)
    updateCheckboxInput(session, 'xaxis', value = T)
    updateCheckboxInput(session, 'yaxis', value = T)
    # print(store)
    if(!is.null(input$data)){
      updateTextInput(session, "hline", value = store[input$cell,1] )
      updateNumericInput(session, "trans1_high",value = store[input$cell,2])
      updateNumericInput(session, "trans1_low",value = store[input$cell,3])
      updateNumericInput(session, "trans2_high",value = store[input$cell,4])
      updateNumericInput(session, "trans2_low",value = store[input$cell,5])
      updateCheckboxInput(session, "Lin.trend", value = store[input$cell,6])
      updateNumericInput(session, "trans3",value = store[input$cell,7])
      updateNumericInput(session, "close",value = store[input$cell,8])

      # print("val")
      # print(input$trans1_high)
      # print(input$trans2_high)
      # print(input$trans2_low)
      if(!is.na(store[input$cell,2]) || !is.na(store[input$cell,4]) || !is.na(store[input$cell,5]) ){
        # print("go true")
        updateCheckboxInput(session, "trans1", value = TRUE)
      }
      else{
        updateCheckboxInput(session, "trans1", value = FALSE)
        # print("go false")
      }
      get.spikes()
    }

  })

  # observe toggle and swap view screen
  observeEvent(input$toggle, {
    if(toggle == TRUE){
      toggle <<- FALSE
      shinyjs::enable("loadcheck") ; shinyjs::enable("load") ; shinyjs::enable("cell") ; shinyjs::enable("yaxis") ; shinyjs::enable("xaxis")
      shinyjs::enable("hline"); shinyjs::enable("yMax") ; shinyjs::enable("yMin") ; shinyjs::enable("xMax") ; shinyjs::enable("xMin")
      shinyjs::enable("Lin.trend") ; shinyjs::enable("trans1_high") ; shinyjs::enable("trans1") ; shinyjs::enable("trans2_high"); shinyjs::enable("trans2_low"); shinyjs::enable("trans3")
   shinyjs::enable("download")
    }
    else{
      toggle <<- TRUE
      shinyjs::disable("loadcheck") ; shinyjs::disable("load") ; shinyjs::disable("cell") ; shinyjs::disable("yaxis") ; shinyjs::disable("xaxis")
      shinyjs::disable("hline"); shinyjs::disable("yMax") ; shinyjs::disable("yMin") ; shinyjs::disable("xMax") ; shinyjs::disable("xMin")
      shinyjs::disable("Lin.trend") ; shinyjs::disable("trans1_high") ; shinyjs::disable("trans1") ; shinyjs::disable("trans2_high"); shinyjs::disable("trans2_low"); shinyjs::disable("trans3")
      shinyjs::disable("download")
    }
  })

  # Function to check whether the inputted threshold text is allowable.
  threshold <- function(thres){
    # Want to check that the expression is of the correct form.
    flag <- TRUE
    mg = ""
    if(is.na(thres)){ return(list(flag = flag, mg = mg))}
    # First check it only contains numbers and ,
    letters <- c("!","@","q","w","e","r","t","y","u","i","o","p",
                 "\\[","\\]","a","s","d","f","g","h","j","k","l",";","'","z","x","c","v","b","n","m",
                 "\\`","\\/","\\~","\\@","\\£","\\$","\\%","\\^","\\&","\\*","\\(","\\)","\\_",
                 "\\+","\\-","\\=","\\{","\\}","\\:","\\|","\\\\", "\\,\\,", "\\.\\.",
                 "Q","W","E","R","T","Y","U","I","O","P","A","S","D","F","G","H","J","K","L","Z","X","C","V","B","N","N")
    for(i in letters){
      if( grepl(i, thres)){
        flag <- FALSE
        mg <- "Threshold input must only contain numbers, seperated by ',' "
        break
      }
    }
     # cat("flag after check letters = ",flag,"\n")
    # check it doesn't end with a comma
    if(endsWith(thres, ",")){flag <- FALSE ; mg = ""}

    # Check the correct number of in the list
    if(flag){
      text <- paste0("c(",thres,")",collapse = "")
      thres <- eval(parse(text = text))
      if((length(thres) %% 2) == 0) {flag <- FALSE; mg = ""}
    }
    return(list(flag = flag, mg = mg))
  }

  # function that converts the data to spikes.
  get.spikes <- eventReactive(c(input$hline,
                                input$Lin.trend,
                                input$trans2_low,
                                input$trans2_high,
                                input$trans1_high,
                                input$trans3,
                                input$ok,
                                input$close),
                              {
                                # Print to show where we are in the console, and print the current inputs.
                                # print("Updating spikes...")
                                # cat("current input vals: ", input$hline,
                                    # input$Lin.trend,
                                    # input$trans2_low,
                                    # input$trans2_high,
                                    # input$trans1_high,
                                    # input$trans3,".\n")

                                # Update the store for the inputs.
                                store[input$cell,1] <<- input$hline
                                store[input$cell,2] <<- input$trans1_high
                                store[input$cell,4] <<- input$trans2_high
                                store[input$cell,5] <<- input$trans2_low
                                store[input$cell,6] <<- input$Lin.trend
                                store[input$cell,7] <<- input$trans3
                                store[input$cell,8] <<- input$close
                                # closeness <- 1
                                closeness <- input$close
                                if(is.na(closeness) || closeness == ''){
                                  closeness = 1
                                }

                                # Check conditions whether we need to compute spikes.
                                if (is.null(input$data) || is.na(store[input$cell,1])) return(NULL)
                                if(!threshold(store[input$cell,1])$flag) return(NULL)


                                # Set the trainsients. If trans1_lower = NA need to set to 0.
                                trans1 <- store[input$cell,2]
                                if(trans1 == "" || is.na(trans1)){
                                  trans1 <- 0
                                }
                                transients <- c(trans1, min(store[input$cell,5],store[input$cell,4]) , max(store[input$cell,5],store[input$cell,4]) )
                                if(any(!is.na(transients))){
                                  transients <- transients[!is.na(transients)]
                                }
                                else{ transients <- 0}

                                # Check the spike.occur is reasonable.
                                text <- paste0("c(",store[input$cell,1],")",collapse = "")
                                spike.occur <- eval(parse(text = text))
                                leng <- length(spike.occur)
                                if (leng > 1){
                                  val.change <- spike.occur[(leng/2+1.5):leng]
                                  spike.occur <-spike.occur[-((leng/2+1.5):leng)]
                                }
                                else{
                                  val.change <- NULL
                                }


                                # Get the time.step used in the experiment.
                                # time.step <- vals$data[4,1] - vals$data[3,1]
                                time.step <- vals$data[,1]

                                ##################################################################
                                # NEED CLEVER WAY OF USING TRANS3 TO CHANGE THE DATA END TIME.
                                ##################################################################
                                # Find the element number corresponding to the maximum value in time.step which is < trans3.
                                if(!is.na(input$trans3)){
                                  max.time <- max(which(time.step < (input$trans3 + 1e-10)))
                                }
                                else{max.time <- length(time.step)}


                                # Get the raw spike positions. (For rugs on time series plot.)
                                spike.pos <- SimSpikes::show_spikes(vals$data[1:max.time,input$cell+1],transients = transients, time.step = time.step[1:max.time], spike.occur = spike.occur,val.change =val.change,
                                                         closeness = closeness, return.ans = TRUE)

                                spikes <- SimSpikes::data_to_spike(vals$data[1:max.time,input$cell+1],transients = transients, time.step = time.step[1:max.time], spike.occur = spike.occur , val.change =val.change, closeness = closeness, rm.trend = FALSE)

                                # Do we need to remove linear trend?
                                if(store[input$cell,6] == TRUE){
                                  # Split into cases where we need two transient periods or only 1.

                                  if(length(transients) == 3){
                                    # Second transient period is active. Print to console.
                                    # print("entered trans > 1")

                                    # Next need to check whether the second transient period overlaps the end time.
                                    if(transients[3] > time.step[length(time.step)]){
                                      # It does so only one trend needs removing.
                                      # print('2nd trans period overlaps')
                                      # print(spikes)
                                      spikes <- SimSpikes::ISI_to_spike(remove_trend(spike_to_ISI(spikes, adjust = FALSE), rm.end.time = TRUE))
                                    }
                                    else{
                                      # It doesn't overlap so need to remove both trends.
                                      # Find how many ISIs occur before the transients.
                                      no.low <- length(spike.pos[spike.pos < store[input$cell,5]])

                                      # Split spike sequence into two
                                      spike.low <- spikes[1:no.low]
                                      spike.high <- spikes[(no.low+1):(length(spikes))]

                                      # print(spikes)
                                      # print(spike.low)
                                      # print(spike.high)

                                      # Convert spikes to ISI and remove trends
                                      ISI.low <- SimSpikes::remove_trend(spike_to_ISI(spike.low), rm.end.time = FALSE)
                                      ISI.high <- SimSpikes::remove_trend(spike_to_ISI(spike.high, adjust = TRUE), rm.end.time = TRUE)

                                      # if(length(spike.high) == 2) ISI.high = NULL
                                      # print("ISI")
                                      # print(ISI.low)
                                      # print(ISI.high)
                                      ## check that there is no linear trend in part 2.
                                      # spike.no <- 1:length(ISI.high)
                                      # fit  <- lm(spike.no~ISI.high)
                                      # print(summary(fit))

                                      # Join the ISIs together.
                                      spikes <- SimSpikes::ISI_to_spike(c(ISI.low, ISI.high))
                                    }
                                  }
                                  else{
                                    # print("entered trans < 1")
                                    # Need to only remove one transient period.
                                    spikes <- SimSpikes::data_to_spike(vals$data[1:max.time,input$cell+1],transients = transients, time.step = time.step[1:max.time], spike.occur = spike.occur , val.change =val.change, closeness = closeness,
                                                            rm.trend = TRUE)
                                  }
                                }

                                spikes.out <- spikes
                                # Next we want to save the spikes into df.

                                # Add NA's to spikes if nrow df larger then length of spikes
                                if(length(spikes) < nrow(df) ){
                                  na.s <- rep(NA,(nrow(df)-length(spikes)))
                                  spikes <- c(spikes, na.s)
                                }

                                # Add NA to end of data frame if length spikes is larger.
                                if(length(spikes) > nrow(df) ){
                                  df[(nrow(df)+1):(length(spikes)),] <<- NA
                                }
                                df[,input$cell] <<- spikes



                                return(list(spikes = spikes.out, spike.pos = spike.pos))
                              })

  # Use brush to update plot axes.
  observeEvent(input$plot_brush,{
    # print('Brushed on plot')
    # print(paste0('xmin = ',input$plot_brush$xmin,' and xmax = ',input$plot_brush$xmax))
    # Update the xmin and xmax values and uncheck the default values
    updateNumericInput(session,'xMin',value = input$plot_brush$xmin)
    updateNumericInput(session,'xMax',value = input$plot_brush$xmax)
    updateCheckboxInput(session, 'xaxis', value = F)

    #  Update ymin and ymuax but DONT uncheck the box.
    tau <- vals$data[,1]
    brushed <-input$plot_brush
    dat <- vals$data[,input$cell+1]

    low <- max(which(tau < input$plot_brush$xmax))
    high <- min(which(tau > input$plot_brush$xmin))
    ymin <- min(vals$data[low:high,input$cell+1])
    ymax <- max(vals$data[low:high,input$cell+1])
    # Check that the values don't go out of bounds.
    if(is.na(ymin) || is.na(ymax)){
      ymin <- min(dat,na.rm = T)
      ymax <- max(dat, na.rm = T)
    }
    updateNumericInput(session,'yMin',value = ymin)
    updateNumericInput(session,'yMax',value = ymax)

    session$resetBrush("plot_brush")
  })

  # Create time series plot.
  output$lineplot <- renderPlot({
    data <- vals$data
    input$toggle ;

    # print("Do lineplot")
    # print(head(data))
    if (is.null(input$data)|| is.na(input$cell) || !input$ok || err || toggle )
      return(NULL)

    if ( !(input$cell %in% c(1:(ncol(data)-1))))
      return(NULL)

    if(input$yMax <= input$yMin || input$xMax <= input$xMin){
      return(NULL)
    }

    get.spikes()
    # decide the y limits.
    if(input$yaxis == FALSE){
      ymin <- input$yMin ; ymax <- input$yMax
    }
    else{
      cur <- data[,input$cell+1]
      ymin <- min(cur[!is.na(cur)]) ; ymax <- max(cur[!is.na(cur)])
    }

    # Decide the x limits.
    if(input$xaxis == FALSE){
      xmin <- input$xMin ; xmax <- input$xMax
    }
    else{
      cur <- data[,1]
      xmin <- min(cur[!is.na(cur)]) ; xmax <- max(cur[!is.na(cur)])
    }

    # Create plot.
    color = "#434343"
    par(mar = c(4, 4, 1, 1))
    plot(x =  data[,1], y =  data[,input$cell+1], type = "l",
         xlab = "Time", ylab = "Calcium conc", xlim = c(xmin,xmax), ylim = c(ymin,ymax), col = color, fg = color, col.lab = color, col.axis = color)

    # Do we want to display transient periods.
    if(!is.null(input$trans1_high)){
      polygon(c(0,input$trans1_high,input$trans1_high,0), c(ymin,ymin,ymax,ymax), border = "red", col = rgb(255,0,0,max = 255, alpha = 50))
    }
    # Only allow to show second region if it happens after the first interval and low < high.
    # if(!is.na(input$trans2_low) && !is.na(input$trans2_high)){
    #   if(is.na(input$trans1_high) && input$trans2_low < input$trans2_high){
    #     print('2nd trans: No 1 entered')
    #     polygon(c(input$trans2_low,input$trans2_high,input$trans2_high,input$trans2_low), c(ymin,ymin,ymax,ymax), border = "red", col = rgb(255,0,0,max = 255, alpha = 50))
    #   }
    #   if(!is.na(input$trans1_high) && (input$trans1_high < input$trans2_low) && (input$trans2_low < input$trans2_high) ){
    #     print('2nd trans: No 2 entered')
    #     polygon(c(input$trans2_low,input$trans2_high,input$trans2_high,input$trans2_low), c(ymin,ymin,ymax,ymax), border = "red", col = rgb(255,0,0,max = 255, alpha = 50))
    #   }
    # }
    if(!is.na(input$trans2_low) && !is.na(input$trans2_high)){
      if(is.na(input$trans1_high)){
        polygon(c(input$trans2_low,input$trans2_high,input$trans2_high,input$trans2_low), c(ymin,ymin,ymax,ymax), border = "red", col = rgb(255,0,0,max = 255, alpha = 50))
      }
      if(!is.na(input$trans1_high) && (input$trans1_high < input$trans2_low) && (input$trans1_high < input$trans2_high) ){
        polygon(c(input$trans2_low,input$trans2_high,input$trans2_high,input$trans2_low), c(ymin,ymin,ymax,ymax), border = "red", col = rgb(255,0,0,max = 255, alpha = 50))
      }
    }
    if(!is.na(input$trans3)){
      polygon(c(input$trans3,max(vals$data[,1]),max(vals$data[,1]),input$trans3), c(ymin,ymin,ymax,ymax), border = "red", col = rgb(255,0,0,max = 255, alpha = 50))
    }


    if(threshold(store[input$cell,1])$flag){
      text <- paste0("c(",store[input$cell,1],")",collapse = "")
      thres <- eval(parse(text = text))
      no.lines <- (length(thres)+1)/2

      # Now get the values which we use for the change on the x-axis.
      if(no.lines == 1){x.axis <- c(0,max(data[,1]))}
      else{
        x.axis <- c(0, thres[(length(thres)-no.lines+2):length(thres)], max(data[,1]))
      }

      # Get the y.axis values.
      y.axis <- thres[1:no.lines]
      for(i in 1:no.lines){
        lines(c(x.axis[i],x.axis[i+1]),c(y.axis[i],y.axis[i]), col = "red", lwd = 2)
      }
    }

    # Add rug off the spikes.
    if(!(is.na(store[input$cell,1]) || store[input$cell,1] == "" )){
      # Only add if transient regions are reasonable.
      if((is.na(input$trans2_high) || is.na(input$trans2_low)) ||
         (!is.na(input$trans1_high) && (input$trans1_high < input$trans2_low) && (input$trans1_high < input$trans2_high))){
        rug(get.spikes()$spike.pos,ticksize = 0.05, col = "red", lwd = 2)
      }
      if(  (!is.na(input$trans2_high) || !is.na(input$trans2_low))  && (input$trans2_high > input$trans2_low)  ){
        rug(get.spikes()$spike.pos,ticksize = 0.05, col = "red", lwd = 2)
      }
  }

  })

  # Create ISI plot.
  output$ISIplot <- renderPlot({
    # First need to call the variables that the ISI plot depends on for update.
    input$hline
    input$trans1_low
    input$trans2_low
    input$trans2_high
    input$Lin.trend

    # Print to inform console we're in the Output ISI part of code.
    # print("Output ISIplot")

  ## ** List of conditions that forces the output to not appear **
    # making sure we have inputted date, taken a cell value, checking the inputs are allowable and toggle is in correct value.
    if (is.null(input$data)|| is.na(input$cell) || !input$ok || err || toggle ) {return(NULL)}
    # Checking we have inputted a threshold value.
    if (is.na(store[input$cell,1]) || store[input$cell,1] == ""){return(NULL)}
    # Checking the transient periods don't overlap and the second one has the correct order.
    if (!is.na(input$trans2_high) && !is.na(input$trans2_low)){
      # print("Why we in here!")
      # print(input$trans2_high)
      # print(input$trans2_low)
      if((!is.na(input$trans1_high) && ((input$trans1_high >= input$trans2_low) || (input$trans1_high >= input$trans2_high)) )){
        return(NULL)
      }
    }
    # Checking we have some ISIs to create the plot for.
    spikes <- get.spikes()$spikes
    if(length(spikes) <1){
      return(NULL)
    }
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    # Convert the spikes to ISI and create the order they appear in.
    ISI <- spike_to_ISI(spikes)
    ISI.no <- 1:(length(ISI)-1)

    # Check values if neccessary.
    # print(ISI)
    # print(ISI.no)

    # Create the plot of the ISI no vs ISI time.
    plot(ISI.no,ISI[ISI.no], xlab = HTML("ISI Number"), ylab = HTML("Time"), pch = 20)

  ## **Now we want to add the line of best fit through the points, However need to split into cases to allow for fit before and after a transient period.**
    # First we check if we have a second transient period.(Not at the start of recording)
    if(!is.na(store[input$cell,4]) && !is.na(store[input$cell,5])){
      # TRUE we have a transient period not at the start.

      # Next need to check we have spikes in the second interval
      # print('HOBGOBLIN')
      # print(spikes)
      # print(spikes[length(spikes)])
      # print(store[input$cell,4])
      if(store[input$cell,4] > spikes[length(spikes[!is.na(spikes)])]){
        xaxis <- 1:(length(ISI)-1)
        fit.before <- lm(ISI[xaxis] ~ xaxis )
        Y <- predict(fit.before, newdata=data.frame(x=xaxis))
        lines(xaxis,Y)
      }
      else{

      # Plot line between transients
      before <- sum(get.spikes()$spike.pos < store[input$cell,4])
      lines(c(before+0.5,before+0.5), c(0,max(ISI)+200), col = "black", lwd = 2)
      text(x = before , y = max(ISI)-10, labels = "1")
      text(x = before + 1, y = max(ISI)-10, labels = "2")

      # Add fitted linear lines.
      # In the first section (1)
      if(before > 2){
        xaxis <- 1:before
        fit.before <- lm(ISI[xaxis] ~ xaxis )
        Y <- predict(fit.before, newdata=data.frame(x=xaxis))
        lines(xaxis,Y)
      }

      # In the second section (2)
      if(length(ISI-before > 2)){
        xaxis <- (before+1):(length(ISI)-1)
        fit.after <- lm(ISI[xaxis] ~ xaxis)
        Y <- predict(fit.after, newdata=data.frame(x=xaxis))
        lines(xaxis,Y)
      }



      }
    }
    else{
      xaxis <- 1:(length(ISI)-1)
      fit.before <- lm(ISI[xaxis] ~ xaxis )
      Y <- predict(fit.before, newdata=data.frame(x=xaxis))
      lines(xaxis,Y)
    }



  })

  # This output informs the viewer what cell number they are looking at.
  output$desc <- renderUI({
    input$cancel
    input$ok
    input$toggle

    if (is.null(input$data) || err || toggle){return(HTML("<p> <font size = '+3'> Welcome to View Data </font> </p>
                                                          <p> This app is for viewing time series data of calcium concentration in a cell, and thresholding the time series data to get spike sequences. </p>

                                                          <p> <font size = '+1'> Using the app </font> </p>

                                                          <p> On the sidebar click 'INPUT DATA', this opens a pop up window where you can select the time series data. This must be of type .csv and NOT an excel spreadsheet, where the first row is column names. You can also choose to load a table with the thresholding details for each cell  (threshold value, transients, remove linear trends), saved as a .csv file. Once loaded you can use 'cell number' to scroll through the time series data. Naturally, default axes are originally chosen, however if you wish to zoom in on a particular part of the plot you can do this by unchecking the default axes and choosing your own.       </p>

                                                          <p> <font size = '+1'> Thresholding </font> </p>
                                                          <p> At the bottom of the sidebar there is a 'Thresholding' section. Enter the value you wish to threshold in the 'Threshold at?' box. If you want to threshold at different values then entering in the form c<sub>1 </sub>, c<sub>2 </sub>, ... , c<sub>k </sub>, t<sub>1 </sub>, ... ,t<sub>k-1 </sub> will threshold at c<sub>1 </sub> between [0, t<sub>1 </sub>] and c<sub>2 </sub> between [t<sub>1 </sub> t<sub>2 </sub>] and so on. This will show as a red horizontal line on the time series plot. Once you have chosen the value to threshold at, a plot showing ISIs, along with a linear fit to show any trends. Also, the spikes will show as rugs under the time series plot and shown at the bottom of the screen
                                                          </p><p>
                                                          You have the option to remove transient periods, to do so check the 'Transient period' box. Then two transient areas will pop up in the sidebar one starting at time 0 and a later transient period. Transient areas will appear as red boxes on the plot. For transient period 1, this works by beginning time from the first spike outside of the transient period. Transient area 2 works by removing all spikes in the red area and glueing the last spike before and first spike after the transient area as one spike time.
                                                          </p><p>
                                                          If you check remove linear trends, this calcuates any linear trend in the ISI of the spike times then removes the trend such that you maintain the experiment time.
                                                          In the ISI plot you can see fitted line/s that show any linear trends in the data. If you check 'remove linear trends' these trends shall be removed (whilst maintaining absolute experiment time). If you have cut a transient section in the middle of the data then linear trends are shown for both sides of where the cut took place.
                                                          </p>


                                                          <p> <font size = '+1'> Downloading the results </font> </p>

                                                          <p>You can download the results by clicking the download button, where you will be prompted to choose a location to download a .zip file. The .zip contains: the spikes are saved as spikes.csv, the thresholding values saved in store.csv and details.txt that tells you the filename of the data used. <p>
                                                          <p> If you experience any issues you can email: pmxjp8@nottingham.ac.uk </p>"))}

    #  Create message that tells them of the file the results are stored in.
    saved <- paste0("<p> Click the download button if you wish to save your spikes!</p>")

    if(is.na(input$cell)){ return(HTML(paste(saved,"Please enter a cell number.")))}
    # if(!(input$cell %in% 1:(ncol(data())-1))){ return(HTML(paste(saved,"Please enter a valid cell number, between 1 and", ncol(data())-1),"."))}
    # if(check()$flag){return(HTML(paste0(check()$message, collapse = "<br/>")))}

    HTML(paste(saved,"You have selected cell", input$cell, ", maximum cell number is", ncol(vals$data)-1,"."))
    })

  # This output is used to show the spikes they get from the real data.
  output$Spikes <- renderUI({
    input$toggle ; input$hline ; input$trans1_high ; input$trans2_high;
    input$trans2_low ; input$Lin.trend

    if (is.null(input$data) || err || toggle){return(HTML(""))}
    if (!is.na(input$trans1_high) && (!is.na(input$trans2_high) && !is.na(input$trans2_low))){
      if ((input$trans1_high >= input$trans2_low) || (input$trans1_high >= input$trans2_high)){
        return(HTML("<font color = 'red'> Transient regions must not overlap! </font>"))
      }
    }
    # print("output Spikes")

    # Check whether the inputted threshold is allowable.
    okay <- threshold(input$hline)
    if(!okay$flag){
      mg <- paste0("<hr>", okay$mg,collapse = "")
      return(HTML(mg))
    }

    # spikes <- df[,input$cell]
    # spikes <- spikes[!is.na(spikes)]
    spikes <- get.spikes()$spikes
    message <- ""
    if(length(spikes) == 0){
      message <- c(message,"There are no spikes.")
    }
    if(length(spikes) == 1){
      message <- c(message,paste("The spike time is: ",spikes[1],"."))
    }
    if(length(spikes) == 2){

      message <- c(message,(paste("The spikes times are: ",spikes[1],",",spikes[2],".")))
    }
    if(length(spikes)>2){
      message <- c(message,"The spikes times are: <font size ='-1'> ")
      for(i in 1:(length(spikes)-1)){message <- c(message, paste(round(spikes[i],digits = 2),", "))}
      message <- c(message, paste(spikes[length(spikes)],". </font>"))
    }
    all <- paste(message, collapse= "")

    all <- c("<hr>", all)
    # print("--------------------")
    HTML(paste(all, collapse = "<br/>"))
  })

  # Download the current files.
  {
    output$download <- downloadHandler(
      filename = function(){
        text <- paste0(c("ViewData"))
        paste0(text,".zip")

      },
      content = function(file){
        #go to a temp dir to avoid permission issues
        owd <- setwd(tempdir())
        on.exit(setwd(owd))
        files <- NULL;

        # Next create the files to put into the zipped folder.

        fileConn<-file(paste0("Details.txt"))
        writeLines(c(paste0("The file the data comes from is ",(input$data)$name)), fileConn)
        close(fileConn)

        # Save the store.
        fileName.store <- paste("store.csv",sep = "")
        write.csv(store, file = fileName.store)

        # Save the spikes in csv file.
        fileName.spikes <- paste("spikes.csv",sep = "")
        write.csv(df, file = fileName.spikes)


        files <- c("Details.txt", fileName.store, fileName.spikes)
        # #loop through the sheets
        # for (i in 1:input$sheet){
        #   #write each sheet to a csv file, save the name
        #   fileName <- paste(input$text,"_0",i,".csv",sep = "")
        #   write.table(data()$wb[i],fileName,sep = ';', row.names = F, col.names = T)
        #   files <- c(fileName,files)
        # }

        #create the zip file
        zip(file,files)
      }
    )
  }


  }
