# server.r
# Mark Hagemann
# 7/16/2014
# For shiny implementation of load estimation for DCR tributaries

require(ggplot2)
require(plyr)
require(reshape2)
library(car)

# Load
source("functions/metaFunctions.r")
source("functions/functions2.r")
source("functions/functions_comparing.R")
load("data/dailyData.rda")
load("data/sampleDaysData.rda")
load("data/selected_models.rda")

shinyServer(function(input, output, session) {
  
  ####
  #### Reactively generate objects ####
  #### 
  
    ## Select particulars of analysis  ####
  
  inputData = reactive({
    makeModelData(caldata = gdat, preddata = preds.c, location = input$station,
                  var = input$variable, startDate = input$dateRange[1],
                  endDate = input$dateRange[2])
  })
  
  selectedModels = reactive({
    modname = paste0(input$station, input$variable)
    allBestModels[[modname]]
  })

    
  #### Results #### 
  
  resultsTS = reactive({
    dat = getPredictions(lm.model=selectedModels()$bestLm,  
                   gam.model=selectedModels()$bestGam,
                   obsdata = inputData()$caldata, 
                   preddata = inputData()$preddata,
                   interpdata = inputData()$interpdata,
                   bound.conc = input$boundconc)
    dat
  })
  
  results.m = reactive({
    r1 = melt(resultsTS(), id.vars = c("Date", "sampleType"), na.rm = T)
    r1$variable = gsub("se.", "se_", r1$variable)
    r1$qty = sapply(strsplit(as.character(r1$variable), ".", fixed = T), "[", 1) 
    r1$method = sapply(strsplit(as.character(r1$variable), ".", fixed = T), 
                       function(x) paste(x[-1], collapse = ".")) 
    
    r2 = dcast(r1, Date + method + sampleType ~ qty)
    
    r2 = within(r1, {
      isFlow = ifelse(variable == "CFS", "Flow", "WQ")
      is.se = ifelse((substr(variable, 1, 3) == "se_"), "SE", "value")
      isPoint = grepl("meas", variable)
      qty = gsub("se_", "", qty)      
    })
    
    r2 = dcast(r2, Date + sampleType + method + qty + isPoint + isFlow ~ is.se)
    
  })

  resultsForGG = reactive({
    require(plyr)
    r1 = results.m()[results.m()$method %in% c("", input$yplot1),]
    r2 = r1[r1$qty %in% c("CFS", input$toPlot),]
  })
  
  tsPlot1 = reactive({
    
    hasMeas = "meas" %in% input$yplot1
    ynames1 = input$yplot1
    if(hasMeas) ynames1 = ynames1[-which(ynames1 == "meas")]
    ynames = paste(input$toPlot, ynames1, sep = ".")
    xdata = resultsTS()$Date
    ydata = as.data.frame(resultsTS()[,ynames])
    names(ydata) = ynames
    maxy = max(ydata, na.rm = T)
    
    plot(xdata, resultsTS()[, paste0(input$toPlot, ".meas")], 
         type = ifelse(hasMeas, "p", "n"), 
         ylim = c(0, maxy), ylab = input$toPlot)
    numind = 1
    for(ycol in ynames){
      lines(xdata, ydata[, ycol], col = numind, lty = numind)
      numind = numind + 1
    }
    
    legend("topleft", legend = ynames1, col = 1:length(ynames), lty = 1:length(ynames))
  })
  
  tsPlot1GG = reactive({
    
    #test = resultsForGG()
    #browser()
    
    g1 = ggplot()
    
    if(sum(!resultsForGG()$isPoint) > 0){
      g1 = g1 + geom_line(data = subset(resultsForGG(), !isPoint), 
                          aes_string(x = "Date", y = "value", 
                                     colour = "method"))  
    }
    
    if(sum(resultsForGG()$isPoint) > 0){
      g1 = g1 + geom_point(data = subset(resultsForGG(), isPoint), 
                           aes_string(x = "Date", y = "value"))
    }
    g1 + facet_wrap(~ isFlow, nrow = 2, scales = "free_y")
  })
  
  loadTable = reactive({
    load_table(resultsTS(), seMethod = input$seMethod, use = input$covmatQuantity,
               totvar = input$addResidVar)
  })
  
### Model Inspection, comparison ####

  modelForms = reactive({
    modelFormsTable(lm.model = selectedModels()$bestLm, 
                    gam.model = selectedModels()$bestGam)
  })
  
  
  lmTermList = reactive({all.vars(terms(selectedModels()$bestLm))[-1]})
  gamTermList = reactive({all.vars(terms(selectedModels()$bestGam))[-1]})
  
  observe({
    updateRadioButtons(session, "termVarLm", label = "LM terms", 
                       choices = lmTermList())
    updateRadioButtons(session, "termVarGam", label = "GAM terms", 
                       choices = gamTermList())
    updateRadioButtons(session, "residXvarLm", label = "LM X variable",
                       choices = lmTermList())
    updateRadioButtons(session, "residXvarGam", label = "GAM X variable",
                       choices = gamTermList())
  })
  
  lmTermPlotData = reactive({
    pdat = termplot(selectedModels()$bestLm, se = T, plot = F) # a list
    resids = residuals(selectedModels()$bestLm, "response")
    presids = residuals(selectedModels()$bestLm, "partial") # a data frame
    rdat = within(selectedModels()$bestLm$model, {
      residuals = resids
      partial.resids = presids})
    list(estimates = pdat, residuals = rdat)
  })
  lmResidualData = reactive({lmTermPlotData()$residuals})
  
  gamTermPlotData = reactive({
    pdat.p = termplot(selectedModels()$bestGam, se = T, plot = F) # a list
    pdat.s = gamTermplot(selectedModels()$bestGam)
    pdat = c(pdat.p, pdat.s)
    
    resids = residuals.glm(selectedModels()$bestGam, "response")
    presids = residuals.glm(selectedModels()$bestGam, "partial") # a data frame
    
    colnames(presids)[grepl("(^s\\()(.*)(\\)$)", colnames(presids), perl = T)] = 
      unlist(regmatches(colnames(presids), gregexpr("(?<=s\\().*?(?=\\))", 
                                        colnames(presids), perl=T),
                        invert = F))
    rdat = within(selectedModels()$bestGam$model, {
      partial.resids = presids
      residuals = resids})
    #browser()
    list(estimates = pdat, residuals = rdat)
  })
  gamResidualData = reactive({
    
    gamTermPlotData()$residuals})

  lmTermPlots = reactive({
    termplotGG(estdat = lmTermPlotData()$estimates[[input$termVarLm]], 
               residdat = lmTermPlotData()$residuals, 
               varname = input$termVarLm)
  })
  
  gamTermPlots = reactive({
    termplotGG(estdat = gamTermPlotData()$estimates[[input$termVarGam]],
               residdat = gamTermPlotData()$residuals,
               varname = input$termVarGam)    
  })
  
  lmResidPlots = reactive({
    residplotGG(residdat = lmResidualData(), 
                xvar = input$residXvarLm, yvar = "residuals")
  })
  
  gamResidPlots = reactive({
    residplotGG(residdat = gamResidualData(),
                xvar = input$residXvarGam, yvar = "residuals")
  })
  
  
  lmResid = reactive({
    residuals(selectedModels()$bestLm)
  })
  gamResid = reactive({
    residuals(selectedModels()$bestGam)
  })
  
  ####
  #### Make Output Objects ####
  ####
  
  output$obsPlot = renderPlot({
    plot(value ~ Date, inputData()$caldata, col = sampleType)})
  
  output$obsTable = renderTable({
    within(inputData()$caldata, {
      Date = format(Date)
      })[c("Date", "value", "sampleType")]
  })
  
  output$TSplot1 = renderPlot({
  #     tsPlot1()
    tsPlot1GG()
    })
  
  output$loadTable = renderTable({loadTable()})
  
  output$loadBarPlot = renderPlot({
    load_barplot(loadTable())
  })

  output$downloadData = downloadHandler(
    filename = function() {
      paste0(input$station, '_', input$variable, '_', 
      paste(format(input$dateRange, "%Y%m%d"), collapse = '-'), '.csv')
      },
    content = function(file) {
      write.csv(resultsTS(), file, row.names = F)
    })
  
  output$resultsTS = renderTable({within(resultsTS(), {Date = format(Date)})})
  
  output$fitPlotConc = renderPlot({
    plot(conc.lm ~ conc.meas, resultsTS(), col = 1, pch = as.numeric(sampleType))
    points(conc.gam ~ conc.meas, resultsTS(), col = 2, pch = as.numeric(sampleType))
    abline(0, 1)
    })
  
  output$fitPlotLoad = renderPlot({
    plot(load.lm ~ load.meas, resultsTS(), col = 1, pch = as.numeric(sampleType))
    points(load.gam ~ load.meas, resultsTS(), col = 2, pch = as.numeric(sampleType))
    abline(0, 1)
    })
  
  output$cvTable = renderTable({data.frame(a = 1:10, b = 2:11)})
  output$modelFormTable = renderTable({modelForms()})
  output$termPlotsLM = renderPlot({lmTermPlots()})
  output$termPlotsGAM = renderPlot({gamTermPlots()})

  output$qqPlotLM = renderPlot({qqPlot(lmResid())})
  output$qqPlotGAM = renderPlot({qqPlot(gamResid())})
  output$boxplotLM = renderPlot({boxplot(lmResid())})
  output$boxplotGAM = renderPlot({boxplot(gamResid())})
  output$acfLM = renderPlot({acf(lmResid())})
  output$acfGAM = renderPlot({acf(gamResid())})
  
  output$lmResidPlot = renderPlot({
    g1 = ggplot()
    g1 = g1 + geom_point(data = lmResidualData(), 
                         aes_string(x = input$residXvarLm, 
                                    y = "residuals"))
    print(g1)
  })
  output$gamResidPlot = renderPlot({
#     gamResidPlots()
    g1 = ggplot()
    g1 = g1 + geom_point(data = gamResidualData(), 
                         aes_string(x = input$residXvarGam, 
                             y = "residuals"))
    print(g1)
    })
  
#   output$test = renderTable({head(gamTermPlotData()$residuals)})
  output$test2 = renderText({cat(str(gamTermPlotData()))})
  output$test3 = renderPlot({tsPlot2()})
  })
