# functions2.r
# Revised versions of functions.r containing only functions needed 
# in reorganization of gam project



makeModelData = function(caldata, preddata, location, var, startDate, endDate, 
                         outs_hi = 0, outs_lo = 0, 
                         sset = c("allSamples", "noWetRoutine", "noWetGrab")){
  #
  # Modified from makeModelData in functions.r
  # returns a list of data.frames for fitting/predicting, including all predictors.
  # Note that some of these are dependent on start and end times, thus
  # it is important that "startDate" and "endDate" are the same as was used
  # in original calibration.
  #
  # Now also includes an interpolation data frame as a returned object.
  #
  
  # Subset data
  sset = match.arg(sset)
  exp1 = ("!(sampleType == 'monthly grab' & precip_in + lagprecip > 0.1)")
  exp2 = ("!(sampleType %in% c('monthly grab', 'falling limb grab') & precip_in + lagprecip > 0.1)")
  
  subSets = c(noWetRoutine = exp1, noWetGrab = exp2, allSamples = NA)
  sset = subSets[sset]
  
  out.dat1 = subset(caldata, Station == location & variable == var)
  if(!is.na(sset)){
    out.dat1 = subset(out.dat1, eval(parse(text = sset)))
  }
  out.dat1 = out.dat1[order(out.dat1$value),]
  out.dat1 = out.dat1[(outs_lo + 1) : (nrow(out.dat1) - outs_hi), ]
  out.dat1 = out.dat1[order(out.dat1$Date), ]
  
  # Add properly scaled time data
  Dates = seq.Date(as.Date("1998-6-6"), as.Date("2015-1-1"), by = 1)
  times = data.frame(Date = Dates, t = as.numeric(Dates))
  meanD = mean(c(Dates[1], Dates[length(Dates)]))
  meant = as.numeric(meanD)
  
  times$t = times$t - meant
  times$t2 = times$t^2
  
  out.dat1 = merge(out.dat1, times, by = "Date", all.x = T)
  out.dat2 = subset(out.dat1, Date >= startDate & Date <= endDate)
  out.dat2$lc = log(out.dat2$value)
  
  out.dat3 = merge(preddata, times, by = "Date", all = T)
  out.dat3 = subset(out.dat3, Station == location & Date >= startDate &
                      Date <= endDate)
  
  # interpolation data:
  interpdata = merge(out.dat3[c("Date", "CFS")], 
                     out.dat1[c("Date", "value", "sampleType")],
                     all = T)
  
  # return list of calibration and prediction data
  out = list(caldata = out.dat2, preddata = out.dat3, interpdata = interpdata)
}


getPvals = function(model){
  if("gam" %in% class(model)){
    sigs1 = summary(model)$p.pv
    s.tab = summary(model)$s.table
    sigs2 = s.tab[,"p-value"]
    names(sigs2) = rownames(s.tab)
    sigs = c(sigs1, sigs2)
  } else if(class(model) == "lm"){
    sigs = summary(model)$coefficients[,"Pr(>|t|)"]
  } else stop("Model must be either 'lm' or 'gam'!")
  sigs[!grepl("ntercept", names(sigs))]
}

tidyModel = function(model){
  # 
  # Performs stepwise removal of predictors until only significant ones remain
  # 
    
  mdat = model$model
  terms = all.vars(terms(model)) 
  yname = terms[1]
  
  pvals = getPvals(model)
  ninsig = sum(pvals > 0.05)
  
  tempmod = model
  
  while(ninsig > 0){
    print(pvals)
    # remove a single term
    xnames = names(pvals)
    to.rm = xnames[pvals == max(pvals)]
    print(to.rm)
    xnames = xnames[names(pvals) != to.rm]
    
    if(length(xnames) == 0){
      warning("No significant terms, use average instead.")
      break 
    }
    
    formRHS = paste(xnames, collapse = " + ")
    form1 = as.formula(paste(yname, formRHS, sep = "~"))
    
    if("gam" %in% class(model)) tempmod = gam(formula = form1, data = mdat) else 
      tempmod = lm(formula = form1, data = mdat, y = T)
    pvals = getPvals(tempmod)
    ninsig = sum(pvals > 0.05)
    
  }
  return(tempmod)
}

gamTermplot = function(x, seWithMean = F){
  # 
  # x is a gam object
  #
  if(length(x$smooth) == 0) return(NULL)
  
  outlist = list()
  
  for(i in 1:length(x$smooth)){
    tmsi = x$smooth[[i]]
    P = mgcv:::plot.mgcv.smooth(tmsi, P = NULL, data = x$model)
    first = tmsi$first.para; last = tmsi$last.para
    p <- x$coefficients[first:last]
    
    P$fit = P$fit = P$X %*% p
    
    if (seWithMean && attr(tmsi, "nCons") > 
          0) {
      if (length(x$cmX) < ncol(x$Vp)) 
        x$cmX <- c(x$cmX, rep(0, ncol(x$Vp) - length(x$cmX)))
      X1 <- matrix(x$cmX, nrow(P$X), ncol(x$Vp), 
                   byrow = TRUE)
      meanL1 <- tmsi$meanL1
      if (!is.null(meanL1)) 
        X1 <- X1/meanL1
      X1[, first:last] <- P$X
      se.fit <- sqrt(pmax(0, rowSums((X1 %*% x$Vp) * 
                                       X1)))
    }
    else se.fit <- sqrt(pmax(0, rowSums((P$X %*% 
                                           x$Vp[first:last, first:last, drop = FALSE]) * 
                                          P$X)))
    P$se.fit = se.fit
    dat.fit = with(P, data.frame(x = x, y = fit, se = se.fit))
    outlist[[i]] = dat.fit
  }
  
  names(outlist) = sapply(x$smooth, "[[", "term")
  outlist
}


termplotGG = function(estdat, residdat, varname, xlab, ylab, title){
  # estdat is a data.frame of the kind returned by termplot(plot = F)[[1]]
#   browser()
  rdat = data.frame(x1 = residdat[,varname], y1 = residdat$partial.resids[,varname])
  g1 = ggplot()
  g1 = g1 + geom_ribbon(data = estdat, aes(x = x, ymin = y - 1.96 * se, 
                                           ymax = y + 1.96 * se, fill = "gray46"))
  g1 = g1 + geom_line(data = estdat, aes(x = x, y = y))
  g1 = g1 + geom_point(data = rdat, aes(x = x1, y = y1))
  g1 = g1 + labs(xlab = varname, ylab = "response")
  g1
}

residplotGG = function(residdat, xvar, yvar){
  # for making residual scatter plots in shiny app
  # residdat is a data.frame of the type returned by model$model, with added residual columns
  rd = residdat[, c(xvar, yvar)]
  localenv = environment()
  g1 = ggplot(environment = localenv)
  
  g1 = g1 + geom_point(data = rd, 
                       aes(x = rd[,xvar], 
                           y = rd[, yvar]), environment = localenv)
  g1
}

formatP = function(pval){
  # formats p-values so that they are easy to read, put in a table
  pout = sprintf("%.3f", pval)
  pout[pval < 0.001] = "<0.001"
  pout
}

predictLoad = function(data, model, start = NULL, end = NULL, na.approx = T, 
                       log = F,
                       bound.conc = T, timecol = "Date"){
  #
  ### modified frm gamLoad function in functions.r
  #### data should have predictor variables at all desired times
  #### bound.conc prevents concentrations from exceeding minimum or maximum in record
  ####
  
  require(zoo)
  require(plyr)
  
  if(is.null(start)) start = min(data$Date)
  if(is.null(end)) end = max(data$Date)
  first = as.Date(start)
  last = as.Date(end)
  dates = data.frame(Date = as.Date(as.numeric(first) : as.numeric(last), 
                                    origin = "1970-1-1"))
  
  # daily, restricted to time period of interest
  newdat = subset(data, Date >= start & Date <= end)
  agcols = c("Date", names(newdat)[sapply(newdat, is.numeric)])
  newdat = aggregate(. ~ Date, data = newdat[agcols], mean, na.action = na.pass)
  names(newdat)[1] = "Date"
  newdat = merge(newdat, dates, all = T)
  newdat$time = as.numeric(newdat$Date)
  
  if(na.approx){
    #     browser()
    predcols = attr(terms(model), "term.labels") # predictor variables
    print(predcols)
    n1 = is.na(newdat[predcols])
    if(length(predcols) > 1){
      newdat[predcols] = apply(newdat[predcols], 2, FUN = na.approx, na.rm = F)
    } else newdat[predcols] = as.data.frame(na.approx(newdat[predcols], na.rm = F))
    
    n2 = is.na(newdat[predcols])
    numInterp = sum(n1 * n2)
    if(numInterp > 0){
      warning(paste("Interpolated", numInterp, "out of", nrow(newdat), 
                    "predictor values."))
    }
  } # interpolate missing values in predictors
  
  pred1 = predict(model, newdata = newdat, se.fit = T)
  meas = model$y
  
  if(log){ # adjust back-transform using smearing estimator
    meas = exp(model$y)
    smear = mean(meas) / mean(exp(model$fitted.values))
    pred1 = as.data.frame(predict(model, newdata = newdat, se.fit = T))
    pred1$fit = smear * exp(pred1$fit)
    pred1$se.fit = sqrt(with(pred1, fit^2 * (exp(se.fit^2) - 1)))
  }
  
  if(bound.conc){ # bound concentrations to range of observed values.
    recmin = min(meas, na.rm = T)
    recmax = max(meas, na.rm = T)
    pred1$fit[pred1$fit > recmax] = recmax
    pred1$fit[pred1$fit < recmin] = recmin
  }
  
  newdat = within(newdat, {
    model.pred = as.numeric(pred1$fit)
    se.pred = as.numeric(pred1$se.fit)
    c1 = model.pred - 1.96 * se.pred
    c2 = model.pred + 1.96 * se.pred
    model.load = model.pred * CFS * 3600 * 24 * 28.317 / 1000000 # convert to kg/day
    load.min = c1 * CFS * 3600 * 24 * 28.317 / 1000000
    load.max = c2 * CFS * 3600 * 24 * 28.317 / 1000000
    se.load = se.pred * CFS * 3600 * 24 * 28.317 / 1000000
  })
  #     "load.min", "load.max", "model.cumuload")]
  dat.out = newdat
  #   if(crossval){
  #     dat.out$cv.resid = newdat$residuals
  #   }
  
  dat.out = dat.out[c("Date", "CFS", "model.pred", "se.pred", "model.load",
                      "se.load")]
  names(dat.out) = mapvalues(names(dat.out), from = c("model.pred", "model.load"),
                             to = c("pred", "load.pred"))
  dat.out
}


interpLoad = function(data, t1, t2, flowcol = "CFS", sampleTypes = "all", 
                      flow.interp = T, valcol = "value"){
  ## data$valcol MUST contain only desired concentrations 
  ## (i.e. subset data BEFORE calling this function)
  ## data$flowcol should have near-daily resolution
  ## sampleTypes can be either "all" or a character vector of values from data$sampleType
  
  names(data) = mapvalues(names(data), from = c(flowcol, valcol),
                          to = c("CFS", "value"))
  
  concdat = data[!is.na(data$value), c("Date", "value", "sampleType")] # only concentration rows
  
  if(sampleTypes != "all"){
    concdat = subset(concdat, sampleType %in% sampleTypes)
  }
  flowdat = data[!is.na(data$CFS), c("Date", "CFS")] # only flow data
  
  # Make sure date range is valid
  first = min(concdat$Date, na.rm = T) # earliest data point
  last = max(concdat$Date, na.rm = T) # latest data point
  if(t1 < first){
    stop(paste("No data before t1 to interpolate. Choose t1 >=", first))
  } else if(t2 > last){
    stop(paste("No data after t2 to interpolate. Choose t2 <=", last))
  }
  
  # Make daily resolution
  dates = data.frame(Date = as.Date(as.numeric(first) : as.numeric(last),
                                    origin = "1970-1-1"))
  cd1 = merge(concdat, dates, all = T)
  fd1 = merge(flowdat, dates, all = T)
  
  cd1 = aggregate(. ~ Date, data = cd1, mean, na.action = na.pass) # make it daily resolution
  fd1 = aggregate(. ~ Date, data = fd1, mean, na.action = na.pass) # make it daily resolution
  
  cd1$value.interp = na.approx(cd1$value, na.rm = F) # now have daily resolution concs
  fd1$CFS.interp = na.approx(fd1$CFS, na.rm = F) # and daily resolution flows
  
  cd2.0 = subset(cd1, Date >= t1 & Date <= t2) # date range of interest
  cd2 = ddply(cd2.0, "Date", function(x) colMeans(x[, -1])) # Average multiple concentrations from single day
  fd2.0 = subset(fd1, Date >= t1 & Date <= t2)
  fd2 = ddply(fd2.0, "Date", function(x) colMeans(x[, -1])) # Average multiple flows from single day
  
  naflow = sum(is.na(fd2$CFS))  # number of missing flow data
  if(naflow > 0){
    warning(paste(ifelse(flow.interp, "Interpolated", "Missing"), 
                  naflow, "out of", nrow(fd2), "flow values"))
  }
  
  if(flow.interp) fd2$CFS = fd2$CFS.interp
  
  dat3 = cbind(cd2, fd2["CFS"])
  dat3 = dat3[!is.na(dat3$CFS),]
  
  dat4 = within(dat3, {
    load.meas = CFS * value * 3600 * 24 * 28.317 / 1000000 # convert to kg/day
    load.pred = CFS * value.interp * 3600 * 24 * 28.317 / 1000000 # convert to kg/day
    cumuload = cumsum(load.pred)
  })
  
  dat.out = dat4[c("Date", "value", "CFS", "load.meas", "value.interp", "load.pred")]
  names(dat.out) = mapvalues(names(dat.out), from = "value.interp", to = "pred")
  
  dat.out
}


stepLoad = function(data, t1, t2, flowcol = "CFS", sampleTypes = "all", 
                      flow.interp = T, valcol = "value"){
  ## data$valcol MUST contain only desired concentrations 
  ## (i.e. subset data BEFORE calling this function)
  ## data$flowcol should have near-daily resolution
  ## sampleTypes can be either "all" or a character vector of values from data$sampleType
  require(plyr)
  
  names(data) = mapvalues(names(data), from = c(flowcol, valcol),
                          to = c("CFS", "value"))
  
  concdat = data[!is.na(data$value), c("Date", "value", "sampleType")] # only concentration rows
  
  if(sampleTypes != "all"){
    concdat = subset(concdat, sampleType %in% sampleTypes)
  }
  flowdat = data[!is.na(data$CFS), c("Date", "CFS")] # only flow data
  
  # Make sure date range is valid
  first = min(concdat$Date, na.rm = T) # earliest data point
  last = max(concdat$Date, na.rm = T) # latest data point
  if(t1 < first){
    stop(paste("No data before t1 to interpolate. Choose t1 >=", first))
  } else if(t2 > last){
    stop(paste("No data after t2 to interpolate. Choose t2 <=", last))
  }
  
  # Make daily resolution
  dates = data.frame(Date = as.Date(as.numeric(first) : as.numeric(last),
                                    origin = "1970-1-1"))
  cd1 = merge(concdat, dates, all = T)
  fd1 = merge(flowdat, dates, all = T)
  
  cd1 = aggregate(. ~ Date, data = cd1, mean, na.action = na.pass) # make it daily resolution
  fd1 = aggregate(. ~ Date, data = fd1, mean, na.action = na.pass) # make it daily resolution
  
#   browser()
  cd1$value.interp = midpointStep(x = cd1$Date, y = cd1$value, xout = dates$Date)$y


  fd1$CFS.interp = na.approx(fd1$CFS, na.rm = F) # and daily resolution flows
  
  cd2.0 = subset(cd1, Date >= t1 & Date <= t2) # date range of interest
  cd2 = ddply(cd2.0, "Date", function(x) colMeans(x[, -1])) # Average multiple concentrations from single day
  fd2.0 = subset(fd1, Date >= t1 & Date <= t2)
  fd2 = ddply(fd2.0, "Date", function(x) colMeans(x[, -1])) # Average multiple flows from single day
  
  naflow = sum(is.na(fd2$CFS))  # number of missing flow data
  if(naflow > 0){
    warning(paste(ifelse(flow.interp, "Interpolated", "Missing"), 
                  naflow, "out of", nrow(fd2), "flow values"))
  }
  
  if(flow.interp) fd2$CFS = fd2$CFS.interp
  
  dat3 = cbind(cd2, fd2["CFS"])
  dat3 = dat3[!is.na(dat3$CFS),]
  
  dat4 = within(dat3, {
    load.meas = CFS * value * 3600 * 24 * 28.317 / 1000000 # convert to kg/day
    load.pred = CFS * value.interp * 3600 * 24 * 28.317 / 1000000 # convert to kg/day
    cumuload = cumsum(load.pred)
  })
  
  dat.out = dat4[c("Date", "value", "CFS", "load.meas", "value.interp", "load.pred")]
  names(dat.out) = mapvalues(names(dat.out), from = "value.interp", to = "pred")
  
  dat.out
}

midpointStep = function(x, y, xout, ...){
  # calls approx function from stats
  # automatically removes NAs
  
  df = data.frame(x = x, y = y)
  df = na.omit(df)
  
  x = df$x; y = df$y
  
  xmids = x[1:(length(x)-1)] + diff(x) / 2
  xmids = c(xmids, xout[length(xout)])
  
  out1 = approx(x = xmids, y = y, xout = xout, method = "constant", rule = 2, f = 1)
  out1
}

# set.seed(73)
# x = runif(15) * 100; x = x[order(x)]
# y = runif(15)
# test = midpointStep(x = x, y = y, xout = 0:100)
# plot(x,y)
# lines(test)

findCV = function(name, indices) {
  # used for getting best CV from a list of CV output
  
  a = get(name)
  indNames = sapply(strsplit(names(indices), "\\."), 
                    function(x) paste0(x[-1], collapse = "."))
  a.out = a[names(a) %in% indNames]
  #       print(names(a.out))
  a.out[[1]]
} 

CrossVal = function(model, valdata = NULL, type = "conc", log = F, plot = F) {
  #
  # Performs a leave-one-out cross-validation for the given model.
  # Valdata, if given, must include columns corresponding to model variables
  #
  # Now option to plot calls a separate function
  
  frm = model$terms
  if("gam" %in% class(model)) fit = gam else 
    if(class(model) == "lm") fit = lm else
      stop("model is of wrong class! must be either 'lm' or 'gam'")
    
  
  if(is.null(valdata)){
    valdata = model$model
    attributes(valdata) = within(attributes(valdata), {
      rm(terms)
    })
  }
  
  nfolds <- nrow(valdata)
  case.folds <- rep(1:nfolds,length.out=nrow(valdata))
  
  resids = rep(0, nfolds)
  
  RSS <- 0
  depvar = names(attr(terms(model), "dataClasses"))[1]
  
  if(log){
    meas = exp(valdata[, depvar])
    fitted = exp(model$fitted.values)
    smear = mean(meas, na.rm = T) / mean(fitted, na.rm = T) # smearing estimator
  } else{
    meas = valdata[, depvar]
  }
  
  if(type == "load"){
    meas = meas * valdata$CFS * 3600 * 24 * 28.317 / 1000000
  }
  
  TSS = sum((meas - mean(meas, na.rm = T))^2)
  
  # iterate over all data points in valdata
  if(log) {
    for (fold in 1:nfolds) {
      # What are the training cases and what are the test cases?
      train <- valdata[case.folds!=fold,]
      test <- valdata[case.folds==fold,]
      train.value = exp(train[, depvar])
      test.value = exp(test[, depvar])
      
#       browser()
      # Fit to the training set
      current.fit <- fit(frm, data=train)
      # Predict on the test set
      train.fitted = exp(current.fit$fitted.values)
      smear = mean(train.value, na.rm = T) / mean(train.fitted, na.rm = T) # smearing estimator
      prediction <- smear * exp(predict(current.fit, newdata=test))
      
      
      if(type == "load"){
        prediction = prediction * test$CFS * 3600 * 24 * 28.317 / 1000000
        test.value = test.value * test$CFS * 3600 * 24 * 28.317 / 1000000
      }
      resids[fold] = test.value - prediction
    }
  } else {
    for (fold in 1:nfolds) {
      # What are the training cases and what are the test cases?
      train <- valdata[case.folds!=fold,]
      test <- valdata[case.folds==fold,]
      # Fit to the training set
      current.fit <- fit(frm, data=train)
      # Predict on the test set
      prediction <- predict(current.fit, newdata=test)
      meas = test[, depvar]
      
      if(type == "load"){
        prediction = prediction * test$CFS * 3600 * 24 * 28.317 / 1000000
        meas = meas * test$CFS * 3600 * 24 * 28.317 / 1000000
      }
      resids[fold] =  meas - prediction
    }  
  }
  
  RSS <-sum((resids)^2)
  R2 = 1 - RSS / TSS
  RMS = sqrt(RSS / nfolds)
  valdata$residuals = resids
  
  cv = list(R2 = R2, RMS = RMS, residuals = resids, val.dat = valdata)
  
  if(plot){
    plotCV(cv)
    invisible(cv)
  }
  else cv
}

plotCV = function(CV, ...){
  # 
  # Makes CV plot including normal QQ, histogram, scatter (vs. index), ACF
  # CV should be output from CrossVal; ... are passed to par()
  require(car)
  cvr = CV$resid
  old.par = par(no.readonly = T)
  on.exit(par(old.par))
  par(bg = "white", mfrow = c(2, 2), mar = c(3, 5, 2, 2), ...)
  qqPlot(cvr)
  hist(cvr, main = "")
  plot(cvr)
  acf(cvr, main = "")
}

kfoldCV = function(model, valdata = NULL, type = "conc", log = F, nfolds = 10, 
                   seed = NULL) {
  #
  # Performs a k-fold cross-validation using model and valdata
  # 
  
  frm = model$terms
  if("gam" %in% class(model)) fit = gam else 
    if(class(model) == "lm") fit = lm else
      stop("model is of wrong class! must be either 'lm' or 'gam'")
  
  if(is.null(valdata)){
    valdata = model$model
  }
  
  if(is.numeric(seed)) set.seed(seed)
  
  case.folds <- rep(1:nfolds, length.out=nrow(valdata))
  #   # divide the cases as evenly as possible
  case.folds <- sample(case.folds) # randomly permute the order
  
  depvar = names(attr(terms(model), "dataClasses"))[1]
  
  out.dat = data.frame(fold = numeric(0), TSS = numeric(0), RSS = numeric(0),
                       RMSE = numeric(0), R2 = numeric(0), n = numeric(0))
  
  if(log){
    # iterate over all data points in valdata
    for (fold in 1:nfolds) {
      # What are the training cases and what are the test cases?
      train <- valdata[case.folds!=fold,]
      test <- valdata[case.folds==fold,]
      train.value = exp(train[, depvar])
      test.value = exp(test[, depvar])
      
      # Fit to the training set
      current.fit <- fit(frm, data=train)
      # Predict on the test set
      train.fitted = exp(current.fit$fitted.values)
      smear = mean(train.value, na.rm = T) / mean(train.fitted, na.rm = T) # smearing estimator
      prediction <- smear * exp(predict(current.fit, newdata=test))
      
      if(type == "load"){
        prediction = prediction * test$CFS * 3600 * 24 * 28.317 / 1000000
        test.value = test.value * test$CFS * 3600 * 24 * 28.317 / 1000000
      }
      
      n = nrow(test)
      TSS = sum((test.value - mean(test.value, na.rm = T))^2)
      RSS = sum((test.value - prediction)^2)
      RMSE = sqrt(RSS / n)
      R2 = 1 - RSS / TSS
      
      outrow = data.frame(fold, TSS, RSS, RMSE, R2, n)
      out.dat = rbind(out.dat, outrow)
    }
    
    RSS <-sum(out.dat$RSS)
    TSS = sum(out.dat$TSS)
    R2 = 1 - RSS / TSS
    RMSE = sqrt(RSS / nrow(valdata))
  } 
  else { 
    # iterate over all data points in valdata
    for (fold in 1:nfolds) {
      # What are the training cases and what are the test cases?
      train <- valdata[case.folds!=fold,]
      test <- valdata[case.folds==fold,]
      # Fit to the training set
      current.fit <- fit(frm, data=train)
      # Predict on the test set
      prediction <- predict(current.fit, newdata=test)
      meas = test[, depvar]
      
      if(type == "load"){
        prediction = prediction * test$CFS * 3600 * 24 * 28.317 / 1000000
        meas = meas * test$CFS * 3600 * 24 * 28.317 / 1000000
      }
#       browser()
      resids = meas - prediction
      
      n = nrow(test)
      TSS = sum((meas - mean(meas, na.rm = T))^2)
      RSS = sum(resids^2)
      RMSE = sqrt(RSS / n)
      R2 = 1 - RSS / TSS
      
      outrow = data.frame(fold, TSS, RSS, RMSE, R2, n)
      out.dat = rbind(out.dat, outrow)
    }
    
    RSS <-sum(out.dat$RSS)
    TSS = sum(out.dat$TSS)
    R2 = 1 - RSS / TSS
    RMSE = sqrt(RSS / nrow(valdata))
  }
  return(list(R2 = R2, RMS = RMSE, cv.dat = out.dat, mean.R2 = mean(out.dat$R2),
              val.dat = valdata))
}


getBestModIndex = function(modList, cvList, metric = "R2", obj = c("max", "min")){
  obj = match.arg(obj)
#   browser()
  metrics = sapply(cvList, "[[", metric)
  best = which(metrics == max(metrics))

  return(best)
}

depvar = function(model){
  return(names(attr(terms(model), "dataClasses"))[1])
}

islog = function(model, logname = "lc", linname = "value"){
  if(depvar(model) == logname) return(TRUE) else
    if(depvar(model) == linname) return(FALSE) else
      stop("Unrecognized response variable")
}


loadMSE = function(conc.pred, conc.se, conc.residvar, flow, method = c("empirical", "se.fit"), 
                   use = c("chat", "lhat"), totvar = T){
  # 
  # Estimates standard error for loads estimated with models.
  # 
  # An earlier method took a data frame with certain named columns as input;
  # now two vectors, load.pred and load.se, are given as input
  #
  # if totvar == T, uses the law of total variance to estimate the variance. 
  # Otherwise variance misses residual variance component. I'm keeping this apparent bug as an option
  # just so that I can compare the difference in standard error estimates. 
  #
  method = match.arg(method)
  use = match.arg(use)
  
  adjflow = flow * 3600 * 24 * 28.317 / 1000000 # put into million liters per day for load conversion
  n = nrow(conc.pred)
  if(use == "chat"){
    if(method == "se.fit"){
      sigmat = outer(conc.se, conc.se) # sigmas
      corvec = acf(na.omit(conc.pred), lag.max = n - 1, plot = F)$acf # estimated autocorrelation
      cormat = toeplitz(as.numeric(corvec))
      covmat1 = cormat * sigmat
      
    } else if(method == "empirical"){
      covvec = acf(na.omit(conc.pred), lag.max = n - 1, plot = F,
                   type = "covariance")$acf 
      covmat1 = toeplitz(as.numeric(covvec))
      if(totvar) {diag(covmat1) = diag(covmat1) + conc.residvar} # add variance component from residuals
    } 
    coefmat = outer(adjflow, adjflow)
    covmat = covmat1 * coefmat
  }
  
  
  else if(use == "lhat") {
    load.pred = conc.pred * adjflow
    load.se = conc.se *adjflow
    
    if(method == "se.fit"){
      corvec = acf(na.omit(load.pred), lag.max = n - 1, plot = F)$acf 
      cormat = toeplitz(as.numeric(corvec)) 
      varmat = outer(load.se, load.se)
      covmat = cormat * varmat
    } else if(method == "empirical"){
      covvec = acf(na.omit(load.pred), lag.max = n - 1, plot = F,
                   type = "covariance")$acf 
      covmat = toeplitz(as.numeric(covvec))
      if(totvar){diag(covmat) = diag(covmat) + conc.residvar * adjflow^2} # add residual component of variance
    }
  }
  #   browser()
  totvar = sum(covmat)
  totsd = sqrt(totvar)
  totsd
}
