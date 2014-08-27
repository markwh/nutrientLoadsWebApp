# metaFunctions.r
# high-level functions used in main program

source("functions/functions2.r")

makeData = function(caldata, preddata, startDate, endDate, station, variable,
                    sset = c("allSamples", "noWetRoutine", "noWetGrab")){
  #
  # Mostly just calls makeModelData, but first gets saved number of outliers
  #
  if(outliers == "remove"){
    outs_hi = read.csv("outliers_high.csv", row.names = 1)[station, variable]
    outs_lo = read.csv("outliers_low.csv", row.names = 1)[station, variable]
  } else outs_hi = outs_lo = 0

  td = makeModelData(caldata, preddata, location = station, var = variable, 
                     startDate = startDate, endDate = endDate, outs_hi = outs_hi, 
                     outs_lo = outs_lo, sset = sset)
  td
}


fitLMs = function(cal.dat, modelForms, ...){
  # fits linear models, 
  # returns a list of linear models INCLUDING response variable $y in output
  slimFit = function(form) lm(formula(form), data = cal.dat, y = T, ...)
  out.lm = lapply(modelForms, slimFit)
  out.lm
}


fitGAMs = function(cal.dat, modelForms, ...){
  # fits GAMs, 
  # returns a list of GAM models
  require(mgcv)
  slimFit = function(form) gam(formula(form), data = cal.dat, ...) 
  out.gam = lapply(modelForms, slimFit)
  out.gam
}

predictFromModels = function(modelList, predData, bound.conc, ...){
  # 
  # ... gets passed to predictLoad (functions2.r), e.g. 
  #
  
  loglist = lapply(modelList, islog)
  
  pl = function(mod, islog) { 
    predictLoad(predData, model = mod, log = islog, 
                bound.conc = bound.conc, ...)
  } 
  
  out = mapply(pl, mod = modelList, islog = loglist, SIMPLIFY = F, ...)
  
  out
}

predInterps = function(data, startDate, endDate, ...){
  # uses interpLoad in functions2.r
  
  int.mon = interpLoad(data, t1 = startDate, t2 = endDate, 
                       sampleTypes = "monthly grab", ...)
  
  int.all = interpLoad(data, t1 = startDate, t2 = endDate, 
                       sampleTypes = "all", ...)
  step.mon = stepLoad(data, t1 = startDate, t2 = endDate, 
                      sampleTypes = "monthly grab", ...)
  step.all = stepLoad(data, t1 = startDate, t2 = endDate, 
              sampleTypes = "all")
  
  return(list(int.mon = int.mon, int.all = int.all, step.mon = step.mon, 
              step.all = step.all))
}

# runLoadests

cv_allkinds = function(modlist, valdat, seed = NULL){
  # 
  # Gets all kinds of cross-validation from a list containing a single Gam and LM object (one of each)
  #
  
  Gam = modlist[sapply(modlist, function(x) is(x, "gam"))][[1]]
  Lm =  modlist[sapply(modlist, function(x) {is(x, "lm") & !is(x, "gam")})][[1]]
  out = list()
  outnames = character(0)
  methods = c("loo", "kf"); modeltypes = c("Lm", "Gam"); types = c("conc", "load")
  
  i = 1
  for(method in c("loo", "kf")){
    k = ifelse(method == "loo", "all", 10)
    for(modeltype in c("Lm", "Gam")){
      for(type in c("conc", "load")){
        thisname = paste0(c(method, modeltype, type), collapse = "")
        thisone = cv(get(modeltype), k = k, valdata = valdat, type = type, 
                  log = NULL, plot = F, seed = seed)
        outnames = c(outnames, thisname)
        out[[i]] = thisone
        i = i + 1
      }
    }
  }
  names(out) = outnames
  out
}

cv = function(model, k = "all", valdata = NULL, type = c("conc", "load"), 
              log = NULL, plot = F, seed = NULL){
  #
  # Wrapper for CrossVal and kfoldCV functions in functions.r, generalizes
  # so that either k-fold or leave-one-out can be used.
  # folds argument can be either "all" (for leave-one-out), or an integer
  # greater than 1.
  
  type = match.arg(type)
  if(is.null(log)){
    log = names(attr(terms(model), "dataClasses"))[1] == "lc"
  }
  if(k == "all") {
    out = CrossVal(model, valdata, type, log, plot)
  } else {
    out = kfoldCV(model, valdata, type, log, nfolds = k, seed = seed)
  } 
}

findBiases = function(obsdata, lmList, gamList){
  obsLoad = with(obsdata, value  * CFS * 3600 * 24 * 28.317 / 1000000)
  meanObsLoad = mean(obsLoad, na.rm = T)
  sdMean = sd(obsLoad, na.rm = T) / sqrt(length(obsLoad))
  
  getBias = function(modelLoad){
    modMean = mean(modelLoad, na.rm = T)
    biasPct = (modMean / meanObsLoad - 1) * 100
    biasDist = (modMean - meanObsLoad) / sdMean # Bias in stdevs from the sample mean
    return(c(biasPct = biasPct, biasDist = biasDist))
  }
  
  lmBiases = lapply(lmList, function(x) getBias(x$load.pred))
  gamBiases = lapply(gamList, function(x) getBias(x$load.pred))
  
  return(list(lmBiases = lmBiases, gamBiases = gamBiases))
}
# 
getBestModels = function(lmlist, gamlist, lmcv, gamcv, metric = "R2", tidy = T){
  #
  # gets best models from each of the model lists, using the CV lists. 
  # uses getBestMod function in functions2.r
  # if tidy is set to true, then non-significant predictors are removed in a stepwise fashion
  
  lm.ind = getBestModIndex(lmlist, lmcv, metric)
  gam.ind = getBestModIndex(gamlist, gamcv, metric)
  
  if(length(lm.ind) > 1){
    message("Multiple best lm models selected. Selecting best one based on df")
    dfs = sapply(lm.ind, function(x) length(lmlist[[x]]$coefficients))
    lm.ind = lm.ind[which(dfs == min(dfs))] 
    
    if(length(lm.ind) > 1){
      message("Multiple best lm models have same df. Selecting first one arbitrarily.")
      lm.ind = lm.ind[1]
    }
  }
  
  if(length(gam.ind) > 1){
    message("Multiple best gam models selected. Selecting best one based on edf")
    edfs = sapply(gam.ind, function(x) sum(gamlist[[x]]$edf))
    gam.ind = gam.ind[which(edfs == min(edfs))]        
    if(length(gam.ind) > 1){
      message("Multiple best gam models have same df. Selecting first one arbitrarily.")
      gam.ind = gam.ind[1]    }
  }
  
  bestLm = lmlist[[lm.ind]]; bestGam = gamlist[[gam.ind]]
  
  if(tidy){bestLm = tidyModel(bestLm); bestGam = tidyModel(bestGam)}
  
  return(list(bestLm = bestLm, bestGam = bestGam, 
              indices = c(bestLm = lm.ind, bestGam = gam.ind)))
}


getPredictions = function(lm.model, gam.model, obsdata, preddata, interpdata,
                          bound.conc = T, ...){
#   browser()
  lmPred = predictLoad(preddata, model = lm.model, log = islog(lm.model), 
                       bound.conc = bound.conc, ...)
  
  gamPred = predictLoad(preddata, model = gam.model, log = islog(gam.model), 
                        bound.conc = bound.conc, ...)
  
  interpPred = predInterps(interpdata, startDate=min(preddata$Date), 
                           endDate = max(preddata$Date))
  
  interpCols = c("int.mon.pred", "int.mon.load.pred", 
                 "int.all.pred", "int.all.load.pred", "step.mon.pred",
                 "step.mon.load.pred", "step.all.pred", "step.all.load.pred")
  #   browser()
  # make data frames to merge (by Date):
  .cfs = preddata[c("Date", "CFS")]
  .obs = obsdata[c("Date", "value", "sampleType")]
  names(.obs)[2] = "conc.meas"
  .lm = lmPred[c("Date", "pred", "se.pred", "load.pred", "se.load")]
  
  .gam = gamPred[c("Date", "pred", "se.pred", "load.pred", "se.load")]
  .interp = within(as.data.frame(interpPred), {
    Date = int.mon.Date
#     load.meas = int.mon.load.meas
    })[c("Date", interpCols)]
#     browser()
  .models = merge(.lm, .gam, by = "Date", all = T, suffixes = c(".lm", ".gam"))
  
  lall = Reduce(function(...) merge(..., all = T), list(.cfs, .obs, .models, .interp)) 
  lall$load.meas = with(lall, conc.meas * CFS * 3600 * 24 * 28.317 / 1000000)
  lall = lall[c("Date", "CFS", "conc.meas", "load.meas", "sampleType",
                names(.models)[-1], interpCols)]
  names(lall) = gsub("load.pred", "load", names(lall))
  names(lall) = gsub("pred", "conc", names(lall))
  mapfrom = outer(c("int.mon", "int.all", "step.mon", "step.all"), 
                  c("conc", "load"), paste, sep = ".")
  mapto = outer(c("conc", "load"), 
                c("int.mon", "int.all", "step.mon", "step.all"), 
                paste, sep = ".")
  names(lall) = mapvalues(names(lall), mapfrom, t(mapto)) 
  lall

}



# compareBestModels
# saveComparisons
# knitReport