# functions_comparing.r
# Functions for comparing load estimates, plotting comparisions
# Originally made for loadDoc.Rmd, originally located in functions.r.

#### Load estimation and comparison ####

modelFormsTable = function(lm.model, gam.model){
  require(odfWeave)
  
  resp.lm = all.vars(terms(lm.model))[1]
  resp.gam = all.vars(terms(gam.model))[1]
  respcol = matrix(c(resp.lm, resp.gam), nr = 2, dimnames = list(c("lm", "gam"), "response"))
  
  lm.p = getPvals(lm.model)
  gam.p = getPvals(gam.model)
  
  termsTable = matrix(rep("", max(length(lm.p), length(gam.p)) * 2), nr = 2)
  termsTable[1, 1:length(lm.p)] = names(lm.p)
  termsTable[2, 1:length(gam.p)] = names(gam.p)
  
  pTable = matrix(rep("", max(length(lm.p), length(gam.p)) * 2), nr = 2)
  pTable[1, 1:length(lm.p)] = formatP(lm.p)
  pTable[2, 1:length(gam.p)] = formatP(gam.p)
  pTable = apply(pTable, 2, function(x) paste0("(", x, ")"))
  
  dof1 = 1 + length(lm.p)
  dof2 = sum(gam.model$edf)
  dofcol = matrix(sprintf("%.4f", c(dof1, dof2)), nr = 2, 
                  dimnames = list(c("lm", "gam"), "(E)DF"))
  
  formtable = matrixPaste(termsTable, pTable)
  colnames(formtable) = paste("term", 1:ncol(formtable))
  
  formtable = cbind(respcol, formtable, dofcol)
}

plotFitTS = function(data, predCol, measCol, is.load = F, ...){
  data$pc = data[predCol]
  data$mc = data[measCol]
  ymax = max(data[c(predCol, measCol)], na.rm = T)
  ymin = min(data[c(predCol, measCol)], na.rm = T)
  ymax = ymax + (ymax - ymin) * 0.1
  ymin = ymin - (ymax - ymin) * 0.1
  ylab = ifelse(is.load, "load (kg/day)", "concentration (mg/L)")
  
  old.par = par(no.readonly = T)
  par(bg = "white", mfrow = c(2, 1), mar = c(2, 4, 1, 2), oma = c(2, 1, 2, 1))
  plot(data[, predCol] ~ Date, data, type = "l", ylab = ylab, 
       ylim = c(ymin, ymax))#, ...)
  points(data[, measCol] ~ Date, data)
  plot(CFS ~ Date, data, type = "l")
  par(old.par)
}

plotFitScatter = function(data, type = c("lm", "gam"), ...){
  # 
  # modified from version in functions.r
  # type should be "lm" or "gam". 
  #
  type = match.arg(type)
  data$conc = data[, paste0("conc.", type)]
  data$load = data[, paste0("load.", type)]
  
  old.par = par(no.readonly = T)
  par(bg = "white", mfrow = c(1, 2), mar = c(4, 4, 1, 2))
  plot(conc.meas ~ conc, data, xlab = "modeled concentration (mg/L)", 
       ylab = "measured concentration (mg/L)", ...)
  abline(0, 1)
  plot(load.meas ~ load, data, xlab = "modeled load (kg/day)", 
       ylab = "measured load (kg/day)", ...)
  abline(0, 1)
  par(old.par)
}

residPlotsConc = function(data, type = c("lm", "gam")){
  
  type = match.arg(type)
  
  data$conc = data[, paste0("conc.", type)]
  data$resid = with(data, conc - conc.meas)
  
  old.par = par(no.readonly = T)
  par(bg = "white", mfrow = c(2, 2), mar = c(3, 3, 2, 2), oma = c(2, 1, 2, 1))
  
  plot(resid ~ Date, data)
  hist(data$resid, main = "")
  qqPlot(data$resid)
  acf(na.omit(data$resid), main = "")
  
  par(old.par)
  
  data$resid
}

residPlotsLoad = function(data, type = c("lm", "gam")){
  type = match.arg(type)
  data$load = data[, paste0("load.", type)]
  data$resid = with(data, load - load.meas)
  
  old.par = par(no.readonly = T)
  par(bg = "white", mfrow = c(2, 2), mar = c(3, 3, 2, 2), oma = c(2, 1, 2, 1))
  
  plot(resid ~ Date, data)
  hist(data$resid, main = "")
  qqPlot(data$resid)
  acf(na.omit(data$resid), main = "")
  par(old.par)
  
  data$resid
}

residScatterConc = function(data, type = c("lm", "gam")){
  type = match.arg(type)
  data$conc = data[, paste0("conc.", type)]
  data$load = data[, paste0("load.", type)]
  
  old.par = par(no.readonly = T)
  par(bg = "white", mfrow = c(1, 2))
  plot(conc - conc.meas ~ conc.meas, data, ylab = "concentration residuals", 
       xlab = "measured concentration (mg/L)")
  abline(0, 0)
  plot(conc - conc.meas ~ CFS, data, ylab = "concentration residuals", 
       xlab = "measured flow (CFS)")
  abline(0, 0)
  par(old.par)
}

residScatterLoad = function(data, type = c("lm", "gam")){
  type = match.arg(type)
  data$load = data[, paste0("load.", type)]
  
  old.par = par(no.readonly = T)
  par(bg = "white", mfrow = c(1, 2))
  plot(load - load.meas ~ load.meas, data, ylab = "load residuals", 
       xlab = "measured load (kg/day)")
  abline(0, 0)
  plot(load - load.meas ~ CFS, data, ylab = "load residuals", 
       xlab = "measured flow (CFS)")
  abline(0, 0)
  par(old.par)
}

relDif = function(series1, series2, log = T){
  # 
  # Creates 2-d filled contour plot for relative differences between two 
  #  data frames, each with 2 columns, "Date" and a numeric value.
  #  Series 1 is reference, to which series2 is compared
  #  Value plotted is percent difference
  #
  require(caTools)
  
  daysOfInterest = c(1, 7, 14, 30, 60, 90, 180, 365)
  base = 1.05
  daysForLoad = base ^ (0:ceiling(log(365, base = base)))
  
  data = merge(series2, series1, by = "Date", all = F)
  data$dif = data[,2] - data[,3]
  outmat = matrix(nrow = 0, ncol = nrow(data))
  for(i in daysForLoad){
    avdif = runmean(data$dif, k = round(i), align = "center", endrule = "NA")
    avload = runmean(data[,3], k = round(i), align = "center", endrule = "NA")
    reldif = avdif / avload * 100
    outmat = rbind(outmat, reldif)
  } # get running mean for different window sizes
  rownames(outmat) = daysForLoad
  colnames(outmat) = data$Date
  
  date1 = data$Date[1]
  date2 = data$Date[nrow(data)]
  dates = seq.Date(date1, date2, by = "month")
  
  #   outmat[is.na(outmat)] = -110
  #   cols = c("black", cm.colors(19), "red")
  cols = c("black", rainbow(19), "black")
  cols[11] = "lightgray"
  allbreaks = c(-11:-1 * 10, 1:10 * 10, max(110, ceiling(max(outmat, na.rm = T))))
  
  par.old = par(no.readonly = T)
  #   par(bg = "gray")
  if(log){
    # log plot
    filled.contour(as.numeric(data$Date), log(daysForLoad), t(outmat), 
                   levels = allbreaks, 
                   xlab = "Date", ylab = "length of load time period (days)",
                   col = cols, key.title = title(main = "percent\ndif."),
                   plot.axes = { axis(1, at = as.numeric(dates), 
                                      labels = format(dates, "%m/%Y"))
                                 axis(2, at = log(daysOfInterest), 
                                      labels = daysOfInterest) })
#                    plot.title = title(main = paste(series1, series2)))  
  } else {
    #linear plot
    filled.contour(as.numeric(data$Date), daysForLoad, t(outmat), 
                   levels = allbreaks, 
                   xlab = "Date", ylab = "length of load time period (days)",
                   col = cols, key.title = title(main = "percent\ndif."),
                   plot.axes = { axis(1, at = as.numeric(dates),
                                      labels = format(dates, "%m/%Y"))
                                 axis(2, at = daysOfInterest, 
                                      labels = daysOfInterest) })
#                    plot.title = title(main = paste(series1, series2)))  
  }
  par(par.old)
  invisible(outmat)
}

relDifHists = function(reldif){
  old.par = par(no.readonly = T)
  par(bg = "white", mfrow = c(2, 3), mar = c(3, 3, 2, 0))
  hist(reldif[1,], main = "1-day loads")
  hist(reldif[41,], main = "7-day loads")
  hist(reldif[71,], main = "30-day loads")
  hist(reldif[94,], main = "90-day loads")
  hist(reldif[108,], main = "180-day loads")
  hist(reldif[122,], main = "1-year loads")
  par(old.par)
}

rdMeanSD = function(reldif){
  old.par = par(no.readonly = T)
  par(bg = "white", mfrow = c(1, 2))
  meandif = apply(reldif, 1, mean, na.rm = T)
  sddif = apply(reldif, 1,sd, na.rm = T)
  plot(as.numeric(names(meandif)), meandif, type = "l", 
       xlab = "length of load interval (days)", 
       ylab = "mean difference (percent)")
  plot(as.numeric(names(sddif)), sddif, type = "l", 
       xlab = "length of load interval (days)", 
       ylab = "st.dev of difference (percent)")
  par(old.par)
}

plotEstDif = function(data, type1, type2, is.load, ...){
  # type1, type2 from c("interp", "model", "usgs")
  # difference take is type1 - type2
  
  pref = ifelse(is.load, "load", "conc")
  
  d1 = data[, paste(pref, type1, sep = ".")] - 
    data[, paste(pref, type2, sep = ".")]
  d2 = data[, paste(pref, type1, sep = ".")] - 
    data[, paste(pref, "meas", sep = ".")]
  ylab = ifelse(is.load, "load difference (kg/d)", "conc difference (mg/L)")
  
  old.par = par(no.readonly = T)
  par(bg = "white", mfrow = c(2, 1), mar = c(0, 4, 0, 2), oma = c(4, 1, 3, 1))
  plot(d1 ~ data$Date, type = "l", ylab = ylab, ...)
  points(d2 ~ data$Date)
  abline(0, 0, lty = 3)
  
  plot(CFS ~ Date, data, type = "l", col = "blue", ylab = "flow (CFS)")
  par(old.par)
}

load_table = function(tsResults, seMethod = c("empirical", "se.fit"), use = c("chat", "lhat"),
                      totvar = TRUE){
  #
  # tsResults is of the form output by getPredictions()
  #
  
  require(dplyr)
  require(reshape2)
  seMethod = match.arg(seMethod)
  use = match.arg(use)
  
  loads = tsResults %>%
    select(contains("load"), -contains("se.load"), -contains("meas")) %>%
    apply(2, sum)

  
  loads.m = melt(loads, value.name = "load_kg")
  loads.m$method = rownames(loads.m)
  loads.m$SE_kg = NA
  
  gam.residvar = with(tsResults, var(conc.meas - conc.gam, na.rm = T))
  lm.residvar = with(tsResults, var(conc.meas - conc.lm, na.rm = T))
  
  gam.se = loadMSE(conc.pred = tsResults$conc.gam, conc.se = tsResults$se.conc.gam, 
                   conc.residvar = gam.residvar,
                   flow = tsResults$CFS, method = seMethod, use = use, totvar = totvar)
  lm.se = loadMSE(conc.pred = tsResults$conc.lm, conc.se = tsResults$se.conc.lm, 
                  conc.residvar = lm.residvar,
                  flow = tsResults$CFS, method = seMethod, use = use, totvar = totvar)
  
  loads.m["load.gam", "SE_kg"] = gam.se
  loads.m["load.lm", "SE_kg"] = lm.se
  
  loads.m
}

load_barplot = function(loadTable){
  
  require(ggplot2)
  
  loads.m = loadTable
  
  # Error bars represent standard error of the mean
  bp = ggplot(loads.m, aes(x=method, y=load_kg)) + 
        geom_bar(position=position_dodge(), stat="identity") +
        geom_errorbar(aes(ymin=load_kg-2*SE_kg, ymax=load_kg+2*SE_kg),
                      width=.2,                    # Width of the error bars
                      position=position_dodge(.9)) # +
    #     theme(axis.text.x = element_text(angle = 90))
  bp
}