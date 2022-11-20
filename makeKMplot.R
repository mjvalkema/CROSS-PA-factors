#Loading packages
library(survival) # function coxph
library(haven)
library(dplyr)
library(tableone)
library(survival)
library(survminer)
library(cowplot)
library(ggplot2)
library(prodlim)
library(extrafont) # useful??
library(ggpubr)
library(lubridate)

makeKMplot <- function(data, timevar, eventvar, groupingvar, titl, gvartitle, grouplabels=NA, save=FALSE, filename='temp'){
  
  tvar <- data[[timevar]]
  evar <- data[[eventvar]]
  gvar <- data[[groupingvar]]
  colors <- c("#4E79A7","#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1") # Tableau colors
  
  ## Make plot
  sdat <- tibble(tvar, evar, gvar)
  #sdat <- subset(sdat, sdat$gvar != '<NA>')
  sdat <- subset(sdat, !is.na(sdat$gvar))
  
  # Define graph text
  subt <- ""
  xlab <- "Months since completion of nCRT"
  ylab <- "Percentage of patients"
  gvar <- droplevels(gvar, exclude = c("other", "88", "<NA>", "NA", "na", NA))
  legend.title <- gvartitle
  if (is.na(grouplabels[1])) {
    kmgroups <- levels(gvar)
  } else {
    kmgroups <- grouplabels
  }
  
  #survival fit function
  sfit <- survfit(Surv(tvar, evar) ~ gvar, data=sdat)
  #print(sfit)
  
  #Median follow-up for survivors
  #print("median follow-up for survivors:")
  #print(quantile(survfit(Surv(tvar, evar==0) ~ gvar, data=sdat))[1])
  #print("median follow-up for survivors based on all data:")
  #print(quantile(prodlim(Hist(tvar,evar==0)~1,data=sdat,reverse=FALSE))) # median follow-up survivors all categories

  km <- ggsurvplot(sfit, data=sdat,
                   fun = "pct", # showing percentage instead of survival probability
                   palette= colors, # colours of lines/palettes
                   legend.title = legend.title, # Tile of  comparison groups
                   legend.labs = kmgroups, #names of comparison groups                 
                   title= titl, #title 
                   font.title = c(15, "plain", "Black"),
                   #surv.median.line = "hv",
                   xlab=xlab,
                   pval = TRUE, #"display p value, state TRUE, FALSE or "chosen text"
                   ylab=ylab, #name y axis lable
                   size=1,
                   font.tickslab = c(15, "plain", "Black"), #font x and y label numbers
                   risk.table=TRUE, fontsize = 5, #Display risk table
                   conf.int = FALSE, #show confidence interval
                   surv.median.line = "hv",# add median line
                   censor.shape=73, #shape censor lines
                   censor.size=2, #size censor lines
                   xlim=c(0,60), #min and max displayed on X axis
                   break.x.by=6, #intervals displayed on x interval
                   axes.offset=TRUE, ##small indent x and y axis
                   risk.table.title= "Numbers at risk" ,#title risk table
                   risk.table.y.text= FALSE, #show colour with/without text in group labels 
                   tables.height=0.2, #height of table (handy when many groups)
                   tables.theme=theme_cleantable()) #overall theme of ggplot risk table
  
  km$table <- km$table + theme(plot.title=element_blank()) # remove white row above numbers at risk table
  print(km)

  if (save) {
    png(filename, units = "cm", width=20, height=15, res=1200)
    print(km, newpage = FALSE)
    dev.off()
  }

  return(list("sfit"=sfit, "pvalue"=surv_pvalue(sfit, sdat)$pval))
}


