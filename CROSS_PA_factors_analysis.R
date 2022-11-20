## CROSS PA factors analysis
# Written by M.J. Valkema, November 2022

setwd("/Users/Maartje/repos/CROSS-PA-factors")
source('CROSS_PA_factors_dataprep.R')
colors <- c("#4E79A7","#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1") # Tableau colors

# Load packages
library(tableone) # baseline table
library(RColorBrewer) # to be used with river plots
library(riverplot) # river plots

# Load packages for survival
library(survival) # function coxph
library(haven)
library(dplyr)
library(tableone)
library(survival)
library(survminer)
library(cowplot)
library(ggplot2)
library(prodlim)
library(extrafont)
library(ggpubr)
library(lubridate)
library(splines) # time series

# Script for exporting Cox summaries
source('makeCoxSum.R') # returns fit_alt

# KM plots
source('makeKMplot.R')
gvars <- read.csv(file="data/gvars_UK.csv", header=TRUE, sep = ",", na.strings = "")

######### Datasets #########
# Data is the entire dataset
modDataRes <- subset(data, data$resection == "yes") # only patients who underwent resection
modDataResResidual <- subset(modDataRes, modDataRes$TRG != "TRG1") # for assessment of post-treatment PA factors, which are assessed relative to the residual tumor area
modDataResResidual$res_diffgrade <- droplevels(modDataResResidual$res_diffgrade, exclude = c("NA, TRG1"))
modDataResResidual$res_src_atypical_010 <- droplevels(modDataResResidual$res_src_atypical_010, exclude = c("88"))
modDataResResidual$changeN <- droplevels(modDataResResidual$changeN, exclude = c("NA NA"))

######### Baseline characteristics #########
# Time period of treatment
range(data[data$Institute.Abbreviation != "RAD",]$end_ncrt, na.rm = TRUE) # EMC + PATH patients
range(data[data$Institute.Abbreviation == "RAD",]$end_ncrt, na.rm = TRUE) # RAD patients

# Time period of diagnosis
range(data[data$Institute.Abbreviation != "RAD",]$date_diagnosis) # EMC + PATH patients
range(data[data$Institute.Abbreviation == "RAD",]$date_diagnosis) # RAD patients
range(data$date_last_fu, na.rm=TRUE) # range of last follow-up
summary(data$potentialFU) # potential follow-up time

# Baseline table
baselineVars <- c("sex", "age", "location", "cT", "cN", "biopsies_diffgrade",
                  "cross_all", "cross_less_ct", "cross_less_rt")
catVars <- c("sex", "location", "cT", "cN", "biopsies_diffgrade",
              "cross_all", "cross_less_ct")
tableBl <- CreateTableOne(vars = baselineVars, factorVars = catVars, data = data)
tableBl <- print(tableBl, nonnormal = c("age", "cross_less_rt"), quote = FALSE, noSpaces = TRUE, digits=NULL)
write.csv(tableBl, "output/Table1.csv", row.names = TRUE, na = "")

# Baseline table split by hospital (Erasmus MC vs Radboud UMC)
baselineVars <- c("sex", "age", "location", "cT_cat", "cN_cat", "cross_all", "cross_less_ct", "cross_less_rt",
                   "biopsies_diffgrade", "biopsies_mucin_any", "biopsies_src_any", "biopsies_src_atypical_any", 
                  "biopsies_srcpcc_any", "biopsies_rfall_any")
catVars <- c("sex", "location", "cT_cat", "cN_cat", "cross_all", "cross_less_ct",
             "biopsies_diffgrade", "biopsies_mucin_any", "biopsies_src_any", "biopsies_src_atypical_any", 
             "biopsies_srcpcc_any", "biopsies_rfall_any")
tableBl <- CreateTableOne(vars = baselineVars, factorVars = catVars, strata = "hospital", data = data)
tableBl <- print(tableBl, nonnormal = c("age", "cross_less_rt"), quote = FALSE, noSpaces = TRUE, digits=NULL)
write.csv(tableBl, "output/Table1Hospital.csv", row.names = TRUE, na = "")

# Surgery characteristics
surgVars <- c("surgery", "resection", "specify_no_surgery", "spec_no_surgery_other", "spec_no_res","interval_ncrt_surgery")
catsurgVars <- c("surgery", "resection", "specify_no_surgery", "spec_no_surgery_other", "specify_no_resection")
tableSurg <- CreateTableOne(vars = surgVars, factorVars = catsurgVars, data = data)
tableSurg <- print(tableSurg, nonnormal = c("interval_ncrt_surgery"), quote = FALSE, noSpaces = TRUE) # interval in weeks
write.csv(tableSurg, "output/TableSurgery.csv", row.names = TRUE, na = "")

# General resection pathology characteristics in all patients
pathVars <- c( "TRG", "ypT", "ypn", "pCR", "prepT", "prepn", "radicality", "Rmargin")
tablePA_res <- CreateTableOne(vars = pathVars, factorVars = pathVars, data = modDataRes)
tablePA_res <- print(tablePA_res, quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tablePA_res, "output/TablePAResection_general.csv", row.names = TRUE, na = "")

table(modDataRes$radicality, modDataRes$biopsies_diffgrade) # assess tumor differentiation grade for patients with R1 resections
table(modDataRes$Rmargin, modDataRes$biopsies_diffgrade)

# Baseline pathology characteristics (PA factors)
pathVars <- c("biopsies_mucin", "biopsies_src", "biopsies_src_atypical", "biopsies_srcpcc", "biopsies_rfall", "biopsies_diff")
tablePA_base <- CreateTableOne(vars = pathVars, factorVars = pathVars, data = data)
tablePA_base <- print(tablePA_base, quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tablePA_base, "output/TableBaselinePAfactors.csv", row.names = TRUE, na = "")

# Post-treatment pathology characteristics (PA factors)
modDataRes <- subset(data, data$resection == "yes")
modDataResResidual <- subset(modDataRes, modDataRes$TRG != "TRG1") 
pathVars <- c("res_mucin", "res_src", "res_src_atypical", "res_srcpcc", "res_rfall", "res_diff")
tablePA_post <- CreateTableOne(vars = pathVars, factorVars = pathVars, data = modDataResResidual)
tablePA_post <- print(tablePA_post, quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tablePA_post, "output/TablePostPAfactors.csv", row.names = TRUE, na = "")

############ Baseline characteristics split out by outcome ###############
# All relevant (pathology) characteristics in patients split by PA factor in biopsies (from 1%)
vars <- c("age", "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "biopsies_mucin_any", data = data)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_biopsy_mucin.csv", row.names = TRUE, na = "")

vars <- c("age", "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "biopsies_src_any", data = data)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_biopsy_src.csv", row.names = TRUE, na = "")

vars <- c("age", "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "biopsies_src_atypical_any", data = data)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_biopsy_src_atypical.csv", row.names = TRUE, na = "")

# All relevant (pathology) characteristics in patients split by PA factor in biopsies (from 10%)
vars <- c("age", "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "biopsies_mucin_010", data = data)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_biopsy_mucin_010.csv", row.names = TRUE, na = "")

vars <- c("age", "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "biopsies_src_010", data = data)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_biopsy_src_010.csv", row.names = TRUE, na = "")

vars <- c("age", "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "resection", "pCR", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "biopsies_src_atypical_010", data = data)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_biopsy_src_atypical_010.csv", row.names = TRUE, na = "")

# All relevant (pathology) characteristics in patients who underwent resection split by PA factor (from 1%)
vars <- c("age", "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "res_mucin_any", data = modDataResResidual)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_res_mucin.csv", row.names = TRUE, na = "")

vars <- c("age", "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "res_src_any", data = modDataResResidual)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_res_src.csv", row.names = TRUE, na = "")

vars <- c("age", "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "res_src_atypical_any", data = modDataResResidual)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_res_src_atypical.csv", row.names = TRUE, na = "")

# All relevant (pathology) characteristics in patients who underwent resection split by PA factor (from 10%)
vars <- c("age", "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "res_mucin_010", data = modDataResResidual)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_res_mucin_010.csv", row.names = TRUE, na = "")

vars <- c("age", "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "res_src_010", data = modDataResResidual)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_res_src_010.csv", row.names = TRUE, na = "")

vars <- c("age", "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "res_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "res_src_atypical_010", data = modDataResResidual)
tableSplit <- print(tableSplit, nonnormal = "age", quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tableSplit, "output/TableSplit_res_src_atypical_010.csv", row.names = TRUE, na = "")

# Baseline table split by src_categories (exploratory only)
vars <- c("age", "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
factors <- c( "sex", "biopsies_diffgrade", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "TRG", "radicality")
tableSplit <- CreateTableOne(vars = vars, factorVars = factors, strata = "src_categories",  data = data)
tableSplit <- print(tableSplit, nonnormal = c("age", "cross_less_rt"), quote = FALSE, noSpaces = TRUE, digits=NULL)
write.csv(tableSplit, "output/TableSplit_SRCcategories.csv", row.names = TRUE, na = "")

####### Other tables #########
# Baseline pathology characteristics per TRG category
modDataRes <- subset(data, data$resection == "yes")
pathVars <- c("biopsies_mucin", "biopsies_src", "biopsies_src_atypical", "biopsies_srcpcc", "biopsies_rfall", "biopsies_diffgrade")
tablePA_base <- CreateTableOne(vars = pathVars, factorVars = pathVars, strata = "TRG", data = modDataRes)
tablePA_base <- print(tablePA_base, quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tablePA_base, "output/TableBaselinePA_perTRG.csv", row.names = TRUE, na = "")

# Baseline table for patients with SRC spit by pCR yes or no, to investigate whether some patients with SRC and good response have different baseline characteristics
# This is exploratory only
modDataSRC <- subset(modDataRes, data$biopsies_src_any == "1")
baselineVars <- c("sex", "age", "location", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "biopsies_diffgrade",
                  "hospital", "cross_all", "cross_less_ct", "cross_less_rt")
catVars <- c("sex", "location", "cT_cat", "cN_cat", "prepT_cat", "prepN_cat", "biopsies_diffgrade",
             "hospital", "cross_all", "cross_less_ct")
tableSRC <- CreateTableOne(vars = baselineVars, factorVars = catVars, strata = "pCR",  data = modDataSRC)
tableSRC <- print(tableSRC, nonnormal = c("age", "cross_less_rt"), quote = FALSE, noSpaces = TRUE, digits=NULL)
write.csv(tableSRC, "output/TableSRC_pCR.csv", row.names = TRUE, na = "")

# All relevant pathology characteristics in patients who underwent resection (TRG-1-2-3-4) split by pCR
# This is exploratory only
pathVars <- c("biopsies_mucin", "biopsies_src", "biopsies_src_atypical", "biopsies_rfall", "biopsies_diffgrade",
              "biopsies_mucin_any", "biopsies_src_any", "biopsies_src_atypical_any", "biopsies_rfall_any", 
              "TRG", "ypT", "ypn", "pCR", "prepT", "prepn", "radicality")
tablePA_res <- CreateTableOne(vars = pathVars, factorVars = pathVars, strata = "pCR", data = modDataRes)
tablePA_res <- print(tablePA_res, quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tablePA_res, "output/TablePAResectionTRGall_pCR.csv", row.names = TRUE, na = "")

# All relevant pathology characteristics in patients who underwent resection (TRG-1-2-3-4) split by TRG1-2 vs TRG3-4
# This is exploratory only
pathVars <- c("biopsies_mucin", "biopsies_src", "biopsies_src_atypical", "biopsies_rfall", "biopsies_diff",
              "biopsies_mucin_any", "biopsies_src_any", "biopsies_src_atypical_any", "biopsies_rfall_any", 
              "TRG", "ypT", "ypn", "pCR", "prepT", "prepn", "radicality")
tablePA_res <- CreateTableOne(vars = pathVars, factorVars = pathVars, strata = "TRG_cat", data = modDataRes)
tablePA_res <- print(tablePA_res, quote = FALSE, noSpaces = TRUE, contDigits=1)
write.csv(tablePA_res, "output/TablePAResectionTRGall_TRG12vs34.csv", row.names = TRUE, na = "")


######### Comparison pre-treatment and post-treatment PA factors #########
rm(edges)

#########  Do patients change category? For patients with resection and residual tumor TRG 2-3-4 ######### 
edges <- modDataResResidual[, c("biopsies_mucin", "res_mucin")]
# this gives the count of all group changes
edges <- edges %>%
  group_by(biopsies_mucin, res_mucin) %>%
  summarise(Value = n())

edges <- modDataResResidual[, c("biopsies_diffgrade", "res_diffgrade")]
edges <- edges %>%
  group_by(biopsies_diffgrade, res_diffgrade) %>%
  summarise(Value = n())

# in all resected patients, including TRG1
edges <- modDataRes[, c("biopsies_src", "res_src")]
edges <- edges %>%
  group_by(biopsies_src, res_src) %>%
  summarise(Value = n())

edges <- modDataRes[, c("biopsies_src_atypical", "res_src_atypical")]
table(modDataRes$TRG, modDataRes$res_src_atypical) # check if all patients with TRG1 have 0% poorly cohesive cells
edges <- edges %>%
  group_by(biopsies_src_atypical, res_src_atypical) %>%
  summarise(Value = n())

edges <- modDataRes[, c("biopsies_srcpcc", "res_srcpcc")]
table(modDataRes$TRG, modDataRes$res_srcpcc) # check if all patients with TRG1 have 0% poorly cohesive cells
edges <- edges %>%
  group_by(biopsies_srcpcc, res_srcpcc) %>%
  summarise(Value = n())

# prepare "edges" for the plot
names(edges) <- c("N1", "N2", "Value")
edges <- as.data.frame(edges)
edges$N1 <- as.character(edges$N1)
edges$N2 <- as.character(edges$N2)
categories <- unique(c(edges$N1, edges$N2))

# for 4 categories risk factors
nodes <- data.frame(ID = outer(categories, c("_start", "_end"), FUN = "paste0")[1:8],
                    x = c(1, 1, 1, 1, 2, 2, 2, 2), 
                    y = c(1, 2, 3, 4, 1, 2, 3, 4))
edges$N1 <- paste0(edges$N1, "_start") # only do once
edges$N2 <- paste0(edges$N2, "_end") # only do once

## for differentiation grade, 3 categories instead of 4
nodes <- data.frame(ID = outer(categories, c("_start", "_end"), FUN = "paste0")[1:4],
                    x = c(1, 1, 2, 2), 
                    y = c(1, 2, 1, 2))
edges$N1 <- paste0(edges$N1, "_start") # only do once
edges$N2 <- paste0(edges$N2, "_end") # only do once

#palette = paste0(brewer.pal(4, "Set1"), "90") other palette
styles <- lapply(nodes$y, function(n) {
  list(col = colors[n], lty = 0, textcol = "black") # if you use n+1 you can loop through colors
  })
names(styles) <- nodes$ID

r <- makeRiver(nodes=nodes, edges=edges, 
               styles = styles, node_labels = rep(categories, 2))

# customized labels
r <- makeRiver(nodes=nodes, edges=edges, 
               styles = styles, node_labels = c("good-moderate", "poor"))

png(filename = "figures/FigRPmucin.png", units = "cm", width=20, height=30, res=300) # keep these sizes
plot(r)
dev.off()

png(filename = "figures/FigRPsrc.png", units = "cm", width=20, height=30, res=300) # keep these sizes
plot(r)
dev.off()

png(filename = "figures/FigRPpcc.png", units = "cm", width=20, height=30, res=300) # keep these sizes
plot(r)
dev.off()

png(filename = "figures/FigRPsrc_plus_pcc.png", units = "cm", width=20, height=30, res=300) # keep these sizes
plot(r)
dev.off()

png(filename = "figures/FigRPdiffgrade.png", units = "cm", width=20, height=30, res=300) # keep these sizes
plot(r)
dev.off()

######### Primary endpoint: TRG 1-2 vs 3-4 including 2 patients with T4b, assigned -> TRG4 ######### 
modData <- subset(data, data$TRG_T4b_cat != 'NA')

# reference model
fit <- glm(TRG_T4b_cat ~ biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
# returns OR and 95%CI by hand
exp(summary(fit)$coefficients["biopsies_diffgrade1",1] + 
        +     qnorm(c(0.5, 0.025 ,0.975)) * summary(fit)$coefficients["biopsies_diffgrade1",2])
exp(-coef(fit))
summary(fit)

# is pre-treatment differentiation grade associated with pCR? pCR defined as "1"
fit <- glm(pCR ~ biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modDataRes, family = binomial())
exp(summary(fit)$coefficients["biopsies_diffgrade1",1] + 
      +     qnorm(c(0.5, 0.025 ,0.975)) * summary(fit)$coefficients["biopsies_diffgrade1",2])
summary(fit)
table(modData$biopsies_diffgrade, modData$pCR)

# are SRC associated with pCR (validation Corsini et al)
table(modDataRes$biopsies_src_any, modDataRes$pCR)

# PA factors
# extracellular mucin in biopsy
fit_null <- glm(TRG_T4b_cat ~ biopsies_mucin_any + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_mucin_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

fit_null <- glm(TRG_T4b_cat ~ biopsies_mucin_010 + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_mucin_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

# src
fit_null <- glm(TRG_T4b_cat ~ biopsies_src_any + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_src_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

fit_null <- glm(TRG_T4b_cat ~ biopsies_src_010 + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_src_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

# poorly cohesive cells
fit_null <- glm(TRG_T4b_cat ~ biopsies_src_atypical_any + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_src_atypical_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

fit_null <- glm(TRG_T4b_cat ~ biopsies_src_atypical_010 + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_src_atypical_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

# signet-ring cells + poorly cohesive cells
fit_null <- glm(TRG_T4b_cat ~ biopsies_srcpcc_any + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_srcpcc_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

fit_null <- glm(TRG_T4b_cat ~ biopsies_srcpcc_010 + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_srcpcc_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

# rfall
fit_null <- glm(TRG_T4b_cat ~  biopsies_rfall_any + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_rfall_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

fit_null <- glm(TRG_T4b_cat ~  biopsies_rfall_010 + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
fit_alt <- glm(TRG_T4b_cat ~ biopsies_rfall_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = modData, family = binomial())
makeLRSum(fit_null, fit_alt)

######### Survival analyses ######### 
# Check survival times
# Plot distribution OS and DFS in months
hist(data$tvar, breaks = 20)
hist(data$dfsvar, breaks = 20)

# Overall survival
summary(data[data$status == 0,]$tvar) # min, median follow-up
KM_fit <- survfit(Surv(tvar, status) ~ 1, data=data) # ~1 defines all categories, i.e. survival curve for whole dataset
KM_fit # gives median survival, i.e. the time by which half of the subjects will experience the event
summary(KM_fit, times = c(12, 24, 36, 48, 60)) # survival probabilities at several follow-up times in months
quantile(KM_fit, probs = 1 - c(0.25, 0.75)) # at how many days the survival probability is 25% and 75%
plot(KM_fit, xlab = "Time to death (months)", ylab = "Survival probability", 
     main = "Kaplan-Meier Estimate of S(t) for all patients (n = 325)")

# Disease-free survival
summary(data[data$dfstatus == 0,]$dfsvar) # min, median follow-up
KM_fit <- survfit(Surv(dfsvar, dfstatus) ~ 1, data=data) # ~1 defines all categories, i.e. survival curve for whole dataset
KM_fit # gives median survival, i.e. the time by which half of the subjects will experience the event
summary(KM_fit, times = c(12, 24, 36, 48, 60)) # survival probabilities at several follow-up times in months
quantile(KM_fit, probs = 1 - c(0.25, 0.75)) # at how many days the survival probability is 25% and 75%
plot(KM_fit, xlab = "Time to death (months)", ylab = "Survival probability", 
     main = "Kaplan-Meier Estimate of S(t) for all patients (n = 325)")

#------------------ KM plots ---------------
## Check assumptions whether the proportional hazards (PH) function is satisfied
# Pre-treatment PA factors
logrank <- survdiff(Surv(tvar, status == 1) ~ biopsies_diffgrade, data = data)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ biopsies_diffgrade, data = data)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ biopsies_mucin_any, data = data)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ biopsies_mucin_any, data = data)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ biopsies_src_any, data = data)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ biopsies_src_any, data = data)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ biopsies_src_atypical_any, data = data)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ biopsies_src_atypical_any, data = data)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ biopsies_srcpcc_any, data = data)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ biopsies_srcpcc_any, data = data)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ biopsies_rfall_any, data = data)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ biopsies_rfall_any, data = data)
plot(KM_fit, fun = "cumhaz")

# Post-treatment PA factors (only for patients with residual tumor in the resection specimen)
logrank <- survdiff(Surv(tvar, status == 1) ~ res_diffgrade, data = modDataResResidual)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ res_diffgrade, data = modDataResResidual)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ res_mucin_any, data = modDataResResidual)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ res_mucin_any, data = modDataResResidual)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ res_src_any, data = modDataResResidual)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ res_src_any, data = modDataResResidual)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ res_src_atypical_any, data = modDataResResidual)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ res_src_atypical_any, data = modDataResResidual)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ res_srcpcc_any, data = modDataResResidual)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ res_srcpcc_any, data = modDataResResidual)
plot(KM_fit, fun = "cumhaz")

logrank <- survdiff(Surv(tvar, status == 1) ~ res_rfall_any, data = modDataResResidual)
logrank
KM_fit <- survfit(Surv(tvar, status == 1) ~ res_rfall_any, data = modDataResResidual)
plot(KM_fit, fun = "cumhaz")

## Choose patient set before making plots
# Extra datasets
# IDs with very atypical tumor cells in patients with TRG 2-3-4 residual tumor (explorative only, for EMC patients)
dataAtyp <- subset(data, data$Institute.Abbreviation != "RAD")
dataAtyp <- subset(dataAtyp, dataAtyp$resection == "yes")
dataAtyp <- subset(dataAtyp, dataAtyp$TRG != "TRG1")
ids <- c("EMC004", "EMC011", "EMC012", "EMC014", "EMC015", "EMC017", "EMC021", "EMC024", "EMC025", "EMC027", "EMC031", "EMC032", "EMC040", "EMC041", "EMC043", "EMC047", "EMC048", "EMC051", "EMC056", "EMC062", "EMC063", "EMC068", "EMC069", "EMC070", "EMC071", "EMC079", "EMC080", "EMC084", "EMC092", "EMC096", "EMC101", "EMC105", "EMC106", "EMC108", "EMC111", "EMC114", "EMC117", "EMC121", "EMC123", "EMC124", "EMC126", "EMC129", "EMC133", "EMC159", "EMC165", "EMC167", "EMC169", "EMC184", "EMC191", "PATH004")
dataAtyp$atypical <- as.factor(ifelse(dataAtyp$Record.Id %in% ids, "yes", "no"))

# IDs with Barrett segment and SRC in resection specimen
dataBarrett <- subset(modDataRes, modDataRes$res_src_any == "1")
ids2 <- c("PATH011", "PATH018", "EMC078", "EMC130", "EMC155", "RAD104")
dataBarrett$Barrett <- as.factor(ifelse(dataBarrett$Record.Id %in% ids2, "yes", "no"))

# subdata for change in N status (plus TRG status)
subdata <- subset(data, (data$resection == "yes")) # only Radboud UMC
subdata <- subset(subdata, (subdata$Institute.Abbreviation == "RAD"))
subdata <- subdata[!(is.na(subdata$prepN_cat) | is.na(subdata$ypN_cat)), ]
subdata$changeN <- as.factor(paste(subdata$prepN_cat, subdata$ypN_cat))

# For plot with TRG and change in N status
modDataResTRG12 <- subset(modDataRes, modDataRes$TRG_cat == "0") # necessary for plot using subdata2
subdata2 <- modDataResTRG12[!(is.na(modDataResTRG12$prepN_cat) | is.na(modDataResTRG12$ypN_cat)), ]
subdata2$TRG_prepN_cat <- as.factor(paste(subdata2$TRG, subdata2$prepN_cat))
subdata2$TRG_changeN <- as.factor(paste(subdata2$TRG, subdata2$prepN_cat, subdata2$ypN_cat))

subdata3 <- subset(data, (data$biopsies_diffgrade == "1")) # only G3 tumors

## Make some plots separately
# use when saving a plot
#png(filename = "figures/FigKM_OS_ .png", units = "cm", width=25, height=15, res=300) # keep these sizes

#Explorative plots
# Patients with SRCs in resection specimen: with and without residual Barrett
fit <- makeKMplot(data = dataBarrett, 
                  timevar = 'tvar', 
                  eventvar = 'status', 
                  groupingvar = 'Barrett', 
                  titl = 'Overall survival', 
                  gvartitle = 'Signet-ring cells: patients with Barrett vs no Barrett',
                  #grouplabels = c("<1%", "≥1%") # c("absent/present", "present/absent", "present/present", "absent/absent")
) #for grouplabels choose NA if you want the function to use levels(groupingvar)
#dev.off()

# Combination pCR and differentiation grade
# Drop unused levels
modDataRes$pCR_diffgrade <- droplevels(modDataRes$pCR_diffgrade)
fit <- makeKMplot(data = modDataRes, 
                  timevar = 'tvar', 
                  eventvar = 'status', 
                  groupingvar = 'pCR_diffgrade', 
                  titl = 'Overall survival', 
                  gvartitle = 'pCR and differentiation grade',
                  #grouplabels = c("<1%", "≥1%") # c("absent/present", "present/absent", "present/present", "absent/absent")
                  ) #for grouplabels choose NA if you want the function to use levels(groupingvar)
#dev.off()

# SRC in biopsies only for patients undergoing surgery
fit <- makeKMplot(data = modDataRes, 
                  timevar = 'tvar', 
                  eventvar = 'status', 
                  groupingvar = 'biopsies_src_any', 
                  titl = 'Overall survival', 
                  gvartitle = 'signet-ring cells in biopsies (pt ondergoing surgery)',
                  #grouplabels = c("<1%", "≥1%") # c("absent/present", "present/absent", "present/present", "absent/absent")
) #for grouplabels choose NA if you want the function to use levels(groupingvar)
#dev.off()

# Key Image: signet-ring cells pre/post-treatment
# Supplemental Figure 2
png(filename = "figures/FigKM_OS_SRCs_pre_vs_post.png", units = "cm", width=25, height=15, res=1200) # keep these sizes
fit <- makeKMplot(data = modDataRes, 
                  timevar = 'tvar', 
                  eventvar = 'status', 
                  groupingvar = 'src_categories',
                  titl = 'Overall survival in patients undergoing resection (n = 284)', 
                  gvartitle = 'signet-ring cells pre- and post-treatment'
) #for grouplabels choose NA if you want the function to use levels(groupingvar)
dev.off()

# Key Image: poorly cohesive cells cells pre/post-treatment
fit <- makeKMplot(data = modDataRes, 
                  timevar = 'tvar', 
                  eventvar = 'status', 
                  groupingvar = 'pcc_categories',
                  titl = 'Overall survival in patients undergoing resection (n = 284)', 
                  gvartitle = 'poorly cohesive cells pre- and post-treatment'
) #for grouplabels choose NA if you want the function to use levels(groupingvar)

# Key Image: any signet-ring cells + poorly cohesive cells pre/post-treatment
fit <- makeKMplot(data = modDataRes, 
                  timevar = 'tvar', 
                  eventvar = 'status', 
                  groupingvar = 'srcpcc_categories', 
                  titl = 'Overall survival in patients undergoing resection (n = 284)', 
                  gvartitle = 'signet-ring cells + poorly cohesive cells pre- and post-treatment'
) #for grouplabels choose NA if you want the function to use levels(groupingvar)

# Make all plots at once
gvars$grouplabels <- sapply(gvars$grouplabels, gsub, pattern = ">0%", replacement = "≥1%")
fit_list_os <- list()
fit_list_dfs <- list()
pvalue_os <- c()
pvalue_dfs <- c()
for (i in 1:nrow(gvars)) {
  print(i)
  gvar <- gvars[i,]
  print(gvar)
  if (gvar$data == 'data') {
    plotdata <- data
  } else if (gvar$data == 'modDataRes') {
    plotdata <- modDataRes
  } else if (gvar$data == 'modDataResResidual') {
    plotdata <- modDataResResidual
  } else if (gvar$data == 'subdata') {
    plotdata <- subdata
  } else if (gvar$data == 'subdata2') {
    plotdata <- subdata2
  } else if (gvar$data == 'dataAtyp') {
    plotdata <- dataAtyp
  } else {
    stop('Specified plot data not found in memory.')
  }
  
  result <- makeKMplot(data = plotdata, 
                                              timevar = 'tvar', 
                                              eventvar = 'status', 
                                              groupingvar = gvar$groupingvar, 
                                              titl = 'Overall survival', 
                                              gvartitle = gvar$gvartitle,
                                              grouplabels = eval(parse(text=gvar$grouplabels)),
                                              save = TRUE,
                                              filename = paste0('figures/', 'FigKM_OS_', gvar$groupingvar, '.png'))
  
  fit_list_os[[gvar$groupingvar]] <- result$sfit
  pvalue_os <- c(pvalue_os, result$pvalue)
  
  result <- makeKMplot(data = plotdata, 
                                              timevar = 'dfsvar', 
                                              eventvar = 'dfstatus', 
                                              groupingvar = gvar$groupingvar, 
                                              titl = 'Disease-free survival', 
                                              gvartitle = gvar$gvartitle,
                                              grouplabels = eval(parse(text=gvar$grouplabels)),
                                              save = TRUE,
                                              filename = paste0('figures/', 'FigKM_DFS_', gvar$groupingvar, '.png'))
  fit_list_dfs[[gvar$groupingvar]] <- result$sfit
  pvalue_dfs <- c(pvalue_dfs, result$pvalue)
}

# Compute OS metrics for all variables
medianOS <- surv_median(fit_list_os, combine=TRUE)
colnames(medianOS)[3:5] <- paste(colnames(medianOS)[3:5],"OS",sep="_")
medianOS$median_OS <- round(medianOS$median_OS, digits = 1)
medianOS$lower_OS <- round(medianOS$lower_OS, digits = 1)
medianOS$upper_OS <- round(medianOS$upper_OS, digits = 1)
medianOS$medianCI_OS <- paste0(medianOS$median_OS, " (", medianOS$lower_OS, " - ", medianOS$upper_OS, ")")
temp <- as.data.frame(pvalue_os)
temp$pvalue_os <- round(temp$pvalue_os, digits = 3)
temp$id <- gvars$groupingvar
metricsOS <- merge(medianOS, temp)

# Compute DFS metrics for all variables
medianDFS <- surv_median(fit_list_dfs, combine=TRUE)
colnames(medianDFS)[3:5] <- paste(colnames(medianDFS)[3:5],"DFS",sep="_")
medianDFS$median_DFS <- round(medianDFS$median_DFS, digits = 1)
medianDFS$lower_DFS <- round(medianDFS$lower_DFS, digits = 1)
medianDFS$upper_DFS <- round(medianDFS$upper_DFS, digits = 1)
medianDFS$medianCI_DFS <- paste0(medianDFS$median_DFS, " (", medianDFS$lower_DFS, " - ", medianDFS$upper_DFS, ")")
temp <- as.data.frame(pvalue_dfs)
temp$pvalue_dfs <- round(temp$pvalue_dfs, digits = 3)
temp$id <- gvars$groupingvar
metricsDFS <- merge(medianDFS, temp)

metricsOS_DFS <- merge(metricsOS, metricsDFS)
write.csv(metricsOS_DFS, "output/TableMedian_OS_DFS.csv", row.names = TRUE, na = "")

#------------------ Cox regression - models with PA factors ---------------------------- 
## Choose outcome
#Cox model for overall survival
sdat <- data
tvar <- sdat$tvar  #time variable
evar <- sdat$status #event variable

#Cox model for disease-free survival
sdat <- data
tvar <- sdat$dfsvar  #time variable
evar <- sdat$dfstatus #event variable

#### Pre-treatment Cox models ####
# Reference model with differentiation grade
fit_null <- coxph(Surv(tvar, evar) ~ biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) 
#fit_alt <- coxph(Surv(tvar, evar) ~ age + sex + prepT_cat + prepN_cat + biopsies_diffgrade, data = sdat) 
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_diffgrade*sex + age + cT_cat + cN_cat, data = sdat) # interaction term

fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_diffgrade + ns(age, 3) + sex + cT_cat + cN_cat, data = sdat) # non-linear effect age
summary(fit_null)
summary(fit_alt)
anova(fit_null, fit_alt) # interaction terms and non-linear terms are NS, proceed with linear model. NB anova test only applicable to nested models

makeCoxSum(fit_null, fit_alt)

# Check proportional hazards assumption (variables and GLOBAL test should not be significant)
cox.zph(fit_null)
plot(cox.zph(fit_null))

# Extracellular mucine, test for the effect of differentiation grade
fit_null <- coxph(Surv(tvar, evar) ~ biopsies_mucin_any + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_mucin_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade

fit_null <- coxph(Surv(tvar, evar) ~ biopsies_mucin_010 + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_mucin_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
cox.zph(fit_alt)

# SRC, test for the effect of differentiation grade
fit_null <- coxph(Surv(tvar, evar) ~ biopsies_src_any + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_src_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade

fit_null <- coxph(Surv(tvar, evar) ~ biopsies_src_010 + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_src_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
cox.zph(fit_alt)

# poorly cohesive cells, test for the effect of differentiation grade
fit_null <- coxph(Surv(tvar, evar) ~ biopsies_src_atypical_any + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_src_atypical_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade

fit_null <- coxph(Surv(tvar, evar) ~ biopsies_src_atypical_010 + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_src_atypical_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade
summary(fit_null)
summary(fit_alt)

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
cox.zph(fit_alt)

# signet-ring cells + poorly cohesive cells, test for the effect of differentiation grade
fit_null <- coxph(Surv(tvar, evar) ~ biopsies_srcpcc_any + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_srcpcc_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade

fit_null <- coxph(Surv(tvar, evar) ~ biopsies_srcpcc_010 + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_srcpcc_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade
summary(fit_null)
summary(fit_alt)

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
cox.zph(fit_alt)

# highest of all risk factors
fit_null <- coxph(Surv(tvar, evar) ~ biopsies_rfall_any + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_rfall_any + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade

fit_null <- coxph(Surv(tvar, evar) ~ biopsies_rfall_010 + age + sex + cT_cat + cN_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_rfall_010 + biopsies_diffgrade + age + sex + cT_cat + cN_cat, data = sdat) # with diff grade
summary(fit_null)
summary(fit_alt)

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
cox.zph(fit_alt)

#### Post-treatment Cox models ####
# Choose data OS or DFS
sdat <- modDataResResidual
tvar <- sdat$tvar  #time variable
evar <- sdat$status #event variable

tvar <- sdat$dfsvar  #time variable
evar <- sdat$dfstatus #event variable

# Make reference post-treatment model
# one observation deleted = patient in whom LN slides were not available and prep-N is thus missing
fit_null <- coxph(Surv(tvar, evar) ~ biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat)
summary(fit_null)

fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_diffgrade*sex + age + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) #interaction term
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_diffgrade + age + sex + prepT_cat + prepN_cat*ypN_cat + ypT_cat + TRG_cat + radicality, data = sdat) #another interaction term
fit_alt <- coxph(Surv(tvar, evar) ~ biopsies_diffgrade + ns(age, 3) + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) #splines for age
anova(fit_null, fit_alt) # interaction term models and model with splines do not improve the fit, proceed with linear model

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_null))
cox.zph(fit_null)

# extracellular mucin
fit_null <- coxph(Surv(tvar, evar) ~ res_mucin_any + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_mucin_any + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # with diff grade

fit_null <- coxph(Surv(tvar, evar) ~ res_mucin_010 + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_mucin_010 + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # with diff grade

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
cox.zph(fit_alt)

# signet-ring cells
fit_null <- coxph(Surv(tvar, evar) ~ res_src_any + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_src_any + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # with diff grade
summary(fit_null)
summary(fit_alt)

fit_null <- coxph(Surv(tvar, evar) ~ res_src_010 + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_src_010 + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # with diff grade
summary(fit_null)
summary(fit_alt)

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
cox.zph(fit_alt)

# poorly cohesive cells
fit_null <- coxph(Surv(tvar, evar) ~ res_src_atypical_any + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_src_atypical_any + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # with diff grade

fit_null <- coxph(Surv(tvar, evar) ~ res_src_atypical_010 + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_src_atypical_010 + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # with diff grade

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
cox.zph(fit_alt)

# signet-ring cells + poorly cohesive cells
fit_null <- coxph(Surv(tvar, evar) ~ res_srcpcc_any + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_srcpcc_any + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # with diff grade

fit_null <- coxph(Surv(tvar, evar) ~ res_srcpcc_010 + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_srcpcc_010 + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat + radicality, data = sdat) # with diff grade

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
cox.zph(fit_alt)

# highest of all risk factors
fit_null <- coxph(Surv(tvar, evar) ~ res_rfall_any + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_rfall_any + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat, data = sdat) # with diff grade

fit_null <- coxph(Surv(tvar, evar) ~ res_rfall_010 + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ res_rfall_010 + biopsies_diffgrade + age + sex + prepT_cat + prepN_cat + ypT_cat + ypN_cat + TRG_cat, data = sdat) # with diff grade

makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))

######### exploratory analysis (extra)
# very atypical cells (M. Doukas), pt with TRG 2-3-4 tumor
dataAtyp <- subset(modDataResResidual, data$Institute.Abbreviation != "RAD")
ids <- c("EMC004", "EMC011", "EMC012", "EMC014", "EMC015", "EMC017", "EMC021", "EMC024", "EMC025", "EMC027", "EMC031", "EMC032", "EMC040", "EMC041", "EMC043", "EMC047", "EMC048", "EMC051", "EMC056", "EMC062", "EMC063", "EMC068", "EMC069", "EMC070", "EMC071", "EMC079", "EMC080", "EMC084", "EMC092", "EMC096", "EMC101", "EMC105", "EMC106", "EMC108", "EMC111", "EMC114", "EMC117", "EMC121", "EMC123", "EMC124", "EMC126", "EMC129", "EMC133", "EMC159", "EMC165", "EMC167", "EMC169", "EMC184", "EMC191", "PATH004")
dataAtyp$atypical <- as.factor(ifelse(dataAtyp$Record.Id %in% ids, "yes", "no"))

sdat <- dataAtyp
tvar <- sdat$tvar  #time variable
evar <- sdat$status #event variable

fit_null <- coxph(Surv(tvar, evar) ~ age + sex + cT_cat + cN_cat + atypical, data = sdat) # without diff grade
fit_alt <- coxph(Surv(tvar, evar) ~ age + sex + cT_cat + cN_cat + atypical + biopsies_diffgrade, data = sdat) # with diff grade
makeCoxSum(fit_null, fit_alt)
plot(cox.zph(fit_alt))
summary(fit_null)
summary(fit_alt)
anova(fit_null, fit_alt)

