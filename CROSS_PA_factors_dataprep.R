# Script analysis CROSS PA factors
# M.J. Valkema, November 2022
library(dplyr) # function coalesce
library(lubridate) # date difference

setwd(dir=getwd())
# Import dataset with numbers instead of labels
setwd("/Users/Maartje/repos/CROSS-PA-factors")
data <- read.csv(file="data/CROSS_PA_factors_export_20220715.csv",
                          header=TRUE, sep = ";", na.strings = "")

# Dates
data$date_diagnosis <- as.Date(data$date_diagnosis, "%d-%m-%Y")
data$end_ncrt <- as.Date(data$end_ncrt, "%d-%m-%Y")
data$date_death <- as.Date(data$date_death, "%d-%m-%Y")
data$date_last_fu <- as.Date(data$date_last_fu, "%d-%m-%Y")
data$date_recurrence <- as.Date(data$date_recurrence, "%d-%m-%Y")

# Change colnames
names(data)[names(data) == "biopsies_mucous"] <- "biopsies_mucin"
names(data)[names(data) == "biopsies_mucous_other"] <- "biopsies_mucin_other"
names(data)[names(data) == "res_mucous"] <- "res_mucin"
names(data)[names(data) == "res_mucous_other"] <- "res_mucin_other"


# Re-code variables
library(car)
data$sex <- as.factor(car::recode(data$sex, "c('1')='male'; c('2')='female'"))
data$cT <- car::recode(data$cT_stage, "c('3')='cT3'; 
                   c('2')='cT2'; 
                   c('1b')='cT1b'; 
                   c('4a', '4')='cT4a';
                   c('99')='cTx'")
data$cT_cat <- as.factor(car::recode(data$cT_stage, "c('4a', '4', '3', '99')='1'; 
                   c('1b', '2')='0'")) #this was cT3-4a and cT1b-2
data$ypT <- car::recode(data$ypt, "c('1', '1a', '1b')='ypT1';
                   c('2')='ypT2';
                   c('3')='ypT3';
                   c('4')='ypT4'")

data$ypT_cat <- as.factor(car::recode(data$ypt, "c('1', '1a', '1b', '2')='0'; c('3', '4')='1'"))

data$prepT <- car::recode(data$prept, "c('1', '1b')='prepT1'; 
                   c('2') = 'prepT2';
                   c('3')='prepT3';
                   c('4', '4a')='prepT4'")

data$prepT_cat <- as.factor(car::recode(data$prepT, "c('prepT1', 'prepT2')='0';
                   c('prepT3', 'prepT4') = '1'")) 

data$cN <- car::recode(data$cN_stage, "c('0')='cN0'; 
                   c('1')='cN1'; 
                   c('2')='cN2'; 
                   c('3')='cN3'")
data$cN_cat <- as.factor(car::recode(data$cN_stage, "c('0')='0'; c('1', '2', '3')='1'")) # cN0 and cN+

data$ypN_cat <- as.factor(car::recode(data$ypn, "c('ypN0')='0'; c('ypN1', 'ypN2', 'ypN3')='1'"))

data$prepN_cat <- as.factor(car::recode(data$prepn, "c('prepN0')='0'; c('prepN1', 'prepN2', 'prepN3') = '1'"))

data$changeN <- as.factor(paste(data$prepN_cat, data$ypN_cat))
data$changeN <- as.factor(car::recode(data$changeN, "c('0 0')='0'; c('1 0')='1'; c('0 1', '1 1', 'NA 1')='2'; c('NA 0') = NA"))

data$cM <- car::recode(data$cM_stage, "c('0')='cM0'; c('1')='cM1'")

data$cross_all <- car::recode(data$cross_all, "c('0')='reduction'; c('1')='yes'; c('99')='unknown'")

data$histology <- car::recode(data$histology, "c('2')='adenocarcinoma'")

data$location <- car::recode(data$tumor_location, "c('1')='proximal'; 
                   c('2')='mid'; 
                   c('3', '4')='distal/GEJ'")

data$hospital <- car::recode(data$Institute.Abbreviation, "c('EMC', 'PATH')='Erasmus MC' ; c('RAD') = 'Radboud UMC'")

data$surgery <- car::recode(data$surgery, "c('0')='no'; c('1')='yes'")
data$resection <- car::recode(data$resection, "c('0')='no'; c('1')='yes'")
data$interval_ncrt_surgery <- data$interval_ncrt_surgery / 7 # to get weeks

data$specify_no_surgery <- car::recode(data$specify_no_surgery, "c('1')='active surveillance'; 
                   c('2')='refusal'; 
                   c('3')='poor condition'; 
                   c('4')='distant metastases';
                   c('88')='other';")

data$spec_no_surgery_other <- car::recode(data$specify_no_surgery_other, "c('died', 'died after CRE-2', 'died before response evaluation')='died'; 
                   c('inoperable due to submucosal metastases in the esophagus', 'solitary thyroid metastasis + SANO avant la lettre') = 'metastatic disease'")

data$spec_no_res <- car::recode(data$specify_no_resection, "c('liver metastases', 'peroperative liver metastasis', 'peroperative metastases')='peroperative metastases'; 
                   c('T4b', 'T4b (ingrowth aorta)') = 'unresectable tumor (T4b)'")

data$TRG <- as.factor(car::recode(data$trg, "c('1')='TRG1'; 
                   c('2')='TRG2'; 
                   c('3')='TRG3'; 
                   c('4')='TRG4'"))
data$TRG_cat <- as.factor(car::recode(data$trg, "c('1', '2')='0'; c('3', '4')='1'"))
data$trg <- as.factor(data$trg)

data$TRG_T4b <- ifelse(data$spec_no_res == 'unresectable tumor (T4b)', 'T4b', NA)
data$TRG_T4b <- coalesce(data$TRG_T4b, data$TRG)
data$TRG_T4b <- as.factor(car::recode(data$TRG_T4b, "c('T4b')='TRG4'"))
data$TRG_T4b_cat <- as.factor(car::recode(data$TRG_T4b, "c('TRG1', 'TRG2')='0'; c('TRG3', 'TRG4')='1'")) # this was TRG1-2 and TRG3-4

data$TRG_ypTcat <- as.factor(paste(data$TRG, data$ypT_cat))
data$TRG_ypTcat <- as.factor(car::recode(data$TRG_ypTcat, "c('NA NA', 'TRG1 NA')='NA'"))

data$residualtumor <- car::recode(data$trg, "c('1')='no'; c('2', '3', '4')='yes'")

data$pCR <- ifelse((data$TRG_T4b == "TRG1" & data$ypn == "ypN0"), "1", "0") # define patiens with pathologically complete response
data$pCR <- as.factor(data$pCR)

data$Rmargin <- data$radicality
data$radicality <- car::recode(data$radicality, "c('0')='R0'; 
                   c('1', '2', '3')='R1'") # 1 is proximal resection margin, 2 is distal resection margin, 3 is circumferential resection margin

# Re-code scores for biopsies
data$biopsies_src_atypical <- car::recode(data$biopsies_src_atypical, "c('88')='0'") ## There is one with score 88, but this score was assigned a 0

## Risk factor data preparation
# Separate risk factors

# Comparison before and after treatment in patients who underwent resection with numeric variables
data$biopsies_mucin_num <- data$biopsies_mucin
data$biopsies_src_num <- data$biopsies_src
data$biopsies_src_atypical_num <- data$biopsies_src_atypical
data$res_mucin_num <- data$res_mucin
data$res_src_num <- data$res_src
data$res_src_atypical_num <- data$res_src_atypical

# Sum of signet-ring cells and poorly cohesive cells
# Re-code categories into percentiles
data$biopsies_src_perc <- car::recode(data$biopsies_src_num, "c('0')='0'; c('1') = '0.1'; c('2') = '0.4'; c('3')='0.5'")
data$biopsies_src_atypical_perc <- car::recode(data$biopsies_src_atypical_num, "c('0')='0'; c('1') = '0.1'; c('2') = '0.4'; c('3')='0.5'")
data$biopsies_srcpcc_perc <- data$biopsies_src_perc + data$biopsies_src_atypical_perc
data$biopsies_srcpcc_sum <- as.factor(car::recode(data$biopsies_srcpcc_perc, "c('0')='0'; c('0.1') = '1'; c('0.2', '0.4') = '2'; c('0.5', '0.6', '0.8', '0.9')='3'"))

# We re-examined the cases in which summation is not evident: 1-10% + 1-10%, 1-10% + 11-50%;, 11-50% + 11-50%
data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "EMC024"], 3)
data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "EMC036"], 1)
data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "EMC169"], 1)

data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "PATH003"], 2)
data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "PATH012"], 2)

data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "RAD004"], 3)
data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "RAD023"], 2)
data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "RAD034"], 3)
data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "RAD054"], 2)
data$biopsies_srcpcc_sum <- replace(data$biopsies_srcpcc_sum, data$biopsies_srcpcc_sum[data$Record.Id == "RAD074"], 1)

data$biopsies_srcpcc_any <- as.factor(car::recode(data$biopsies_srcpcc_sum, "c('0')='0'; c('1', '2', '3')='1'"))
data$biopsies_srcpcc_010 <- as.factor(car::recode(data$biopsies_srcpcc_sum, "c('0', '1')='0'; c('2', '3')='1'"))
data$biopsies_srcpcc <- as.factor(car::recode(data$biopsies_srcpcc_sum, "c('0')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'"))

# Define overall risk factor based on the max in each of the categories for every patient
data$biopsies_rfall <- apply(data[, c('biopsies_mucin', 'biopsies_src', 'biopsies_src_atypical')], MARGIN=1, max)
data$biopsies_rfall_any <- as.factor(car::recode(data$biopsies_rfall, "c('0')='0'; c('1', '2', '3')='1'"))
data$biopsies_rfall_010 <- as.factor(car::recode(data$biopsies_rfall, "c('0', '1')='0'; c('2', '3')='1'"))
data$biopsies_rfall <- as.factor(car::recode(data$biopsies_rfall, "c('0')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'"))

data$biopsies_mucin_any <- as.factor(car::recode(data$biopsies_mucin, "c('0')='0'; c('1', '2', '3')='1'"))
data$biopsies_mucin_010 <- as.factor(car::recode(data$biopsies_mucin, "c('0', '1')='0'; c('2', '3')='1'"))
data$biopsies_mucin <- as.factor(car::recode(data$biopsies_mucin, "c('0')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'"))

data$biopsies_src_any <- as.factor(car::recode(data$biopsies_src, "c('0')='0'; c('1', '2', '3')='1'"))
data$biopsies_src_010 <- as.factor(car::recode(data$biopsies_src, "c('0', '1')='0'; c('2', '3')='1'"))
data$biopsies_src <- as.factor(car::recode(data$biopsies_src, "c('0')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'"))

data$biopsies_src_atypical_any <- as.factor(car::recode(data$biopsies_src_atypical, "c('0')='0'; c('1', '2', '3')='1'"))
data$biopsies_src_atypical_010 <- as.factor(car::recode(data$biopsies_src_atypical, "c('0', '1')='0'; c('2', '3')='1'"))
data$biopsies_src_atypical <- as.factor(car::recode(data$biopsies_src_atypical, "c('0')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'"))

# Differentiation grade biopsies
data$biopsies_diff <- as.factor(car::recode(data$biopsies_diff, "c('1')='good'; 
                   c('2')='moderate'; 
                   c('3', '0')='poor'")) # including category 0, which is no (well-formed) glands
data$biopsies_diffgrade <- as.factor(car::recode(data$biopsies_diff, "c('good', 'moderate')='0'; c('poor')='1'"))

data$pCR_diffgrade <- paste(data$pCR, data$biopsies_diffgrade)
data$pCR_diffgrade <- as.factor(car::recode(data$pCR_diffgrade, "c('0 0')='non-pCR & good/mod'; c('0 1')='non-pCR & poor' ; c('1 0')='pCR & good/mod'; c('1 1')='pCR & poor'"))


# Re-code scores for resection specimen
# Define overall risk factor based on the max in each of the categories for every patient
data$res_rfall <- apply(data[, c('res_mucin', 'res_src', 'res_src_atypical')], MARGIN=1, max)
data$res_rfall <- as.factor(car::recode(data$res_rfall, "c('88')='0'")) # patient with TRG1 get lowest category
data$res_rfall_any <- as.factor(car::recode(data$res_rfall, "c('0')='0'; c('1', '2', '3')='1'"))
data$res_rfall_010 <- as.factor(car::recode(data$res_rfall, "c('0', '1')='0'; c('2', '3')='1'"))
data$res_rfall <- as.factor(car::recode(data$res_rfall, "c('0')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'"))

# Sum of signet-ring cells and poorly cohesive cells
# Re-code categories into percentiles
data$res_src_perc <- car::recode(data$res_src_num, "c('0')='0'; c('1') = '0.1'; c('2') = '0.4'; c('3')='0.5'")
data$res_src_atypical_perc <- car::recode(data$res_src_atypical_num, "c('0')='0'; c('1') = '0.1'; c('2') = '0.4'; c('3')='0.5'")
data$res_srcpcc_perc <- data$res_src_perc + data$res_src_atypical_perc
data$res_srcpcc_perc <- ifelse(data$trg == '1', '88', data$res_srcpcc_perc) # patients in TRG1 get assigned "other -> 0% category"
data$res_srcpcc_sum <- as.factor(car::recode(data$res_srcpcc_perc, "c('0', '88')='0'; c('0.1') = '1'; c('0.2', '0.4') = '2'; c('0.5', '0.6', '0.8', '0.9')='3'"))
# We re-examined the cases in which summation is not evident: 1-10% + 1-10%, 1-10% + 11-50%;, 11-50% + 11-50%
data$res_srcpcc_sum <- replace(data$res_srcpcc_sum, data$res_srcpcc_sum[data$Record.Id == "EMC130"], 2)
data$res_srcpcc_sum <- replace(data$res_srcpcc_sum, data$res_srcpcc_sum[data$Record.Id == "EMC137"], 3)

data$res_srcpcc_sum <- replace(data$res_srcpcc_sum, data$res_srcpcc_sum[data$Record.Id == "PATH011"], 1)
data$res_srcpcc_sum <- replace(data$res_srcpcc_sum, data$res_srcpcc_sum[data$Record.Id == "PATH014"], 2)

data$res_srcpcc_sum <- replace(data$res_srcpcc_sum, data$res_srcpcc_sum[data$Record.Id == "EMC155"], 1)
data$res_srcpcc_sum <- replace(data$res_srcpcc_sum, data$res_srcpcc_sum[data$Record.Id == "RAD006"], 2)

data$res_srcpcc_any <- as.factor(car::recode(data$res_srcpcc_sum, "c('0')='0'; c('1', '2', '3')='1'"))
data$res_srcpcc_010 <- as.factor(car::recode(data$res_srcpcc_sum, "c('0', '1')='0'; c('2', '3')='1'"))
data$res_srcpcc <- as.factor(car::recode(data$res_srcpcc_sum, "c('0')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'"))

data$res_mucin <- ifelse(data$trg == '1', '88', data$res_mucin) # patients with TRG 1 get assigned "other -> 0% category"
data$res_mucin_any <- as.factor(car::recode(data$res_mucin, "c('0', '88')='0'; c('1', '2', '3')='1'"))
data$res_mucin_010 <- as.factor(car::recode(data$res_mucin, "c('0', '88', '1')='0'; c('2', '3')='1'"))
data$res_mucin <- as.factor(car::recode(data$res_mucin, "c('0', '88')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'")) # There is one of "other" which was TRG1

data$res_src_any <- as.factor(car::recode(data$res_src, "c('0', '88')='0'; c('1', '2', '3')='1'"))
data$res_src_010 <- as.factor(car::recode(data$res_src, "c('0', '88', '1')='0'; c('2', '3')='1'"))
data$res_src <- as.factor(car::recode(data$res_src, "c('0', '88')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'")) # There is one of "other"/88 which was TRG1

# signet-ring cells plus poorly cohesive cells before and after nCRT
data$srcpcc_categories <- as.factor(paste(data$biopsies_srcpcc_any, data$res_srcpcc_any))
data$srcpcc_categories <- as.factor(car::recode(data$srcpcc_categories, "c('0 0')='absent/absent'; 
                   c('0 1')='absent/present'; 
                   c('0 NA')='absent, no res';
                   c('1 0')='present/absent';
                   c('1 1')='present/present';
                   c('1 NA')='present, no res'")) 

data$src_categories <- as.factor(paste(data$biopsies_src_any, data$res_src_any))
data$src_categories <- as.factor(car::recode(data$src_categories, "c('0 0')='absent/absent'; 
                   c('0 1')='absent/present'; 
                   c('0 NA')='absent, no res';
                   c('1 0')='present/absent';
                   c('1 1')='present/present';
                   c('1 NA')='present, no res'")) 

data$pcc_categories <- as.factor(paste(data$biopsies_src_atypical_any, data$res_src_atypical_any))
data$pcc_categories <- as.factor(car::recode(data$pcc_categories, "c('0 0')='absent/absent'; 
                   c('0 1')='absent/present'; 
                   c('0 NA')='absent, no res';
                   c('1 0')='present/absent';
                   c('1 1')='present/present';
                   c('1 NA')='present, no res'")) 

data$res_src_atypical_any <- as.factor(car::recode(data$res_src_atypical, "c('0', '88')='0'; c('1', '2', '3')='1'"))
data$res_src_atypical_010 <- as.factor(car::recode(data$res_src_atypical, "c('0', '1')='0'; c('2', '3')='1'"))
data$res_src_atypical <- as.factor(car::recode(data$res_src_atypical, "c('0', '88')='<1%'; 
                   c('1')='1-10%'; 
                   c('2')='11-50%'; 
                   c('3')='51-100%'")) # There is one of "other"/88 which was TRG1

# Differentiation grade
data$res_diff <- ifelse(data$trg == '1', 'NA, TRG1', data$res_diff) 
data$res_diff <- as.factor(car::recode(data$res_diff, "c('1')='good'; 
                   c('2')='moderate'; 
                   c('3', '0')='poor'")) # including category 0, which is no (well-formed) glands
data$res_diffgrade <- as.factor(car::recode(data$res_diff, "c('good', 'moderate')='0'; c('poor')='1'"))

# Re-code outcomes
data$death_cause <- car::recode(data$death_cause, "c('1')='esophageal cancer'; c('2')='other cause'; c('99')='unknown'")

data$status <- as.integer(car::recode(data$status, "c('1')='1'; c('2')='0'")) #1 is dead, 0 (was 2) is alive

data$dfstatus <- as.numeric(data$recurrence) + as.numeric(data$status)
data$dfstatus <- as.integer(car::recode(data$dfstatus, "c('0')='0'; c('1', '2')='1'")) #1 is recurrence-dead, 0 is no recurrence

# Define overall and disease-free survival
data$tvar <- as.numeric(coalesce(data$ostime_death, data$ostime_last_fu)) # join date of death and last follow-up date into one variable
data$dfsvar <- as.numeric(coalesce(data$dfstime_recurrence, data$tvar)) # join date of recurrence with death/last FU into one variable

data$tvar <- data$tvar / (365/12) # convert overall survival days to months
data$dfsvar <- data$dfsvar / (365/12) # convert disease-free survival days to months

data$potentialFU <- lubridate::interval(data$date_diagnosis, as.Date("15-07-2022", format = '%d-%m-%Y')) %/% months(1)

write.csv(data, "output/dataClean.csv", row.names = TRUE, na = "")
