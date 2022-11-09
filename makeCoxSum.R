library(broom)

makeCoxSum <- function(fit_null, fit_alt) {
  
  # Export Cox summaries of null model (without differentiation grade) and alternative model (with differentiation grade)
  coxmodel <- tidy(fit_null, conf.int = TRUE, exponentiate = TRUE)
  coxmodel[2:7] <- round(coxmodel[2:7], digits =2)
  coxmodel$estCI <- paste0(coxmodel$estimate, " (", coxmodel$conf.low, " - ", coxmodel$conf.high, ")")
  coxmodel
  print("Cox null:")
  print(coxmodel[1, c("term", "p.value", "estCI")])
  summary(fit_alt)
  coxmodel <- tidy(fit_alt, conf.int = TRUE, exponentiate = TRUE)
  coxmodel[2:7] <- round(coxmodel[2:7], digits =2)
  coxmodel$estCI <- paste0(coxmodel$estimate, " (", coxmodel$conf.low, " - ", coxmodel$conf.high, ")")
  coxmodel
  print("Cox alt with differentiation grade")
  print(coxmodel[1, c("term", "p.value", "estCI")])
  return(coxmodel)
  
}

makeLRSum <- function(fit_null, fit_alt) {
  
  # Export Cox summaries of null model (without differentiation grade) and alternative model (with differentiation grade)
  coxmodel <- tidy(fit_null, conf.int = TRUE, exponentiate = TRUE)
  coxmodel[2:7] <- round(coxmodel[2:7], digits =2)
  coxmodel$estCI <- paste0(coxmodel$estimate, " (", coxmodel$conf.low, " - ", coxmodel$conf.high, ")")
  coxmodel
  print("LR null:")
  print(coxmodel[2, c("term", "p.value", "estCI")])
  summary(fit_alt)
  coxmodel <- tidy(fit_alt, conf.int = TRUE, exponentiate = TRUE)
  coxmodel[2:7] <- round(coxmodel[2:7], digits =2)
  coxmodel$estCI <- paste0(coxmodel$estimate, " (", coxmodel$conf.low, " - ", coxmodel$conf.high, ")")
  coxmodel
  print("LR alt with differentiation grade")
  print(coxmodel[2, c("term", "p.value", "estCI")])
  return(coxmodel)
  
}

