suppressPackageStartupMessages({library(mice);library(sandwich);library(jsonlite);library(splines)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5) stop("Usage: Rscript fit_historical_falls.R PREPARED HISTORY_INPUTS ORIGINAL_CACHE OUTPUT HISTORY_CACHE")
script <- sub("^--file=", "", commandArgs()[grep("^--file=",commandArgs())][1])
code <- normalizePath(file.path(dirname(script),".."),winslash="/")
source(file.path(code,"workflow/path_guard.R"))
for(i in 4:5) {assert_external(args[i]);dir.create(args[i],recursive=TRUE,showWarnings=FALSE)}
# Import only these frozen helpers; do not execute the original analysis driver.
needed <- c("primary_full","pool_scalar","load_cohort","fit_one","pool_model")
loaded <- character()
for(e in parse(file.path(code,"charls_elsa_models/fit_models.R"))) {
 if(is.call(e) && identical(e[[1]],as.name("<-")) && is.symbol(e[[2]]) && as.character(e[[2]]) %in% needed) {
   eval(e);loaded <- c(loaded,as.character(e[[2]]))
 }
}
stopifnot(setequal(loaded,needed))
rows <- list();diagnostics <- list();capacity <- list();verification <- list();changes <- list()
for(co in c("CHARLS","ELSA")) {
 cut <- if(co=="CHARLS") 10 else 4
 obj <- load_cohort(co,cut,file.path(args[1],paste0(co,"_analysis.csv")),file.path(args[3],paste0(co,"_fall_",cut,"_mi.rds")))
 hr <- read.csv(file.path(args[2],paste0(co,"_history_falls.csv")),stringsAsFactors=FALSE,check.names=FALSE)
 stopifnot(!anyDuplicated(hr$id),nrow(hr)==nrow(obj$raw))
 stopifnot(isTRUE(all.equal(hr[names(obj$raw)],obj$raw,check.attributes=TRUE)))
 hf <- hr$history_fall[match(obj$d$id,hr$id)]
 # Verify the locked primary estimates and interval restriction.
 for(label in c("frozen_primary_check","both_intervals_18_30")) {
  z <- obj
  if(label=="both_intervals_18_30") z$sets <- lapply(obj$sets,function(x) x[x$exposure_gap_months>=18 & x$exposure_gap_months<=30 & x$followup_gap_months>=18 & x$followup_gap_months<=30,,drop=FALSE])
  a <- pool_model(z,primary_full,label);rows[[length(rows)+1]] <- a$row;diagnostics[[length(diagnostics)+1]] <- a$diag
  if(label=="frozen_primary_check") {
   expected <- read.csv(file.path(code,"charls_elsa_models/aggregate_results/sequential_models.csv"));expected <- expected[expected$cohort==co & expected$model=="full_adjustment",]
  } else {
   expected <- read.csv(file.path(code,"charls_elsa_models/aggregate_results/robustness_sensitivity.csv"));expected <- expected[expected$cohort==co & expected$sensitivity==label,]
  }
  cols <- intersect(c("n","events","rr","ci_low","ci_high","risk0","risk1","risk_difference","rd_low","rd_high"),names(expected))
  pass <- all(vapply(cols,function(v) isTRUE(all.equal(a$row[[v]],expected[[v]],tolerance=1e-8)),logical(1)))
  verification[[length(verification)+1]] <- data.frame(cohort=co,model=label,pass=pass,columns=paste(cols,collapse=";"))
  stopifnot(pass)
 }
 keep <- !is.na(hf)
 d <- obj$d[keep,,drop=FALSE];d$history_fall <- factor(hf[keep],levels=0:1)
 im <- d[,c(all.vars(primary_full),"history_fall"),drop=FALSE]
 cache <- file.path(args[5],paste0(co,"_history_falls_mi.rds"))
 if(file.exists(cache)) imp <- readRDS(cache) else {
  meth <- mice::make.method(im);pred <- mice::make.predictorMatrix(im)
  for(v in names(im)) {
   if(!anyNA(im[[v]])) meth[v] <- "" else if(is.factor(im[[v]])) meth[v] <- if(nlevels(im[[v]])==2) "logreg" else "polyreg" else meth[v] <- "pmm"
  }
  message(co,": starting paired imputations")
  imp <- mice::mice(im,m=30,maxit=20,method=meth,predictorMatrix=pred,seed=if(co=="CHARLS")202611151 else 202611152,printFlag=FALSE)
  saveRDS(imp,cache)
 }
 stopifnot(imp$m==30,imp$iteration==20,isTRUE(all.equal(imp$data,im,check.attributes=TRUE)))
 paired <- obj;paired$sets <- lapply(seq_len(imp$m),function(i) {
  z <- mice::complete(imp,i);z$id <- as.character(d$id);z$cluster <- as.character(d$cluster);z$cohort <- co;z
 })
 ev <- imp$loggedEvents
 if(is.null(ev)) ev <- data.frame(it=integer(),im=integer(),dep=character(),meth=character(),out=character())
 write.csv(ev,file.path(args[4],paste0(co,"_imputation_logged_events.csv")),row.names=FALSE)
 cv <- mice::convergence(imp)
 write.csv(cv,file.path(args[4],paste0(co,"_imputation_convergence.csv")),row.names=FALSE)
 before <- pool_model(paired,primary_full,"same_sample_original")
 after <- pool_model(paired,update(primary_full,.~.+history_fall),"same_sample_plus_history_fall")
 stopifnot(before$row$n==after$row$n,before$row$events==after$row$events)
 for(a in list(before,after)) {rows[[length(rows)+1]]<-a$row;diagnostics[[length(diagnostics)+1]]<-a$diag}
 changes[[length(changes)+1]] <- data.frame(cohort=co,n=before$row$n,events=before$row$events,rr_before=before$row$rr,rr_after=after$row$rr,log_rr_change=after$row$log_rr-before$row$log_rr,rr_ratio=after$row$rr/before$row$rr)
 nr <- d[d$history_fall=="0",,drop=FALSE]
 ns <- table(factor(nr$hist,levels=0:1));es <- tapply(nr$y,factor(nr$hist,levels=0:1),sum)
 feasible <- all(ns>=100) && all(es>=20) && sum(nr$y)>=200 && length(unique(nr$cluster))>=40
 capacity[[length(capacity)+1]] <- data.frame(cohort=co,missing_history=sum(!keep),available_n=sum(keep),restricted_n=nrow(nr),restricted_events=sum(nr$y),hist0_n=ns[1],hist1_n=ns[2],hist0_events=es[1],hist1_events=es[2],restriction_feasible=feasible)
 if(feasible) {
  two <- paired;two$sets <- lapply(paired$sets,function(z)z[z$history_fall=="0",,drop=FALSE])
  a <- pool_model(two,primary_full,"no_falls_at_history_and_index")
  rows[[length(rows)+1]]<-a$row;diagnostics[[length(diagnostics)+1]]<-a$diag
 }
 message(co,": completed")
}
rr <- do.call(rbind,rows);dd <- do.call(rbind,diagnostics)
stopifnot(all(dd$converged),all(is.finite(rr$rr)),all(is.finite(rr$se_log_rr)))
write.csv(rr,file.path(args[4],"historical_falls_models.csv"),row.names=FALSE)
write.csv(dd,file.path(args[4],"model_diagnostics.csv"),row.names=FALSE)
write.csv(do.call(rbind,changes),file.path(args[4],"paired_point_estimate_changes.csv"),row.names=FALSE)
write.csv(do.call(rbind,capacity),file.path(args[4],"capacity_and_missingness.csv"),row.names=FALSE)
write.csv(do.call(rbind,verification),file.path(args[4],"locked_result_verification.csv"),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[4],"sessionInfo.txt"))
print(rr[,c("cohort","model","n","events","rr","ci_low","ci_high","risk0","risk1","risk_difference","rd_low","rd_high")],row.names=FALSE)
