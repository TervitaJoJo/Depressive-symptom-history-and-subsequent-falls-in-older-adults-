suppressPackageStartupMessages({
  library(mice)
  library(sandwich)
  library(jsonlite)
  library(splines)
})

args <- commandArgs(trailingOnly=TRUE)
script_path <- sub("^--file=", "", commandArgs()[grep("^--file=",commandArgs())][1])
source(file.path(dirname(script_path),"../workflow/path_guard.R"))
if(length(args)>=2) assert_external(args[2])
if(length(args)>=3) assert_external(args[3])

if(length(args)<2) stop("Usage: Rscript fit_models.R INPUT_DIRECTORY OUTPUT_DIRECTORY [CACHE_DIRECTORY]")
base <- normalizePath(args[1], winslash="/", mustWork=TRUE)
outdir <- normalizePath(args[2], winslash="/", mustWork=FALSE)
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(outdir,"logs"), recursive=TRUE, showWarnings=FALSE)

cache_dir <- if(length(args)>=3) args[3] else file.path(outdir,"local_cache")
dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
charls_csv <- file.path(base,"CHARLS_analysis.csv")
elsa_csv <- file.path(base,"ELSA_analysis.csv")
charls_mi <- file.path(cache_dir,"CHARLS_fall_10_mi.rds")
elsa_mi <- file.path(cache_dir,"ELSA_fall_4_mi.rds")

primary_full <- y ~ hist + ns(age,df=3) + female + education + window + exposure_gap_months + followup_gap_months + current_score + partnered + smoking + alcohol + stroke + hearte + diabe + arthre + stair_difficulty + adl_difficulty
base_temporal <- y ~ ns(age,df=3) + female + education + window + exposure_gap_months + followup_gap_months
current_only <- update(base_temporal, . ~ . + current_score)
current_history <- update(current_only, . ~ . + hist)
history_function <- update(current_history, . ~ . + stair_difficulty + adl_difficulty)

pool_scalar <- function(q,u,dfcom) {
  q <- as.numeric(q); u <- as.numeric(u); m <- length(q)
  qb <- mean(q); ub <- mean(u); b <- if(m>1) var(q) else 0
  tt <- ub + (1 + 1/m)*b
  se <- sqrt(tt)
  lambda <- if(tt > 0) (1 + 1/m)*b/tt else 0
  if(m==1 || lambda < 1e-12) {
    df <- dfcom
  } else {
    old <- (m-1)/lambda^2
    obs <- (dfcom+1)/(dfcom+3)*dfcom*(1-lambda)
    df <- 1/(1/old + 1/obs)
  }
  crit <- qt(.975,df)
  c(est=qb,se=se,lo=qb-crit*se,hi=qb+crit*se,
    p=2*pt(-abs(qb/se),df),df=df,lambda=lambda)
}

load_cohort <- function(co, cut, csv, cache, endpoint="fall") {
  raw <- read.csv(csv, stringsAsFactors=FALSE, check.names=FALSE)
  d <- raw[raw$cesd_current < cut & !is.na(raw[[endpoint]]),,drop=FALSE]
  d$y <- d[[endpoint]]
  d$hist <- factor(as.integer(d$cesd_prior >= cut), levels=0:1)
  d$current_score <- factor(d$cesd_current, levels=0:(cut-1))
  for(v in c("education","smoking","window")) d[[v]] <- factor(d[[v]])
  for(v in c("partnered","alcohol","stroke","hearte","diabe","arthre","stair_difficulty","adl_difficulty")) d[[v]] <- factor(d[[v]],levels=0:1)
  d$female <- factor(d$female,levels=0:1)
  vars <- all.vars(primary_full)
  im <- d[,vars,drop=FALSE]
  if(file.exists(cache)) {
    imp <- readRDS(cache)
  } else {
    meth <- mice::make.method(im); pred <- mice::make.predictorMatrix(im)
    for(v in names(im)) {
      if(!anyNA(im[[v]])) meth[v] <- "" else if(is.factor(im[[v]])) meth[v] <- if(nlevels(im[[v]])==2) "logreg" else "polyreg" else meth[v] <- "pmm"
    }
    seed_value <- if(co=="CHARLS") 20261331 else 20261230+cut+match(endpoint,c("fall","recurrent_fall","injurious_fall"))
    imp <- mice::mice(im,m=30,maxit=20,method=meth,predictorMatrix=pred,seed=seed_value,printFlag=FALSE)
    saveRDS(imp,cache)
  }
  stopifnot(imp$m==30,imp$iteration==20)
  if(!isTRUE(all.equal(imp$data, im, check.attributes=TRUE))) {
    stop(paste0(co," cached imputation data do not match the frozen primary analysis input"))
  }
  sets <- lapply(seq_len(imp$m), function(i) {
    z <- mice::complete(imp,i)
    z$cluster <- as.character(d$cluster)
    z$id <- as.character(d$id)
    z$cohort <- co
    if("survey_weight" %in% names(d)) z$survey_weight <- d$survey_weight
    z$original_adl <- d$adl_difficulty
    z$current_score_num <- as.numeric(as.character(z$current_score))
    z
  })
  list(cohort=co,cut=cut,raw=raw,d=d,imp=imp,sets=sets)
}

charls <- load_cohort("CHARLS",10,charls_csv,charls_mi)
elsa <- load_cohort("ELSA",4,elsa_csv,elsa_mi)
cohorts <- list(CHARLS=charls,ELSA=elsa)

fit_one <- function(dat, form, weight_col=NULL) {
  if("window" %in% all.vars(form) && length(unique(dat$window))==1) form <- update(form,.~.-window)
  if(is.null(weight_col)) {
    ww <- rep(1,nrow(dat))
  } else {
    ww <- dat[[weight_col]]
    if(any(!is.finite(ww)) || any(ww < 0)) stop("Invalid analysis weights")
    if(sum(ww) <= 0) stop("No positive analysis weights")
    ww <- ww / mean(ww)
  }
  dat$.analysis_weight <- ww
  fit <- glm(form, data=dat, family=poisson(link="log"), weights=.analysis_weight,
             control=glm.control(epsilon=1e-12,maxit=100))
  if(anyNA(coef(fit))) stop("Aliased model coefficients")
  vv <- sandwich::vcovCL(fit, cluster=dat$cluster, type="HC1")
  out <- list(fit=fit,vv=vv,dat=dat,weights=ww,
              n=nrow(dat),events=sum(dat$y),clusters=length(unique(dat$cluster)),
              p=length(coef(fit)),dfcom=max(1,length(unique(dat$cluster))-length(coef(fit))),
              converged=fit$converged,max_fitted=max(fitted(fit)),n_fitted_gt1=sum(fitted(fit)>1))
  if("hist1" %in% names(coef(fit))) {
    b <- coef(fit); k <- which(names(b)=="hist1")
    d0 <- dat; d1 <- dat
    d0$hist <- factor("0",levels=levels(dat$hist))
    d1$hist <- factor("1",levels=levels(dat$hist))
    x0 <- model.matrix(delete.response(terms(fit)),d0,contrasts.arg=fit$contrasts,xlev=fit$xlevels)
    x1 <- model.matrix(delete.response(terms(fit)),d1,contrasts.arg=fit$contrasts,xlev=fit$xlevels)
    p0 <- as.vector(exp(x0 %*% b)); p1 <- as.vector(exp(x1 %*% b))
    ss <- ww/sum(ww)
    grad <- colSums(x1*(ss*p1))-colSums(x0*(ss*p0))
    out$q <- unname(b[k]); out$u <- unname(vv[k,k])
    out$risk0 <- sum(ss*p0); out$risk1 <- sum(ss*p1)
    out$rd <- sum(ss*(p1-p0)); out$rdu <- as.numeric(t(grad)%*%vv%*%grad)
  } else {
    out$q <- NA_real_; out$u <- NA_real_; out$risk0 <- NA_real_; out$risk1 <- NA_real_; out$rd <- NA_real_; out$rdu <- NA_real_
  }
  out
}

pool_model <- function(obj, form, label, weight_col=NULL) {
  vals <- lapply(obj$sets, function(z) fit_one(z,form,weight_col))
  has_hist <- all(vapply(vals,function(v)!is.na(v$q),logical(1)))
  if(has_hist) {
    rr <- pool_scalar(vapply(vals,`[[`,numeric(1),"q"),vapply(vals,`[[`,numeric(1),"u"),min(vapply(vals,`[[`,numeric(1),"dfcom")))
    rd <- pool_scalar(vapply(vals,`[[`,numeric(1),"rd"),vapply(vals,`[[`,numeric(1),"rdu"),min(vapply(vals,`[[`,numeric(1),"dfcom")))
    row <- data.frame(cohort=obj$cohort,model=label,m=length(vals),n=vals[[1]]$n,events=vals[[1]]$events,clusters=vals[[1]]$clusters,
      log_rr=unname(rr["est"]),se_log_rr=unname(rr["se"]),rr=exp(unname(rr["est"])),ci_low=exp(unname(rr["lo"])),ci_high=exp(unname(rr["hi"])),p=unname(rr["p"]),df=unname(rr["df"]),lambda=unname(rr["lambda"]),
      risk0=mean(vapply(vals,`[[`,numeric(1),"risk0")),risk1=mean(vapply(vals,`[[`,numeric(1),"risk1")),risk_difference=unname(rd["est"]),rd_low=unname(rd["lo"]),rd_high=unname(rd["hi"]),weighting=ifelse(is.null(weight_col),"unweighted",weight_col),stringsAsFactors=FALSE)
  } else {
    row <- data.frame(cohort=obj$cohort,model=label,m=length(vals),n=vals[[1]]$n,events=vals[[1]]$events,clusters=vals[[1]]$clusters,
      log_rr=NA,se_log_rr=NA,rr=NA,ci_low=NA,ci_high=NA,p=NA,df=NA,lambda=NA,risk0=NA,risk1=NA,risk_difference=NA,rd_low=NA,rd_high=NA,weighting=ifelse(is.null(weight_col),"unweighted",weight_col),stringsAsFactors=FALSE)
  }
  diag <- data.frame(cohort=obj$cohort,model=label,imputation=seq_along(vals),n=vapply(vals,`[[`,numeric(1),"n"),events=vapply(vals,`[[`,numeric(1),"events"),clusters=vapply(vals,`[[`,numeric(1),"clusters"),converged=vapply(vals,`[[`,logical(1),"converged"),max_fitted=vapply(vals,`[[`,numeric(1),"max_fitted"),n_fitted_gt1=vapply(vals,`[[`,numeric(1),"n_fitted_gt1"),weighting=ifelse(is.null(weight_col),"unweighted",weight_col))
  list(row=row,diag=diag,vals=vals)
}

seq_forms <- list(
  current_symptoms_only=current_only,
  current_symptoms_plus_history=current_history,
  plus_current_function=history_function,
  full_adjustment=primary_full
)
seq <- list(); diags <- list(); k <- 0L
for(co in names(cohorts)) for(nm in names(seq_forms)) {
  k <- k+1L; ans <- pool_model(cohorts[[co]],seq_forms[[nm]],nm)
  seq[[k]] <- ans$row; diags[[k]] <- ans$diag
}
seq_df <- do.call(rbind,seq); diag_df <- do.call(rbind,diags)
write.csv(seq_df,file.path(outdir,"sequential_models.csv"),row.names=FALSE)
write.csv(diag_df,file.path(outdir,"sequential_model_diagnostics.csv"),row.names=FALSE)

# Joint Wald tests for the current-symptom block and for the added history term.
joint_test <- function(obj, form, block, label) {
  vals <- lapply(obj$sets,function(z) fit_one(z,form))
  cn <- names(coef(vals[[1]]$fit))
  terms <- if(block=="current_score") grep("^current_score",cn,value=TRUE) else grep("^hist",cn,value=TRUE)
  if(!length(terms)) stop(paste("No terms found for",block))
  Q <- do.call(cbind,lapply(vals,function(v)coef(v$fit)[terms]))
  U <- array(NA_real_,dim=c(length(terms),length(terms),length(vals)))
  for(i in seq_along(vals)) U[,,i] <- vals[[i]]$vv[terms,terms,drop=FALSE]
  tt <- mitml::testConstraints(qhat=Q,uhat=U,constraints=terms,method="D1",df.com=min(vapply(vals,`[[`,numeric(1),"dfcom")))
  tr <- tt$test[1,]
  data.frame(cohort=obj$cohort,model=label,block=block,terms=paste(terms,collapse=";"),m=length(vals),F=unname(tr[1]),df1=unname(tr[2]),df2=unname(tr[3]),p=unname(tr[4]),riv=unname(tr[5]),stringsAsFactors=FALSE)
}
joint_rows <- list(); k <- 0L
for(co in names(cohorts)) {
  k <- k+1L; joint_rows[[k]] <- joint_test(cohorts[[co]],current_only,"current_score","current_symptoms_only")
  k <- k+1L; joint_rows[[k]] <- joint_test(cohorts[[co]],current_history,"current_score","current_symptoms_plus_history")
  k <- k+1L; joint_rows[[k]] <- joint_test(cohorts[[co]],current_history,"history","current_symptoms_plus_history")
}
write.csv(do.call(rbind,joint_rows),file.path(outdir,"incremental_joint_wald.csv"),row.names=FALSE)

# Formal cross-cohort heterogeneity based on the same frozen full model.
primary_rows <- seq_df[seq_df$model=="full_adjustment",]
stopifnot(nrow(primary_rows)==2)
theta <- primary_rows$log_rr; se <- primary_rows$se_log_rr; names(theta)<-primary_rows$cohort; names(se)<-primary_rows$cohort
w <- 1/se^2; theta_fe <- sum(w*theta)/sum(w); se_fe <- sqrt(1/sum(w));
qhet <- sum(w*(theta-theta_fe)^2); phet <- pchisq(qhet,df=1,lower.tail=FALSE); i2 <- ifelse(qhet>0,max(0,(qhet-1)/qhet)*100,0)
diff <- theta["CHARLS"]-theta["ELSA"]; sediff <- sqrt(se["CHARLS"]^2+se["ELSA"]^2); z <- diff/sediff; pdiff <- 2*pnorm(-abs(z))
het <- rbind(
  data.frame(row_type="cohort_specific",cohort="CHARLS",log_rr=theta["CHARLS"],se_log_rr=se["CHARLS"],rr=exp(theta["CHARLS"]),ci_low=exp(theta["CHARLS"]-1.96*se["CHARLS"]),ci_high=exp(theta["CHARLS"]+1.96*se["CHARLS"]),p=primary_rows$p[primary_rows$cohort=="CHARLS"],Q=NA,Q_df=NA,Q_p=NA,I2=NA,heterogeneity_contrast=NA,contrast_se=NA,contrast_p=NA),
  data.frame(row_type="cohort_specific",cohort="ELSA",log_rr=theta["ELSA"],se_log_rr=se["ELSA"],rr=exp(theta["ELSA"]),ci_low=exp(theta["ELSA"]-1.96*se["ELSA"]),ci_high=exp(theta["ELSA"]+1.96*se["ELSA"]),p=primary_rows$p[primary_rows$cohort=="ELSA"],Q=NA,Q_df=NA,Q_p=NA,I2=NA,heterogeneity_contrast=NA,contrast_se=NA,contrast_p=NA),
  data.frame(row_type="fixed_effect",cohort="CHARLS+ELSA",log_rr=theta_fe,se_log_rr=se_fe,rr=exp(theta_fe),ci_low=exp(theta_fe-1.96*se_fe),ci_high=exp(theta_fe+1.96*se_fe),p=2*pnorm(-abs(theta_fe/se_fe)),Q=qhet,Q_df=1,Q_p=phet,I2=i2,heterogeneity_contrast=diff,contrast_se=sediff,contrast_p=pdiff)
)
write.csv(het,file.path(outdir,"charls_elsa_heterogeneity.csv"),row.names=FALSE)

# Follow-up observation rates and inverse-probability-of-observation sensitivity analysis.
prepare_ipw <- function(obj) {
  r <- obj$raw
  e <- r[r$cesd_current < obj$cut,,drop=FALSE]
  e$observed <- as.integer(!is.na(e$fall))
  e$hist <- factor(as.integer(e$cesd_prior >= obj$cut),levels=0:1)
  e$current_score <- factor(e$cesd_current,levels=0:(obj$cut-1))
  for(v in c("education","smoking","window")) e[[v]] <- factor(e[[v]])
  for(v in c("partnered","alcohol","stroke","hearte","diabe","arthre","stair_difficulty","adl_difficulty")) e[[v]] <- factor(e[[v]],levels=0:1)
  e$female <- factor(e$female,levels=0:1)
  group_label <- ifelse(e$hist=="1","prior_high_current_low","prior_low_current_low")
  rate_df <- do.call(rbind,lapply(unique(group_label),function(g){
    x <- e$observed[group_label==g]
    data.frame(cohort=obj$cohort,group=g,n=length(x),observed=sum(x),missing=sum(!x),observed_fraction=mean(x),stringsAsFactors=FALSE)
  }))
  # Complete the observation model by explicit missing categories for categorical covariates and medians for numeric predictors.
  pred_vars <- c("hist","age","female","education","window","exposure_gap_months","current_score","partnered","smoking","alcohol","stroke","hearte","diabe","arthre","stair_difficulty","adl_difficulty")
  for(v in pred_vars) {
    if(is.numeric(e[[v]])) {
      med <- median(e[[v]],na.rm=TRUE); e[[v]][is.na(e[[v]])] <- med
    } else {
      x <- as.character(e[[v]]); x[is.na(x)] <- "Missing"; e[[v]] <- factor(x)
    }
  }
  obs_form <- observed ~ hist + ns(age,df=3) + female + education + window + exposure_gap_months + current_score + partnered + smoking + alcohol + stroke + hearte + diabe + arthre + stair_difficulty + adl_difficulty
  ofit <- glm(obs_form,data=e,family=binomial())
  phat <- predict(ofit,type="response")
  phat <- pmax(pmin(phat,1-1e-6),1e-6)
  pbar <- mean(e$observed)
  e$ipw <- pbar/phat
  e$group_label <- group_label
  wsummary <- data.frame(cohort=obj$cohort,eligible=nrow(e),observed=sum(e$observed),pbar=pbar,min_p=min(phat),max_p=max(phat),min_ipw=min(e$ipw[e$observed==1]),p50_ipw=median(e$ipw[e$observed==1]),p99_ipw=quantile(e$ipw[e$observed==1],.99),max_ipw=max(e$ipw[e$observed==1]),ess=(sum(e$ipw[e$observed==1])^2/sum(e$ipw[e$observed==1]^2)),observation_model_n=ofit$df.null+1)
  idmap <- setNames(e$ipw,as.character(e$id))
  list(e=e,rate=rate_df,wsum=wsummary,idmap=idmap,fit=ofit)
}

ipw_objs <- lapply(cohorts,prepare_ipw)
followup_df <- do.call(rbind,lapply(ipw_objs,`[[`,"rate"))
write.csv(followup_df,file.path(outdir,"followup_observation.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(ipw_objs,`[[`,"wsum")),file.path(outdir,"followup_ipw_diagnostics.csv"),row.names=FALSE)

ipw_results <- list(); k <- 0L
for(co in names(cohorts)) {
  obj <- cohorts[[co]]; idmap <- ipw_objs[[co]]$idmap
  sets_ipw <- lapply(obj$sets,function(z){z$ipw <- unname(idmap[as.character(z$id)]);z})
  obj2 <- obj; obj2$sets <- sets_ipw
  ans <- pool_model(obj2,primary_full,"followup_IPW",weight_col="ipw")
  k <- k+1L; ipw_results[[k]] <- ans$row
}
ipw_df <- do.call(rbind,ipw_results)
write.csv(ipw_df,file.path(outdir,"followup_ipw_results.csv"),row.names=FALSE)

# ELSA survey-weight sensitivity; CHARLS has no directly corresponding weight in the frozen harmonized extract.
wt_rows <- list()
ans_elsa_w <- pool_model(elsa,primary_full,"survey_weighted",weight_col="survey_weight")
wt_rows[[1]] <- ans_elsa_w$row
wt_rows[[2]] <- data.frame(cohort="CHARLS",model="survey_weighted",m=NA,n=NA,events=NA,clusters=NA,log_rr=NA,se_log_rr=NA,rr=NA,ci_low=NA,ci_high=NA,p=NA,df=NA,lambda=NA,risk0=NA,risk1=NA,risk_difference=NA,rd_low=NA,rd_high=NA,weighting="unavailable_in_frozen_extract",stringsAsFactors=FALSE)
write.csv(do.call(rbind,wt_rows),file.path(outdir,"weighted_sensitivity.csv"),row.names=FALSE)


# Restricted analyses retain the same frozen primary imputation objects.
robust <- list(); kk <- 0L
for(co in names(cohorts)) {
  ob <- cohorts[[co]]
  cc <- ob$d[complete.cases(ob$d[,all.vars(primary_full)]),,drop=FALSE]
  cc$original_adl <- cc$adl_difficulty
  oc <- ob; oc$sets <- list(cc)
  kk <- kk+1L; robust[[kk]] <- pool_model(oc,primary_full,"complete_case")$row
  oi <- ob; oi$sets <- lapply(ob$sets,function(z) z[z$exposure_gap_months>=18 & z$exposure_gap_months<=30 & z$followup_gap_months>=18 & z$followup_gap_months<=30,,drop=FALSE])
  kk <- kk+1L; robust[[kk]] <- pool_model(oi,primary_full,"both_intervals_18_30")$row
  oa <- ob; oa$sets <- lapply(ob$sets,function(z) z[!is.na(z$original_adl) & z$original_adl=="0",,drop=FALSE])
  kk <- kk+1L; robust[[kk]] <- pool_model(oa,update(primary_full,.~.-adl_difficulty),"no_current_adl")$row
}
write.csv(do.call(rbind,robust),file.path(outdir,"restricted_analyses.csv"),row.names=FALSE)
# Optional extra models can require fresh, computationally expensive imputations.
if(tolower(Sys.getenv("RUN_EXTRA_MODELS","false"))=="true") {
  ex <- list()
  o3 <- load_cohort("ELSA",3,elsa_csv,file.path(cache_dir,"ELSA_fall_3_mi.rds"))
  ex[[1]] <- pool_model(o3,primary_full,"alternative_threshold_3")$row
  for(ep in c("recurrent_fall","injurious_fall")) {
    oo <- load_cohort("ELSA",4,elsa_csv,file.path(cache_dir,paste0("ELSA_",ep,"_4_mi.rds")),endpoint=ep)
    ex[[length(ex)+1L]] <- pool_model(oo,primary_full,ep)$row
  }
  write.csv(do.call(rbind,ex),file.path(outdir,"elsa_extra_models.csv"),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(outdir,"R_sessionInfo.txt"))
message("Completed main, restricted, observation-weighted and ELSA response-weighted models.")
