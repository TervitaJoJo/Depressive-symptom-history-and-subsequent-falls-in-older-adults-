suppressPackageStartupMessages({library(mice);library(sandwich);library(splines);library(mitml)})
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

