suppressPackageStartupMessages({library(mice);library(sandwich);library(jsonlite);library(splines);library(mitml)})
args <- commandArgs(trailingOnly=TRUE)
script_path <- sub("^--file=", "", commandArgs()[grep("^--file=",commandArgs())][1])
source(file.path(dirname(script_path),"../workflow/path_guard.R"))
if(length(args)>=2) {assert_external(args[1]); assert_external(args[2])}

if(length(args)<2)stop('Usage: Rscript fit_four_visit.R INPUT_DIRECTORY OUTPUT_DIRECTORY [preflight]')
base <- normalizePath(args[1], winslash='/', mustWork=TRUE)
resdir <- normalizePath(args[2],winslash='/',mustWork=FALSE); dir.create(resdir,recursive=TRUE,showWarnings=FALSE)
private <- base
mode <- if(length(args)>2) args[3] else 'fit'
d <- read.csv(file.path(private,'CHARLS_four_visit_observed.csv'),colClasses=c(id='character',household='character',cluster='character'))
stopifnot(nrow(d)==1443,sum(d$y)==230,!anyDuplicated(d$id))
d$history <- factor(d$history,levels=0:2)
d$current_score <- factor(d$current_score,levels=0:9)
d$education <- factor(d$education,levels=1:3)
d$smoking <- factor(d$smoking,levels=0:2)
for(v in c('female','partnered','alcohol','stroke','hearte','diabe','arthre','stair_difficulty','adl_difficulty'))d[[v]] <- factor(d[[v]],levels=0:1)
age_basis <- ns(d$age,df=3)
age_knots <- attr(age_basis,'knots'); age_boundary <- attr(age_basis,'Boundary.knots')
full <- y ~ history + ns(age,knots=age_knots,Boundary.knots=age_boundary) + female + education + current_score + partnered + smoking + alcohol + stroke + hearte + diabe + arthre + stair_difficulty + adl_difficulty + history_gap_months + landmark_gap_months + followup_gap_months
vars <- all.vars(full); vars <- setdiff(vars,c('age_knots','age_boundary'))
im <- d[,vars,drop=FALSE]
meth <- mice::make.method(im); pred <- mice::make.predictorMatrix(im)
for(v in names(im)){
  if(!anyNA(im[[v]]))meth[v] <- '' else if(is.factor(im[[v]]))meth[v] <- if(nlevels(im[[v]])==2) 'logreg' else 'polyreg' else meth[v] <- 'pmm'
}
stopifnot(all(meth[c('y','history','current_score')]==''))
cc <- d[complete.cases(d[,vars]),,drop=FALSE]
xx_pre <- model.matrix(full,cc)
stopifnot(qr(xx_pre)$rank==ncol(xx_pre))
preflight <- list(mode='reproduction_preflight',formula=paste(deparse(full),collapse=' '),seed=2026103001,m=30,maxit=20,observed_n=nrow(d),complete_case_n=nrow(cc),n_clusters=length(unique(d$cluster)),model_parameters=ncol(xx_pre),complete_case_design_rank=qr(xx_pre)$rank,age_knots=unname(age_knots),age_boundary=unname(age_boundary),imputation_methods=as.list(meth),imputation_predictors='All other model variables except self; outcome and categorical exposure included. Age enters imputation as observed continuous age. Household IDs are excluded.',reference_history='0',overall_test='mitml::testConstraints D1, 2 df, supplied clustered covariance matrices, df.com=G-p',scalar_pooling='Rubin total variance with Barnard-Rubin finite degrees of freedom; dfcom=G-p')
if(mode=='preflight'){
  write_json(preflight,file.path(resdir,'model_preflight.json'),pretty=TRUE,auto_unbox=TRUE)
  write.csv(pred,file.path(resdir,'imputation_predictor_matrix.csv'))
  writeLines(capture.output(sessionInfo()),file.path(resdir,'preflight_R_sessionInfo.txt'))
  cat(toJSON(preflight,pretty=TRUE,auto_unbox=TRUE),'\n');quit(save='no',status=0)
}
# This reproduction script does not imply preregistration or a prospective study freeze.
if(file.exists(file.path(resdir,'run_completion.json')))stop('Completed run already exists; do not overwrite without an explicit versioned rerun.')
warns <- character()
pool_scalar <- function(q,u,dfcom){
  m <- length(q); qb <- mean(q); ub <- mean(u); b <- if(m>1)var(q) else 0
  tt <- ub+(1+1/m)*b; se <- sqrt(tt); lambda <- if(tt>0)(1+1/m)*b/tt else 0
  if(m==1 || lambda<1e-12){df <- dfcom}else{
    old <- (m-1)/lambda^2; obs <- (dfcom+1)/(dfcom+3)*dfcom*(1-lambda); df <- 1/(1/old+1/obs)
  }
  crit <- qt(.975,df)
  c(estimate=qb,se=se,low=qb-crit*se,high=qb+crit*se,p=2*pt(-abs(qb/se),df),df=df,lambda=lambda,within=ub,between=b,total=tt,mcse=sqrt(b/m))
}
fit_one <- function(z,label,i){
  fit <- glm(full,data=z,family=poisson(link='log'),control=glm.control(epsilon=1e-12,maxit=100),na.action=na.fail)
  X <- model.matrix(fit); b <- coef(fit); V <- sandwich::vcovCL(fit,cluster=z$cluster,type='HC1',cadjust=TRUE)
  stopifnot(fit$converged,!anyNA(b),qr(X)$rank==ncol(X),all(is.finite(V)))
  jj <- match(c('history1','history2'),names(b)); stopifnot(!anyNA(jj))
  risk <- numeric(3); grad <- matrix(NA_real_,3,length(b)); max_cf <- 0; n_cf_over1 <- integer(3)
  for(g in 0:2){
    zg <- z; zg$history <- factor(g,levels=0:2)
    xg <- model.matrix(delete.response(terms(fit)),zg,contrasts.arg=fit$contrasts,xlev=fit$xlevels)
    pg <- as.vector(exp(xg %*% b)); risk[g+1] <- mean(pg); grad[g+1,] <- colMeans(xg*pg)
    max_cf <- max(max_cf,max(pg)); n_cf_over1[g+1] <- sum(pg>1)
  }
  rv <- grad %*% V %*% t(grad)
  if(i==1){
    write.csv(data.frame(y=z$y,cluster=z$cluster,X,check.names=FALSE),file.path(private,paste0(label,'_audit_matrix_private.csv')),row.names=FALSE)
    write.csv(data.frame(term=names(b),estimate=b,se=sqrt(diag(V))),file.path(resdir,paste0(label,'_audit_coefficients.csv')),row.names=FALSE)
    write.csv(V,file.path(resdir,paste0(label,'_audit_covariance.csv')))
  }
  list(q=b[jj],u=V[jj,jj,drop=FALSE],risk=risk,risk_u=rv,n=nrow(z),events=sum(z$y),G=length(unique(z$cluster)),p=length(b),dfcom=length(unique(z$cluster))-length(b),diagnostics=data.frame(analysis=label,imputation=i,n=nrow(z),events=sum(z$y),clusters=length(unique(z$cluster)),parameters=length(b),rank=qr(X)$rank,converged=fit$converged,iterations=fit$iter,condition_number=kappa(X,exact=TRUE),max_fitted=max(fitted(fit)),n_fitted_gt1=sum(fitted(fit)>1),max_counterfactual=max_cf,n_cf0_gt1=n_cf_over1[1],n_cf1_gt1=n_cf_over1[2],n_cf2_gt1=n_cf_over1[3]))
}
run_models <- function(sets,label){
  vals <- lapply(seq_along(sets),function(i)fit_one(sets[[i]],label,i))
  m <- length(vals); dfcom <- min(sapply(vals,`[[`,'dfcom'))
  Q <- sapply(vals,`[[`,'q'); U <- array(unlist(lapply(vals,`[[`,'u')),dim=c(2,2,m))
  dimnames(U) <- list(c('history1','history2'),c('history1','history2'),NULL)
  if(m>1){
    joint <- mitml::testConstraints(qhat=Q,uhat=U,constraints=c('history1','history2'),method='D1',df.com=dfcom)
    jointrow <- data.frame(analysis=label,method='D1 finite-df pooled Wald',m=m,F=joint$test[1,1],df1=joint$test[1,2],df2=joint$test[1,3],p=joint$test[1,4],riv=joint$test[1,5],dfcom=dfcom)
  }else{
    fv <- as.numeric(t(Q[,1]) %*% solve(U[,,1],Q[,1])/2)
    jointrow <- data.frame(analysis=label,method='Household-robust Wald F',m=m,F=fv,df1=2,df2=dfcom,p=pf(fv,2,dfcom,lower.tail=FALSE),riv=0,dfcom=dfcom)
  }
  contrasts <- list('1_vs_0'=c(1,0),'2_vs_0'=c(0,1),'2_vs_1_supportive'=c(-1,1))
  rr <- lapply(names(contrasts),function(nm){
    a <- contrasts[[nm]]
    qs <- as.vector(t(a)%*%Q); us <- sapply(seq_len(m),function(i)as.numeric(t(a)%*%U[,,i]%*%a))
    pooled <- pool_scalar(qs,us,dfcom)
    data.frame(analysis=label,contrast=nm,m=m,n=vals[[1]]$n,events=vals[[1]]$events,clusters=vals[[1]]$G,parameters=vals[[1]]$p,log_rr=unname(pooled['estimate']),se_log_rr=unname(pooled['se']),rr=exp(unname(pooled['estimate'])),ci_low=exp(unname(pooled['low'])),ci_high=exp(unname(pooled['high'])),p=unname(pooled['p']),df=unname(pooled['df']),lambda=unname(pooled['lambda']),mcse_log_rr=unname(pooled['mcse']))
  })
  risks <- sapply(vals,`[[`,'risk'); risks_u <- array(unlist(lapply(vals,`[[`,'risk_u')),dim=c(3,3,m))
  # Standardization is descriptive/conditional, not a causal intervention.
  rs <- lapply(0:2,function(g){
    z <- pool_scalar(risks[g+1,],sapply(seq_len(m),function(i)risks_u[g+1,g+1,i]),dfcom)
    data.frame(analysis=label,history=g,risk=unname(z['estimate']),low=unname(z['low']),high=unname(z['high']))
  })
  rd <- lapply(c('1_vs_0','2_vs_0'),function(nm){
    a <- if(nm=='1_vs_0')c(-1,1,0) else c(-1,0,1)
    z <- pool_scalar(as.vector(t(a)%*%risks),sapply(seq_len(m),function(i)as.numeric(t(a)%*%risks_u[,,i]%*%a)),dfcom)
    data.frame(analysis=label,contrast=nm,rd=unname(z['estimate']),low=unname(z['low']),high=unname(z['high']),p=unname(z['p']))
  })
  component_rows <- do.call(rbind,lapply(seq_len(m),function(i)data.frame(analysis=label,imputation=i,q1=Q[1,i],q2=Q[2,i],u11=U[1,1,i],u12=U[1,2,i],u22=U[2,2,i],dfcom=dfcom)))
  write.csv(component_rows,file.path(resdir,paste0(label,'_pooling_components.csv')),row.names=FALSE)
  list(rr=do.call(rbind,rr),joint=jointrow,risk=do.call(rbind,rs),rd=do.call(rbind,rd),diagnostics=do.call(rbind,lapply(vals,`[[`,'diagnostics')))
}
withCallingHandlers({
  cat('START MICE: 30 imputations x 20 iterations; seed 2026103001\n');flush.console()
  cache <- file.path(private,'primary_mi_private.rds')
  if(file.exists(cache)) {
    imp <- readRDS(cache)
    stopifnot(imp$m==30,imp$iteration==20,isTRUE(all.equal(imp$data,im,check.attributes=TRUE)))
  } else {
    imp <- mice(im,m=30,maxit=20,method=meth,predictorMatrix=pred,seed=2026103001,printFlag=FALSE)
    saveRDS(imp,cache)
  }
  event_count <- if(is.null(imp$loggedEvents))0L else nrow(imp$loggedEvents)
  if(event_count>0)write.csv(imp$loggedEvents,file.path(resdir,'mi_logged_events.csv'),row.names=FALSE)
  cv <- mice::convergence(imp)
  write.csv(cv,file.path(resdir,'mi_convergence.csv'),row.names=FALSE)
  cm <- imp$chainMean; chain <- list(); ci <- 0L
  for(v in dimnames(cm)[[1]])for(it in seq_len(dim(cm)[2]))for(ch in seq_len(dim(cm)[3])){
    value <- cm[v,it,ch]
    if(!is.na(value)){ci<-ci+1L;chain[[ci]]<-data.frame(variable=v,iteration=it,chain=ch,mean=value)}
  }
  if(length(chain))write.csv(do.call(rbind,chain),file.path(resdir,'mi_chain_means.csv'),row.names=FALSE)
  sets <- lapply(seq_len(30),function(i){z<-mice::complete(imp,i);z$cluster<-d$cluster;z})
  stopifnot(all(sapply(sets,function(z)!anyNA(z[,vars]))))
  for(z in sets)stopifnot(identical(z$y,d$y),identical(z$history,d$history),identical(z$current_score,d$current_score))
  cat('MICE complete. Fitting fixed primary model.\n');flush.console()
  primary <- run_models(sets,'primary')
  complete <- run_models(list(cc),'complete_case')
  outputs <- list(primary=primary,complete_case=complete)
  for(nm in names(primary))write.csv(do.call(rbind,lapply(outputs,`[[`,nm)),file.path(resdir,paste0(nm,'.csv')),row.names=FALSE)
  groups <- aggregate(cbind(n=rep(1,nrow(d)),events=d$y),by=list(history=d$history),FUN=sum)
  groups$crude_risk <- groups$events/groups$n
  write.csv(groups,file.path(resdir,'group_events.csv'),row.names=FALSE)
  writeLines(capture.output(sessionInfo()),file.path(resdir,'R_sessionInfo.txt'))
  writeLines(capture.output(print(primary$rr),print(primary$joint)),file.path(resdir,'main_result_console.txt'))
  write_json(list(status='completed',finished_at=as.character(Sys.time()),primary_imputations=30,complete_case_models=1,all_converged=all(primary$diagnostics$converged)&&all(complete$diagnostics$converged),mi_logged_events=event_count,warning_messages=warns),file.path(resdir,'run_completion.json'),pretty=TRUE,auto_unbox=TRUE)
  print(groups);print(primary$rr);print(primary$joint)
},warning=function(w){warns <<- c(warns,conditionMessage(w));cat('WARNING:',conditionMessage(w),'\n');invokeRestart('muffleWarning')})