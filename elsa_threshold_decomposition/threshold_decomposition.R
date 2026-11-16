# ELSA two-rule threshold decomposition. Post hoc plan: analysis_plan.md.
# Rscript threshold_decomposition.R <input.csv> <private_dir> <results_dir> [B=1000] [workers=10]
suppressPackageStartupMessages({library(mice); library(sandwich); library(splines); library(jsonlite); library(parallel)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)>=3)
script_path <- sub("^--file=", "", commandArgs()[grep("^--file=",commandArgs())][1])
source(file.path(dirname(script_path),"../workflow/path_guard.R"))
assert_external(args[2]); assert_external(args[3])
infile <- normalizePath(args[1],winslash="/",mustWork=TRUE)
private <- args[2]; out <- args[3]
B <- if(length(args)>=4) as.integer(args[4]) else 1000L
workers <- if(length(args)>=5) as.integer(args[5]) else 10L
for(p in c(private,out)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
raw <- read.csv(infile,stringsAsFactors=FALSE,colClasses=c(id="character",cluster="character"))
d <- raw[!is.na(raw$cesd_current)&raw$cesd_current<4&!is.na(raw$fall),,drop=FALSE]
stopifnot(nrow(d)==4779, !anyDuplicated(d$id),all(d$cesd_prior %in% 0:8),all(d$cesd_current %in% 0:3),all(d$fall %in% 0:1))
d$y <- d$fall; d$prior_score <- factor(d$cesd_prior,levels=0:8);d$current_score <- factor(d$cesd_current,levels=0:3)
for(v in c("education","smoking","window")) d[[v]] <- factor(d[[v]])
for(v in c("female","partnered","alcohol","stroke","hearte","diabe","arthre","stair_difficulty","adl_difficulty")) d[[v]] <- factor(d[[v]],levels=0:1)
vars <- c("y","prior_score","age","female","education","window","exposure_gap_months","followup_gap_months","current_score","partnered","smoking","alcohol","stroke","hearte","diabe","arthre","stair_difficulty","adl_difficulty")
im <- d[,vars,drop=FALSE]
cache <- file.path(private,"common_mi.rds")
if(file.exists(cache)) {imp<-readRDS(cache);stopifnot(identical(imp$data,im))} else {
 method<-mice::make.method(im);predictor<-mice::make.predictorMatrix(im)
 for(v in names(im)) method[v]<-if(!anyNA(im[[v]])) "" else if(is.factor(im[[v]])) {if(nlevels(im[[v]])==2) "logreg" else "polyreg"} else "pmm"
 cat("Starting common MICE: 30 imputations, 20 iterations\n");flush.console()
 imp<-mice(im,m=30,maxit=20,method=method,predictorMatrix=predictor,seed=20261101,printFlag=FALSE)
 saveRDS(imp,cache)
}
stopifnot(imp$m==30,imp$iteration==20)
write.csv(if(is.null(imp$loggedEvents)) data.frame() else imp$loggedEvents,file.path(out,"imputation_logged_events.csv"),row.names=FALSE)
if(!is.null(imp$loggedEvents)) stop("Inspect imputation logged events before proceeding")
conv <- tryCatch(mice::convergence(imp),error=function(e) data.frame(error=conditionMessage(e)))
write.csv(conv,file.path(out,"imputation_convergence.csv"),row.names=FALSE)
write.csv(data.frame(variable=names(im),n_missing=vapply(im,function(x)sum(is.na(x)),integer(1))),file.path(out,"missingness.csv"),row.names=FALSE)
age_basis <- ns(d$age,df=3)
age_knots <- attr(age_basis,"knots"); age_boundary <- attr(age_basis,"Boundary.knots")
groups <- sort(unique(d$cluster));G <- length(groups)
form <- y ~ hist + age1+age2+age3 + female + education + window + exposure_gap_months + followup_gap_months + current_score + partnered + smoking + alcohol + stroke + hearte + diabe + arthre + stair_difficulty + adl_difficulty
cell_rules <- data.frame(cell=c("A","B","C","D"),current_cut=c(4,3,4,3),history_cut=c(4,4,3,3))
sets <- lapply(1:30,function(i) {
 z<-complete(imp,i);z$cluster<-d$cluster;z$current_num<-d$cesd_current;z$prior_num<-d$cesd_prior
 for(j in 1:3)z[[paste0("age",j)]]<-age_basis[,j]
 z
})
make_design <- function(z,cc,hc) {
 z<-z[z$current_num<cc,,drop=FALSE];z$current_score<-droplevels(z$current_score)
 z$hist<-as.integer(z$prior_num>=hc)
 x<-model.matrix(form,z);y<-z$y
 f<-glm.fit(x,y,family=poisson("log"),control=glm.control(epsilon=1e-10,maxit=100))
 stopifnot(f$converged, !anyNA(f$coefficients))
 k<-which(colnames(x)=="hist");idx<-match(z$cluster,groups)
 list(x=x,y=y,idx=idx,k=k,beta=f$coefficients,q=unname(f$coefficients[k]),n=nrow(z),events=sum(y),exposed=sum(z$hist),clusters=length(unique(z$cluster)),p=ncol(x),max_fitted=max(f$fitted.values),n_fitted_gt1=sum(f$fitted.values>1),data=z)
}
all_designs <- lapply(sets,function(z) lapply(1:4,function(j)make_design(z,cell_rules$current_cut[j],cell_rules$history_cut[j])))
cat("Four-cell fits complete; starting paired household covariance estimation\n");flush.console()
# Analytic joint sandwich cross-check using stacked household influence functions.
analytic <- lapply(all_designs,function(ds) {
 inf<-matrix(0,nrow=G,ncol=4);diag_rows<-list()
 for(j in 1:4) {
  v<-ds[[j]];f<-glm(form,data=v$data,family=poisson("log"),control=glm.control(epsilon=1e-10,maxit=100))
  stopifnot(max(abs(coef(f)-v$beta))<1e-8)
  vc<-sandwich::vcovCL(f,cluster=v$data$cluster,type="HC1")
  mu<-f$fitted.values;bread_inv<-solve(crossprod(v$x,v$x*as.vector(mu)))
  infl_row<-as.vector((v$x*(v$y-mu))%*%bread_inv[,v$k])
  hi<-rowsum(infl_row,v$idx,reorder=FALSE);hidx<-as.integer(rownames(hi))
  correction<-(v$clusters/(v$clusters-1))*((v$n-1)/(v$n-v$p))
  estvar<-sum(hi^2)*correction
  stopifnot(abs(estvar-vc[v$k,v$k])<1e-8)
  inf[hidx,j]<-hi[,1]*sqrt(correction)
  diag_rows[[j]]<-data.frame(cell=cell_rules$cell[j],n=v$n,events=v$events,exposed=v$exposed,clusters=v$clusters,p=v$p,max_fitted=v$max_fitted,n_fitted_gt1=v$n_fitted_gt1,analytic_se=sqrt(estvar))
 }
 list(u=crossprod(inf),diag=do.call(rbind,diag_rows))
})
qmat<-do.call(rbind,lapply(all_designs,function(ds)vapply(ds,`[[`,numeric(1),"q")))
colnames(qmat)<-cell_rules$cell
write.csv(data.frame(imputation=1:30,qmat),file.path(out,"imputation_estimates.csv"),row.names=FALSE)
write.csv(do.call(rbind,lapply(1:30,function(m)cbind(imputation=m,analytic[[m]]$diag))),file.path(out,"fit_diagnostics.csv"),row.names=FALSE)
# Fixed designs and household frequency weights give exactly the same coefficient as replicated rows.
checkv<-all_designs[[1]][[1]];set.seed(831109);hh<-sample.int(G,G,replace=TRUE);freq<-tabulate(hh,nbins=G);ww<-freq[checkv$idx]
f_w<-glm.fit(checkv$x,checkv$y,weights=ww,family=poisson("log"),start=checkv$beta,control=glm.control(epsilon=1e-10,maxit=100))
rr<-rep(seq_along(ww),ww);f_r<-glm.fit(checkv$x[rr,,drop=FALSE],checkv$y[rr],family=poisson("log"),control=glm.control(epsilon=1e-10,maxit=100))
stopifnot(max(abs(f_w$coefficients-f_r$coefficients))<1e-8)
boot_worker <- function(m) {
 dest<-file.path(private,sprintf("bootstrap_m%02d_B%d.rds",m,B))
 if(file.exists(dest)) return(readRDS(dest))
 ds<-all_designs[[m]];set.seed(20261101+1000L*m)
 boot<-matrix(NA_real_,nrow=B,ncol=4);messages<-list();starttime<-Sys.time()
 for(b in seq_len(B)) {
  freq<-tabulate(sample.int(G,G,replace=TRUE),nbins=G)
  for(j in 1:4) {
   v<-ds[[j]];ww<-freq[v$idx]
   f<-tryCatch(withCallingHandlers(glm.fit(v$x,v$y,weights=ww,family=poisson("log"),start=v$beta,control=glm.control(epsilon=1e-10,maxit=100)),warning=function(w){messages[[length(messages)+1L]]<<-data.frame(imputation=m,replicate=b,cell=j,type="warning",message=conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e)e)
   if(inherits(f,"error") || !isTRUE(f$converged) || anyNA(f$coefficients)) {
    messages[[length(messages)+1L]]<-data.frame(imputation=m,replicate=b,cell=j,type="failure",message=if(inherits(f,"error"))conditionMessage(f) else "nonconverged or singular")
   } else boot[b,j]<-f$coefficients[v$k]
  }
 }
 obj<-list(m=m,boot=boot,u=if(anyNA(boot))matrix(NA_real_,4,4) else cov(boot),messages=if(length(messages))do.call(rbind,messages) else data.frame(),elapsed=as.numeric(difftime(Sys.time(),starttime,units="secs")))
 saveRDS(obj,dest);obj
}
cl<-makeCluster(min(workers,30))
clusterExport(cl,c("private","B","G","all_designs","boot_worker"),envir=environment())
bootres<-tryCatch(parLapplyLB(cl,1:30,boot_worker),finally=stopCluster(cl))
messages<-do.call(rbind,lapply(bootres,`[[`,"messages"));write.csv(messages,file.path(out,"bootstrap_messages.csv"),row.names=FALSE)
write.csv(data.frame(imputation=1:30,replicates=B,failures=vapply(bootres,function(v)sum(!complete.cases(v$boot)),integer(1)),seconds=vapply(bootres,`[[`,numeric(1),"elapsed")),file.path(out,"bootstrap_diagnostics.csv"),row.names=FALSE)
stopifnot(all(vapply(bootres,function(v)all(is.finite(v$boot)),logical(1))))
contrast<-rbind(A=c(1,0,0,0),B=c(0,1,0,0),C=c(0,0,1,0),D=c(0,0,0,1),total=c(-1,0,0,1),eligibility=c(-.5,.5,-.5,.5),history_definition=c(-.5,-.5,.5,.5),path_interaction=c(1,-1,-1,1),eligibility_first=c(-1,1,0,0),history_after_eligibility=c(0,-1,0,1),history_first=c(-1,0,1,0),eligibility_after_history=c(0,0,-1,1))
stopifnot(all(contrast["eligibility",]+contrast["history_definition",]==contrast["total",]))
dfcom<-min(vapply(all_designs[[1]],function(v)v$clusters-v$p,numeric(1)))
pool_joint<-function(ulist,label) {
 W<-Reduce(`+`,ulist)/30;between<-cov(qmat);T<-W+(1+1/30)*between;point<-colMeans(qmat)
 stopifnot(min(eigen(T,symmetric=TRUE,only.values=TRUE)$values)>0)
 rows<-lapply(1:nrow(contrast),function(k) {
  l<-contrast[k,];q<-sum(l*point);u<-as.numeric(t(l)%*%W%*%l);bv<-as.numeric(t(l)%*%between%*%l);tv<-u+(1+1/30)*bv;se<-sqrt(tv);lambda<-(1+1/30)*bv/tv
  df<-if(lambda<1e-12)dfcom else 1/(lambda^2/29+1/(((dfcom+1)/(dfcom+3))*dfcom*(1-lambda)))
  ci<-q+c(-1,1)*qt(.975,df)*se
  data.frame(term=rownames(contrast)[k],method=label,estimate=q,se=se,ci_low=ci[1],ci_high=ci[2],exp_estimate=exp(q),exp_low=exp(ci[1]),exp_high=exp(ci[2]),df=df,lambda=lambda,within_variance=u,between_variance=bv)
 })
 for(nm in c("W","between","T"))write.csv(get(nm),file.path(out,paste0(label,"_covariance_",nm,".csv")),row.names=FALSE)
 do.call(rbind,rows)
}
pooled<-rbind(pool_joint(lapply(bootres,`[[`,"u"),"MI_Boot"),pool_joint(lapply(analytic,`[[`,"u"),"joint_HC1"))
write.csv(pooled,file.path(out,"pooled_contrasts.csv"),row.names=FALSE)
cells<-merge(cell_rules,analytic[[1]]$diag,by="cell",sort=FALSE)
cells<-merge(cells,pooled[pooled$method=="MI_Boot"&pooled$term%in%cell_rules$cell,],by.x="cell",by.y="term",sort=FALSE)
write.csv(cells,file.path(out,"four_cell_results.csv"),row.names=FALSE)
# Before/after reclassification counts with observed outcome, no individual-level export.
prior_group<-ifelse(d$cesd_prior<3,"0-2",ifelse(d$cesd_prior==3,"3","4-8"));current_group<-ifelse(d$cesd_current<3,"0-2","3")
counts<-aggregate(cbind(n=rep(1,nrow(d)),events=d$y),list(current=current_group,prior=prior_group),sum);counts$risk<-counts$events/counts$n
write.csv(counts,file.path(out,"observed_score_cells.csv"),row.names=FALSE)
checks<-list(date=as.character(Sys.Date()),input_md5=unname(tools::md5sum(infile)),n_broad=nrow(d),households=G,n_narrow=sum(d$cesd_current<3),excluded_current3=sum(d$cesd_current==3),reclassified_prior3_broad=sum(d$cesd_prior==3),reclassified_prior3_narrow=sum(d$cesd_prior==3&d$cesd_current<3),imputations=30,iterations=20,bootstrap_per_imputation=B,bootstrap_fits=30*B*4,bootstrap_failed=0,spline_knots=age_knots,spline_boundaries=age_boundary,dfcom=dfcom,weighted_vs_replicated_max_abs=max(abs(f_w$coefficients-f_r$coefficients)),component_identity_max_abs=max(abs(qmat%*%contrast["eligibility",]+qmat%*%contrast["history_definition",]-qmat%*%contrast["total",])))
write_json(checks,file.path(out,"validation.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
print(pooled[pooled$method=="MI_Boot",c("term","estimate","ci_low","ci_high","exp_estimate","exp_low","exp_high")]);cat("DONE\n")
