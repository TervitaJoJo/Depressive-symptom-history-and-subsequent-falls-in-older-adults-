args <- commandArgs(trailingOnly=TRUE)
script_path <- sub("^--file=", "", commandArgs()[grep("^--file=",commandArgs())][1])
source(file.path(dirname(script_path),"../workflow/path_guard.R"))
if(length(args)>=3) {assert_external(args[2]); assert_external(args[3])}

if(length(args)<3)stop('Usage: Rscript targeted_analyses.R EXTENDED_INPUT OUTPUT FROZEN_CACHE')
indir<-args[1];outdir<-args[2];cachedir<-args[3];dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
scriptfile<-sub('^--file=','',commandArgs()[grep('^--file=',commandArgs())]);source(file.path(dirname(scriptfile),'model_helpers.R'))
dir.create(cachedir,recursive=TRUE,showWarnings=FALSE)
newcache<-file.path(outdir,'private_cache');dir.create(newcache,showWarnings=FALSE)
rows<-list();diags<-list();scores<-list();curves<-list();tests<-list();evs<-list();flow<-list();meds<-list();domainmiss<-list();knotslog<-list()
add <- function(o,f,label) {
 a<-pool_model(o,f,label);rows[[length(rows)+1]]<<-a$row;diags[[length(diags)+1]]<<-a$diag
 cat(o$cohort,label,'n=',a$row$n,'RR=',sprintf('%.4f [%.4f,%.4f]',a$row$rr,a$row$ci_low,a$row$ci_high),'\n');flush.console();a
}
subset_obj<-function(o,mask){z<-o;z$sets<-lapply(o$sets,function(x)x[mask,,drop=FALSE]);z}
fmt<-function(d,co,cut){d$y<-d$fall;d$hist<-factor(as.integer(d$cesd_prior>=cut),levels=0:1);d$current_score<-factor(d$cesd_current,levels=0:(cut-1));for(v in c('education','smoking','window'))d[[v]]<-factor(d[[v]]);for(v in c('female','partnered','alcohol','stroke','hearte','diabe','arthre','stair_difficulty','adl_difficulty'))d[[v]]<-factor(d[[v]],levels=0:1);d}
new_mi<-function(d,form,label,co,extras=character(0)){
 im<-d[,unique(c(all.vars(form),extras)),drop=FALSE];cache<-file.path(newcache,paste0(co,'_',label,'.rds'))
 if(file.exists(cache)){imp<-readRDS(cache);stopifnot(isTRUE(all.equal(imp$data,im)))}else{
  method<-mice::make.method(im);pred<-mice::make.predictorMatrix(im)
  for(v in names(im)){if(!anyNA(im[[v]]))method[v]<-'' else if(is.factor(im[[v]]))method[v]<-if(nlevels(im[[v]])==2)'logreg' else 'polyreg' else method[v]<-'pmm'}
  imp<-mice(im,m=30,maxit=20,method=method,predictorMatrix=pred,seed=20261031+match(co,c('CHARLS','ELSA'))+ncol(im),printFlag=FALSE);saveRDS(imp,cache)
 }
 write.csv(imp$loggedEvents,file.path(outdir,paste0(co,'_',label,'_imputation_events.csv')),row.names=FALSE)
 list(cohort=co,sets=lapply(1:imp$m,function(i){z<-complete(imp,i);z$cluster<-as.character(d$cluster);z}),d=d,imp=imp)
}
rcs3<-function(x,k){(pmax(x-k[1],0)^3-pmax(x-k[2],0)^3*(k[3]-k[1])/(k[3]-k[2])+pmax(x-k[3],0)^3*(k[2]-k[1])/(k[3]-k[2]))/(k[3]-k[1])^2}
pooled_test<-function(vals,terms,co,label){
 Q<-do.call(cbind,lapply(vals,function(v)coef(v$fit)[terms]));U<-array(NA_real_,c(length(terms),length(terms),length(vals)));for(i in seq_along(vals))U[,,i]<-vals[[i]]$vv[terms,terms,drop=FALSE]
 tt<-mitml::testConstraints(qhat=Q,uhat=U,constraints=terms,method='D1',df.com=min(vapply(vals,`[[`,numeric(1),'dfcom')))$test[1,]
 data.frame(cohort=co,test=label,F=unname(tt[1]),df1=unname(tt[2]),df2=unname(tt[3]),p=unname(tt[4]))
}
for(co in c('CHARLS','ELSA')){
 cut<-if(co=='CHARLS')10 else 4
 csv<-file.path(indir,paste0(co,'_extended.csv'));raw<-read.csv(csv,stringsAsFactors=FALSE)
 o<-load_cohort(co,cut,csv,file.path(cachedir,paste0(co,'_fall_',cut,'_mi.rds')))
 d<-o$d;stopifnot(nrow(d)==if(co=='CHARLS')2475 else 4779)
 a<-add(o,primary_full,'original_primary_check');ev<-function(rr)if(rr<=1)1 else rr+sqrt(rr*(rr-1))
 evs[[co]]<-data.frame(cohort=co,rr=a$row$rr,ci_low=a$row$ci_low,ci_high=a$row$ci_high,evalue_point=ev(a$row$rr),evalue_ci=ev(a$row$ci_low))
 if(co=='CHARLS'){
  add(subset_obj(o,d$antidepressant_report==0),primary_full,'exclude_reported_antidepressants')
  add(subset_obj(o,d$antidepressant_report==0 & d$sedative_report==0),primary_full,'exclude_reported_antidepressants_or_sedatives')
 }else{
  eligible<-!is.na(d$elsa_prior_med_observed)&d$elsa_prior_med_observed==1
  add(subset_obj(o,eligible),primary_full,'prior_nurse_medication_subset')
  add(subset_obj(o,eligible & !is.na(d$elsa_prior_antidepressant) & d$elsa_prior_antidepressant==0),primary_full,'exclude_prior_nurse_antidepressants')
 }
 domains<-c('recall_total','poor_eyesight','log_walk_time',if(co=='CHARLS')'balance_score')
 ok<-complete.cases(d[,domains])
 formplus<-update(primary_full,as.formula(paste('. ~ . +',paste(domains,collapse=' + '))))
 dd<-fmt(d[ok,,drop=FALSE],co,cut);dd$poor_eyesight<-factor(dd$poor_eyesight,levels=0:1)
 if(co=='CHARLS')dd$balance_score<-factor(dd$balance_score,levels=1:4)
 # Include all added measured domains in the imputation model; compare on identical sets.
 aug<-new_mi(dd,formplus,'measured_domains',co)
 add(aug,primary_full,'measured_domain_subset_base')
 add(aug,formplus,'plus_measured_domains')
 domainmiss[[co]]<-data.frame(cohort=co,variable=domains,missing=vapply(domains,function(v)sum(is.na(d[[v]])),integer(1)),n=nrow(d))
 # Continuous historical scores, with observed nonlinear terms included in imputation.
 fd<-fmt(d,co,cut);fd$prior_score<-fd$cesd_prior;k<-as.numeric(quantile(fd$prior_score,c(.1,.5,.9),type=7));if(length(unique(k))<3){stopifnot(co=='ELSA');k<-c(0,2,5)}
 fd$prior_rcs<-rcs3(fd$prior_score,k);fs<-update(primary_full,.~.-hist+prior_score+prior_rcs);fl<-update(fs,.~.-prior_rcs)
 so<-new_mi(fd,fs,'score_shape',co);knotslog[[co]]<-data.frame(cohort=co,k1=k[1],k2=k[2],k3=k[3],reference=0)
 vl<-lapply(so$sets,function(z)fit_one(z,fl));v<-lapply(so$sets,function(z)fit_one(z,fs))
 stopifnot(all(vapply(c(vl,v),function(a)a$converged,logical(1))))
 for(label in c('prior_score_linear','prior_score_spline')){
  vv<-if(label=='prior_score_linear')vl else v
  diags[[length(diags)+1]]<-data.frame(cohort=co,model=label,imputation=seq_along(vv),n=vapply(vv,`[[`,numeric(1),'n'),events=vapply(vv,`[[`,numeric(1),'events'),clusters=vapply(vv,`[[`,numeric(1),'clusters'),converged=vapply(vv,`[[`,logical(1),'converged'),max_fitted=vapply(vv,`[[`,numeric(1),'max_fitted'),n_fitted_gt1=vapply(vv,`[[`,numeric(1),'n_fitted_gt1'),weighting='unweighted')
 }
 step<-if(co=='CHARLS')5 else 1
 pl<-pool_scalar(vapply(vl,function(a)coef(a$fit)['prior_score']*step,numeric(1)),vapply(vl,function(a)a$vv['prior_score','prior_score']*step^2,numeric(1)),min(vapply(vl,`[[`,numeric(1),'dfcom')))
 scores[[co]]<-data.frame(cohort=co,n=nrow(d),events=sum(d$y),increment=step,rr=exp(pl['est']),ci_low=exp(pl['lo']),ci_high=exp(pl['hi']),p=pl['p'])
 tests[[paste0(co,'overall')]]<-pooled_test(v,c('prior_score','prior_rcs'),co,'overall_2df');tests[[paste0(co,'nonlinear')]]<-pooled_test(v,'prior_rcs',co,'nonlinear_1df')
 upper<-as.numeric(quantile(fd$prior_score,.99));grid<-seq(0,upper,length.out=151)
 curves[[co]]<-do.call(rbind,lapply(grid,function(x){dx<-c(x,rcs3(x,k)-rcs3(0,k));tt<-pool_scalar(vapply(v,function(a)sum(dx*coef(a$fit)[c('prior_score','prior_rcs')]),numeric(1)),vapply(v,function(a)as.numeric(t(dx)%*%a$vv[c('prior_score','prior_rcs'),c('prior_score','prior_rcs')]%*%dx),numeric(1)),min(vapply(v,`[[`,numeric(1),'dfcom')));data.frame(cohort=co,score=x,rr=exp(tt['est']),ci_low=exp(tt['lo']),ci_high=exp(tt['hi']))}))
 write.csv(as.data.frame(table(score=fd$prior_score)),file.path(outdir,paste0(co,'_score_counts.csv')),row.names=FALSE)
 # Record target-wave mortality and outcome ascertainment among all current-low participants.
 low<-raw[raw$cesd_current<cut,,drop=FALSE];low$history<-as.integer(low$cesd_prior>=cut)
 low$status<-ifelse(!is.na(low$fall),'fall_response',ifelse(low$recorded_death==1,'recorded_death',ifelse(low$target_status%in%c(1,4),'known_alive_no_valid_fall_response','unknown_vital_or_interview_status')))
 flow[[co]]<-aggregate(rep(1,nrow(low)),list(cohort=rep(co,nrow(low)),history=low$history,status=low$status),sum);names(flow[[co]])[4]<-'n'
 fdeath<-update(primary_full,.~.-followup_gap_months);add(o,fdeath,'falls_only_no_future_duration')
 cd<-low[!is.na(low$fall)|low$recorded_death==1,,drop=FALSE];cd$fall[cd$recorded_death==1]<-1;cd<-fmt(cd,co,cut)
 mo<-new_mi(cd,fdeath,'fall_or_recorded_death',co);add(mo,fdeath,'fall_or_recorded_death')
 write.csv(do.call(rbind,rows),file.path(outdir,'targeted_models.csv'),row.names=FALSE)
 write.csv(do.call(rbind,evs),file.path(outdir,'evalues.csv'),row.names=FALSE)
 write.csv(do.call(rbind,scores),file.path(outdir,'continuous_score_models.csv'),row.names=FALSE)
 write.csv(do.call(rbind,tests),file.path(outdir,'spline_tests.csv'),row.names=FALSE)
 write.csv(do.call(rbind,curves),file.path(outdir,'spline_curves.csv'),row.names=FALSE)
 write.csv(do.call(rbind,flow),file.path(outdir,'mortality_ascertainment.csv'),row.names=FALSE)
 write.csv(do.call(rbind,domainmiss),file.path(outdir,'domain_missingness.csv'),row.names=FALSE)
 write.csv(do.call(rbind,knotslog),file.path(outdir,'spline_knots.csv'),row.names=FALSE)
 write.csv(do.call(rbind,diags),file.path(outdir,'targeted_diagnostics.csv'),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(outdir,'R_session.txt'))
cat('DONE\n')
