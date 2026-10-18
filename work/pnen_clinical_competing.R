library(data.table);library(survival);library(splines)
base<-Sys.getenv('PNEN_PROJECT',getwd());setwd(base)
src<-readLines('work/pnen_clinical_standardization.R')
stopat<-grep('^estimates<-list',src)[1]-1L
prefix<-paste(src[seq_len(stopat)],collapse='\n')
prefix<-gsub("'outputs/pnen_clinical'","'outputs/pnen_clinical_competing'",prefix,fixed=TRUE)
eval(parse(text=prefix));set.seed(20261017);B<-1000L
out<-'outputs/pnen_clinical_competing';dir.create(out,recursive=TRUE,showWarnings=FALSE)
ajsrc<-readLines('work/pnen_formal_competing.R')
code<-paste(ajsrc[grep('^calcAJ<-function',ajsrc):(grep('^results<-list',ajsrc)[1]-1)],collapse='\n')
code<-sub('fit<-glm.fit(ca$X[ok,,drop=FALSE],cbind(cases[ok],count[ok]-cases[ok]),family=binomial())',"X<-ca$X; q<-qr(X[ok,,drop=FALSE]); X<-X[,sort(q$pivot[seq_len(q$rank)]),drop=FALSE]; fit<-glm.fit(X[ok,,drop=FALSE],cbind(cases[ok],count[ok]-cases[ok]),family=binomial(),control=glm.control(maxit=100))",code,fixed=TRUE)
code<-sub('ca$X%*%fit$coefficients','X%*%fit$coefficients',code,fixed=TRUE)
code<-sub('if(any(ps<=0|ps>=1))','if(any((ps<=0|ps>=1) & mult>0))',code,fixed=TRUE)
eval(parse(text=code));aj_impl<-calcAJ
warn<-list();fail<-list();est<-list();direct<-list();events<-list();stage_now<-'';model_now<-'';L_now<-0;rep_now<-0
calcAJ<-function(...)withCallingHandlers(aj_impl(...),warning=function(w){warn[[length(warn)+1L]]<<-data.table(stage=stage_now,model=model_now,landmark=L_now,replicate=rep_now,message=conditionMessage(w));invokeRestart('muffleWarning')})
for(st in c('Localized','Regional')) {
 stage_now<-st;z<-copy(d[stage3==st]);z[,extent_group:=factor(extent_group)];z[,size_group:=droplevels(size_group)];n<-nrow(z)
 caches<-list();points<-list()
 for(mm in names(forms)) {
  model_now<-mm;model_formula<-forms[[mm]];covars<-all.vars(model_formula);rep_now<-0
  caches[[mm]]<-lapply(LM,function(L){L_now<<-L;make_cache(z,L)})
  points[[mm]]<-lapply(seq_along(LM),function(j){L_now<<-LM[j];calcAJ(caches[[mm]][[j]],rep(1,n),TRUE)})
 }
 for(L in LM){a<-z[days>L*365.25 | (L==0 & days==0)];a[,cause:=fcase(event==0,0L,cancer==1,1L,other==1,2L,default=3L)];events[[length(events)+1L]]<-a[,.(n=.N,cancer_events=sum(event==1 & cause==1 & days<=(L+3)*365.25),other_events=sum(event==1 & cause==2 & days<=(L+3)*365.25),unknown_events=sum(event==1 & cause==3 & days<=(L+3)*365.25),endpoint_n=sum(days>=(L+3)*365.25)),by=g][,`:=`(stage=st,landmark=L)]}
 bb<-array(NA_real_,c(B,3,5,3,2,4),dimnames=list(NULL,names(forms),as.character(LM),methods,c('SCLC','LCNEC'),c('OS','cancer','other','unknown')))
 i0<-which(z$g==0);i1<-which(z$g==1)
 for(b in seq_len(B)) {
  rep_now<-b;m<-integer(n);m[i0]<-tabulate(sample.int(length(i0),length(i0),TRUE),length(i0));m[i1]<-tabulate(sample.int(length(i1),length(i1),TRUE),length(i1))
  for(mm in names(forms))for(j in seq_along(LM)) {
   model_now<-mm;L_now<-LM[j];ans<-tryCatch(calcAJ(caches[[mm]][[j]],m),error=function(e)e)
   if(inherits(ans,'error'))fail[[length(fail)+1L]]<-data.table(stage=st,model=mm,landmark=LM[j],replicate=b,error=conditionMessage(ans)) else bb[b,mm,j,,,]<-ans
  }
  if(b%%100==0)cat(st,b,'/',B,'\n')
 }
 # Identical patient resampling and weights must reproduce every prior clinical OS replicate.
 oldboot<-readRDS(file.path('outputs/pnen_clinical',paste0('bootstrap_',st,'.rds')))
 stopifnot(max(abs(bb[,,,,,'OS']-oldboot[,,,,1:2]))<1e-8)
 for(mm in names(forms))for(j in seq_along(LM))for(met in methods){
  for(oc in c('OS','cancer','other','unknown')) {
   b0<-bb[,mm,j,met,'SCLC',oc];b1<-bb[,mm,j,met,'LCNEC',oc];bd<-b1-b0;ci<-quantile(bd,c(.025,.975))
   pt<-points[[mm]][[j]][met,,oc]
   est[[length(est)+1L]]<-data.table(stage=st,model=mm,landmark=LM[j],method=met,outcome=oc,SCLC=pt[1],SCLC_lower=quantile(b0,.025),SCLC_upper=quantile(b0,.975),LCNEC=pt[2],LCNEC_lower=quantile(b1,.025),LCNEC_upper=quantile(b1,.975),difference=pt[2]-pt[1],lower=ci[1],upper=ci[2],valid_B=sum(is.finite(bd)))
  }
  pt<-points[[mm]][[j]][met,,];dd<-pt['LCNEC',]-pt['SCLC',]
  diffs<-bb[,mm,j,met,'LCNEC',]-bb[,mm,j,met,'SCLC',];v<-diffs[,'other']-diffs[,'cancer']
  direct[[length(direct)+1L]]<-data.table(stage=st,model=mm,landmark=LM[j],method=met,contrast='Other-cause difference minus cancer-related difference',estimate=dd['other']-dd['cancer'],lower=quantile(v,.025),upper=quantile(v,.975),valid_B=sum(is.finite(v)))
 }
 stopifnot(max(abs(apply(bb,c(1,2,3,4,5),sum)-1))<1e-8)
 saveRDS(bb,file.path(out,paste0('bootstrap_',st,'.rds')))
}
rr<-rbindlist(est);dc<-rbindlist(direct);fwrite(rr,file.path(out,'clinical_competing_estimates.csv'));fwrite(dc,file.path(out,'clinical_component_direct_contrasts.csv'));fwrite(rbindlist(events),file.path(out,'competing_horizon_support.csv'))
if(length(warn))fwrite(rbindlist(warn),file.path(out,'clinical_competing_warnings.csv'))
if(length(fail))fwrite(rbindlist(fail),file.path(out,'clinical_competing_failures.csv')) else writeLines('No failed bootstrap estimation cells.',file.path(out,'clinical_competing_failures.txt'))
stopifnot(all(rr$valid_B==B),all(dc$valid_B==B))
prior<-fread('outputs/pnen_formal_competing/conditional_competing_differences.csv')
prior[,outcome:=sub('_CIF$','',outcome)]
ck<-merge(rr[model=='C0'],prior[stage %in% c('Localized','Regional')],by=c('stage','landmark','method','outcome'),suffixes=c('_new','_old'))
stopifnot(nrow(ck)==120,max(abs(ck$difference_new-ck$difference_old))<1e-8)
writeLines(c('Clinical competing-risk extension, 2026-10-18. Seed20261017;1000 joint bootstrap samples shared with the prior clinical OS extension.', 'C0 demographic; C1 adds grouped AJCC6 stage; C2 adds size groups. Same cohort, landmarks and weights as clinical OS analysis.', 'All bootstrap OS group probabilities exactly reproduced the saved clinical OS bootstrap within1e-8. C0 point estimates reproduced original competing-risk results within1e-8.', 'All point-estimate Aalen-Johansen state probabilities independently matched survival::survfit within1e-8. OS+cancer+other+unknown=1 in every bootstrap cell.', 'Direct signed contrast: (LCNEC-SCLC other-cause CIF) minus (LCNEC-SCLC cancer-related CIF), using the same bootstrap replicate. All component and direct intervals are exploratory pointwise intervals without study-wide multiplicity control.', 'No data-selected endpoint suppression, trimming or significance-based choice. Full support/weight/balance outputs retained; no clinical-subtype causal attribution.',capture.output(sessionInfo())),file.path(out,'clinical_competing_validation_and_session.txt'))
print(rr[stage=='Regional' & landmark==5 & model=='C2' & method %in% c('IPTW','overlap')]);print(dc[stage=='Regional' & landmark==5 & method %in% c('IPTW','overlap')])
