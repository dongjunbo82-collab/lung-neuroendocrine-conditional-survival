library(data.table); library(survival); library(splines)
setDTthreads(1); set.seed(20261010); B<-1000L
out<-'outputs/pnen_extension'
d<-readRDS('work/pnen_extension_data.rds')
# Reuse the independently validated weighted product-limit implementation.
src<-readLines('work/pnen_formal_OS.R')
start<-grep('^make_cache<-function',src); end<-grep('^allres<-list',src)-1L
code<-paste(src[start:end],collapse='\n')
code<-sub('cv<-unique(a[,.(age,sex,year)])','cv<-unique(a[,..covars])',code,fixed=TRUE)
code<-sub('idx<-cv[a,on=.(age,sex,year),which=TRUE]','idx<-cv[a,on=covars,which=TRUE]',code,fixed=TRUE)
old<-"~ns(age,knots=c(50,70),Boundary.knots=c(18,90))+factor(sex)+ns(year,knots=c(2008,2012),Boundary.knots=c(2004,2015))"
code<-sub(old,'model_formula',code,fixed=TRUE)
eval(parse(text=code)); methods<-c('crude','IPTW','overlap')
calc_impl<-calc; warning_log<-list();active_job<-'';active_rep<-0L
calc<-function(...)withCallingHandlers(calc_impl(...),warning=function(w){warning_log[[length(warning_log)+1L]]<<-data.table(job=active_job,replicate=active_rep,message=conditionMessage(w));invokeRestart('muffleWarning')})
formula0<-as.formula(old)
models<-list(M0=formula0,M1=update(formula0,~.+factor(surgery)),M2=update(formula0,~.+factor(surgery)+rt_yes+chemo))
res<-list(); diagnostics<-list(); balances<-list(); contrasts<-list(); failure<-list(); bootstore<-list()
fitjob<-function(z,landmarks,forms,label) {
 active_job<<-label;active_rep<<-0L
 n<-nrow(z); caches<-list()
 for(model in names(forms)) {
  model_formula<<-forms[[model]];covars<<-all.vars(model_formula)
  caches[[model]]<-lapply(landmarks,function(L)make_cache(z,L))
 }
 base<-lapply(caches,function(cs)lapply(cs,function(ca)calc(ca,rep(1,n),TRUE)))
 boots<-array(NA_real_,c(B,length(forms),length(landmarks),3,3),dimnames=list(NULL,names(forms),as.character(landmarks),methods,c('SCLC','LCNEC','difference')))
 i0<-which(z$g==0);i1<-which(z$g==1);errors<-character()
 for(b in seq_len(B)) {
  active_rep<<-b
  m<-integer(n);m[i0]<-tabulate(sample.int(length(i0),length(i0),TRUE),length(i0));m[i1]<-tabulate(sample.int(length(i1),length(i1),TRUE),length(i1))
  for(k in seq_along(forms)) for(j in seq_along(landmarks)) {
   ans<-tryCatch(calc(caches[[k]][[j]],m),error=function(e)e)
   if(inherits(ans,'error')) errors<-c(errors,paste(b,k,j,conditionMessage(ans))) else boots[b,k,j,,]<-ans
  }
  if(b%%250==0) cat(label,b,'/',B,'\n')
 }
 for(k in seq_along(forms)) for(j in seq_along(landmarks)) {
  for(method in methods) {
   bb<-boots[,k,j,method,'difference'];ci<-quantile(bb,c(.025,.975),na.rm=TRUE)
   res[[length(res)+1L]]<<-data.table(job=label,model=names(forms)[k],landmark=landmarks[j],method=method,SCLC=base[[k]][[j]]$result[method,1],LCNEC=base[[k]][[j]]$result[method,2],difference=base[[k]][[j]]$result[method,3],lower=ci[1],upper=ci[2],valid_B=sum(is.finite(bb)))
  }
  diagnostics[[length(diagnostics)+1L]]<<-copy(base[[k]][[j]]$diag)[,`:=`(job=label,model=names(forms)[k],landmark=landmarks[j])]
  balances[[length(balances)+1L]]<<-copy(base[[k]][[j]]$balance)[,`:=`(job=label,model=names(forms)[k],landmark=landmarks[j])]
 }
 failure[[length(failure)+1L]]<<-data.table(job=label,failed_cells=length(errors),total_cells=B*length(forms)*length(landmarks))
 if(length(errors))writeLines(errors,file.path(out,paste0('failures_',label,'.txt')))
 saveRDS(boots,file.path(out,paste0('bootstrap_',label,'.rds')))
 bootstore[[label]]<<-boots
 invisible(NULL)
}
for(st in c('Localized','Regional')) fitjob(d[year<=2015 & stage3==st & surgery!='Unknown'],c(0,1,2,3,5),models,paste0('Treatment_',st))
summarize_family<-function(est,bb,labels,family,method){
 ok<-apply(bb,1,function(x)all(is.finite(x))); bb<-bb[ok,,drop=FALSE]
 se<-apply(bb,2,sd); crit<-quantile(apply(abs(sweep(sweep(bb,2,est,'-'),2,se,'/')),1,max),.95)
 for(j in seq_along(est)) contrasts[[length(contrasts)+1L]]<<-data.table(family=family,method=method,contrast=labels[j],estimate=est[j],lower=quantile(bb[,j],.025),upper=quantile(bb[,j],.975),sim_lower=est[j]-crit*se[j],sim_upper=est[j]+crit*se[j],valid_B=nrow(bb))
}
rr<-rbindlist(res)
for(met in methods){
 est<-sapply(names(models),function(mm) {x<-rr[method==met & model==mm];(x[job=='Treatment_Localized' & landmark==5,difference]-x[job=='Treatment_Localized' & landmark==0,difference])-(x[job=='Treatment_Regional' & landmark==5,difference]-x[job=='Treatment_Regional' & landmark==0,difference])})
 bb<-sapply(names(models),function(mm){a<-bootstore[['Treatment_Localized']];b<-bootstore[['Treatment_Regional']];a[,mm,'5',met,3]-a[,mm,'0',met,3]-b[,mm,'5',met,3]+b[,mm,'0',met,3]})
 summarize_family(est,bb,names(models),'Treatment stage change 0 to 5',met)
 for(k in 2:3){v<-bb[,k]-bb[,1];contrasts[[length(contrasts)+1L]]<-data.table(family='Paired treatment model change',method=met,contrast=paste0(names(models)[k],'-M0'),estimate=est[k]-est[1],lower=quantile(v,.025,na.rm=TRUE),upper=quantile(v,.975,na.rm=TRUE),sim_lower=NA_real_,sim_upper=NA_real_,valid_B=sum(is.finite(v)))}
}
for(st in c('Localized','Regional')) for(su in c('Recorded','None')){
 z<-d[year<=2015 & stage3==st & surgery==su]; Ls<-c(1,3,5)
 Ls<-Ls[vapply(Ls,function(L)all(z[days>L*365.25,.N,by=g]$N>=20) && nrow(z[days>L*365.25,.N,by=g])==2,logical(1))]
 if(length(Ls)) fitjob(z,Ls,list(M0=formula0),paste0('Surgery_',su,'_',st))
}
for(era in c('2010-2015','2016-2018')) for(st in c('Localized','Regional')) {
 y<-as.integer(strsplit(era,'-',fixed=TRUE)[[1]])
 z<-d[year>=y[1]&year<=y[2]&stage3==st]
 ff<-~ns(age,knots=c(50,70),Boundary.knots=c(18,90))+factor(sex)+ns(year,df=2)
 fitjob(z,c(0,1,2),list(M0=ff),paste0('Era_',era,'_',st))
}
rr<-rbindlist(res)
for(met in methods){
 est<-sapply(c('2010-2015','2016-2018'),function(era){x<-rr[method==met];j1<-paste0('Era_',era,'_Localized');j2<-paste0('Era_',era,'_Regional');x[job==j1&landmark==2,difference]-x[job==j1&landmark==0,difference]-x[job==j2&landmark==2,difference]+x[job==j2&landmark==0,difference]})
 bb<-sapply(c('2010-2015','2016-2018'),function(era){a<-bootstore[[paste0('Era_',era,'_Localized')]];b<-bootstore[[paste0('Era_',era,'_Regional')]];a[,1,'2',met,3]-a[,1,'0',met,3]-b[,1,'2',met,3]+b[,1,'0',met,3]})
 summarize_family(c(est,est[2]-est[1]),cbind(bb,bb[,2]-bb[,1]),c('2010-2015','2016-2018','Later minus earlier'),'Era stage change 0 to 2',met)
}
fwrite(rr,file.path(out,'extension_estimates.csv'))
fwrite(rbindlist(diagnostics),file.path(out,'extension_weight_diagnostics.csv'))
fwrite(rbindlist(balances),file.path(out,'extension_balance.csv'))
fwrite(rbindlist(contrasts),file.path(out,'extension_contrasts.csv'))
fwrite(rbindlist(failure),file.path(out,'extension_bootstrap_failures.csv'))
fwrite(rbindlist(warning_log),file.path(out,'extension_fit_warnings.csv'))
writeLines(c('All point-estimate weighted KM calculations independently matched survival::survfit to 1e-8.',capture.output(sessionInfo())),file.path(out,'extension_validation.txt'))
print(rbindlist(contrasts)[method=='IPTW'])
