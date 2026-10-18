library(data.table)
library(survival)
library(splines)
setDTthreads(1)
set.seed(20261005)
B <- 2000L
out <- 'outputs/pnen_formal_competing'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- fread(Sys.getenv('SEER_INPUT'),sep='\t')
setnames(d,names(d),c('id','record','sequence','total','sex','race','year','marital','site','laterality','hist','behavior','confirmation','source','grade_old','grade_new','stage','surgery_old','surgery_new','no_surgery','radiation','chemo','months','days','flag','vital','cancer','other','cod','age'))
d[,`:=`(year=year+1800L,hist=hist+8000L)]
setorder(d,id,year,sequence,record)
d<-d[!duplicated(id)]
d<-d[days!=32765 & !is.na(days) & !source %in% c(6,7) & flag!=8 & vital %in% c(0,1) & year<=2015 & hist %in% c(8013,8041)]
stopifnot(all(d$sex %in% c(1,2)),all(d$age %in% 18:90))
d[,`:=`(g=as.integer(hist==8013),event=as.integer(vital==0),stage3=fcase(stage==1,'Localized',stage %in% 2:4,'Regional',stage==7,'Distant',default='Unknown'))]
LM<-c(0,1,2,3,5)
methods<-c('crude','IPTW','overlap')
make_cache<-function(z,L) {
 ix<-which(z$days>L*365.25 | (L==0 & z$days==0)); a<-z[ix]
 cv<-unique(a[,.(age,sex,year)])
 idx<-cv[a,on=.(age,sex,year),which=TRUE]
 X<-model.matrix(~ns(age,knots=c(50,70),Boundary.knots=c(18,90))+factor(sex)+ns(year,knots=c(2008,2012),Boundary.knots=c(2004,2015)),cv)
 qrX<-qr(X)
 X<-X[,sort(qrX$pivot[seq_len(qrX$rank)]),drop=FALSE]
 grp<-lapply(0:1,function(g) {
  j<-which(a$g==g); j<-j[order(a$days[j])]
  ends<-c(which(diff(a$days[j])!=0),length(j))
  list(j=j,ends=ends,times=a$days[j][ends])
 })
 list(ix=ix,a=a,X=X,idx=idx,grp=grp,L=L)
}
calc<-function(ca,m,diagnostics=FALSE) {
 a<-ca$a; mult<-m[ca$ix]
 count<-as.numeric(rowsum(mult,ca$idx,reorder=FALSE)[match(seq_len(nrow(ca$X)),unique(ca$idx)),1])
 cases<-as.numeric(rowsum(mult*a$g,ca$idx,reorder=FALSE)[match(seq_len(nrow(ca$X)),unique(ca$idx)),1])
 ok<-count>0
 fit<-glm.fit(ca$X[ok,,drop=FALSE],cbind(cases[ok],count[ok]-cases[ok]),family=binomial())
 if(!fit$converged || any(!is.finite(fit$coefficients))) stop('PS fit failed')
 ps<-plogis(drop(ca$X%*%fit$coefficients))[ca$idx]
 if(any(ps<=0|ps>=1)) stop('Extreme fitted propensity')
 ws<-list(crude=rep(1,nrow(a)),IPTW=ifelse(a$g==1,1/ps,1/(1-ps)),overlap=ifelse(a$g==1,1-ps,ps))
 result<-matrix(NA_real_,3,3,dimnames=list(methods,c('SCLC','LCNEC','difference')))
 diag<-list()
 for(method in methods) {
  w<-ws[[method]]*mult
  for(g in 0:1) {
   gr<-ca$grp[[g+1]]; j<-gr$j; ends<-gr$ends
   cum<-cumsum(w[j]); previous<-c(0,head(cum[ends],-1))
   ev<-cumsum(w[j]*a$event[j]); deaths<-diff(c(0,ev[ends]))
   risk<-sum(w[j])-previous
   valid<-gr$times<=(ca$L+3)*365.25 & risk>1e-10
   if(sum(w[j][a$days[j]>=(ca$L+3)*365.25])==0) stop('No endpoint support')
   result[method,g+1]<-prod(1-pmin(1,deaths[valid]/risk[valid]))
   if(diagnostics) {
    ww<-w[j]; norm<-ww/mean(ww)
    diag[[length(diag)+1L]]<-data.table(method=method,group=c('SCLC','LCNEC')[g+1],n=length(j),ess=sum(ww)^2/sum(ww^2),max_normalized_weight=max(norm),ps_min=min(ps[j]),ps_max=max(ps[j]),early_censor_n=sum(a$event[j]==0 & a$days[j]<(ca$L+3)*365.25))
   }
  }
  result[method,3]<-result[method,2]-result[method,1]
 }
 if(!diagnostics) return(result)
 features<-cbind(age=a$age,year=a$year,female=as.numeric(a$sex==2),ca$X[ca$idx,-1,drop=FALSE])
 balance<-rbindlist(lapply(methods,function(method) {
  w<-ws[[method]]
  smds<-apply(features,2,function(v) {
   den<-sqrt((var(v[a$g==0])+var(v[a$g==1]))/2)
   if(den==0) return(0)
   (weighted.mean(v[a$g==1],w[a$g==1])-weighted.mean(v[a$g==0],w[a$g==0]))/den
  })
  data.table(method=method,feature=names(smds),smd=smds)
 }))
 # Independently check manual weighted KM against survival::survfit.
 for(method in methods) for(g in 0:1) {
  selected<-which(a$g==g)
  aa<-copy(a[selected]); ww<-ws[[method]][selected]
  sf<-survfit(Surv(days/365.25-ca$L,event)~1,aa,weights=ww)
  ref<-summary(sf,times=3,extend=FALSE)$surv[1]
  stopifnot(abs(ref-result[method,g+1])<1e-8)
 }
 list(result=result,diag=rbindlist(diag),balance=balance)
}
calcAJ<-function(ca,m,diagnostics=FALSE) {
 a<-ca$a; mult<-m[ca$ix]
 count<-as.numeric(rowsum(mult,ca$idx,reorder=FALSE)[match(seq_len(nrow(ca$X)),unique(ca$idx)),1])
 cases<-as.numeric(rowsum(mult*a$g,ca$idx,reorder=FALSE)[match(seq_len(nrow(ca$X)),unique(ca$idx)),1])
 ok<-count>0
 fit<-glm.fit(ca$X[ok,,drop=FALSE],cbind(cases[ok],count[ok]-cases[ok]),family=binomial())
 if(!fit$converged || any(!is.finite(fit$coefficients))) stop('PS fit failed')
 ps<-plogis(drop(ca$X%*%fit$coefficients))[ca$idx]
 if(any(ps<=0|ps>=1)) stop('Extreme fitted propensity')
 ws<-list(crude=rep(1,nrow(a)),IPTW=ifelse(a$g==1,1/ps,1/(1-ps)),overlap=ifelse(a$g==1,1-ps,ps))

 result<-array(NA_real_,c(3,2,4),dimnames=list(methods,c('SCLC','LCNEC'),c('OS','cancer','other','unknown')))
 cause<-fcase(a$event==0,0L,a$cancer==1,1L,a$other==1,2L,default=3L)
 for(method in methods) for(g in 0:1) {
  gr<-ca$grp[[g+1]];j<-gr$j; ends<-gr$ends; w<-ws[[method]]*mult
  cum<-cumsum(w[j]); risk<-sum(w[j])-c(0,head(cum[ends],-1))
  deaths<-sapply(1:3,function(k)diff(c(0,cumsum(w[j]*(cause[j]==k))[ends])))
  valid<-gr$times<=(ca$L+3)*365.25 & risk>1e-10
  hz<-deaths[valid,,drop=FALSE]/risk[valid]
  S<-cumprod(1-pmin(1,rowSums(hz)))
  cif<-colSums(hz*c(1,head(S,-1)))
  result[method,g+1,]<-c(tail(S,1),cif)
  stopifnot(abs(sum(result[method,g+1,])-1)<1e-8)
  if(diagnostics) {
   aa<-copy(a[j]); aa[,cause:=factor(cause[j],levels=0:3,labels=c('censor','cancer','other','unknown'))]
   sf<-survfit(Surv(days/365.25-ca$L,cause)~1,aa,weights=w[j])
   ref<-summary(sf,times=3,extend=FALSE)$pstate[1,]
   stopifnot(max(abs(ref-result[method,g+1,]))<1e-8)
  }
 }
 result
}
results<-list()
for(st in c('Localized','Regional','Distant')) {
 z<-d[stage3==st];n<-nrow(z);caches<-lapply(LM,function(L)make_cache(z,L))
 base<-lapply(caches,function(ca)calcAJ(ca,rep(1,n),TRUE))
 boots<-array(NA_real_,c(B,length(LM),3,2,4))
 idx0<-which(z$g==0);idx1<-which(z$g==1)
 for(b in seq_len(B)) {
  m<-integer(n)
  m[idx0]<-tabulate(sample.int(length(idx0),length(idx0),replace=TRUE),length(idx0))
  m[idx1]<-tabulate(sample.int(length(idx1),length(idx1),replace=TRUE),length(idx1))
  for(j in seq_along(LM)) boots[b,j,,,]<-calcAJ(caches[[j]],m)
  if(b%%100==0)cat(st,b,'/',B,'\n')
 }
 saveRDS(boots,file.path(out,paste0('bootstrap_',st,'.rds')))
 for(j in seq_along(LM))for(k in seq_along(methods))for(o in 1:4) {
  v<-base[[j]][k,,o];dif<-boots[,j,k,2,o]-boots[,j,k,1,o]
  results[[length(results)+1]]<-data.table(stage=st,landmark=LM[j],method=methods[k],outcome=c('OS','cancer_CIF','other_CIF','unknown_CIF')[o],SCLC=v[1],LCNEC=v[2],difference=v[2]-v[1],lower=quantile(dif,.025),upper=quantile(dif,.975))
 }
}
fwrite(rbindlist(results),file.path(out,'conditional_competing_differences.csv'))
writeLines(c('All point estimates verified against survival::survfit multistate weighted Aalen-Johansen, tolerance1e-8.','Unknown death cause retained as separate competing event. Bootstrap2000 joint across landmarks, seed20261005.',capture.output(sessionInfo())),file.path(out,'validation_and_session.txt'))
cat('FORMAL_AJ_COMPLETE\n')
