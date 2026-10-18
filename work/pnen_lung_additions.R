src<-readLines('work/pnen_formal_competing.R');eval(parse(text=src[1:(grep('^results<-list',src)[1]-1)]))
out<-'outputs/pnen_lung_additions';dir.create(out,recursive=TRUE,showWarnings=FALSE)
B<-1000L;set.seed(20261014)
rmst<-function(ca,m,check=FALSE){
 a<-ca$a;mult<-m[ca$ix]
 count<-as.numeric(rowsum(mult,ca$idx,reorder=FALSE)[match(seq_len(nrow(ca$X)),unique(ca$idx)),1]);cases<-as.numeric(rowsum(mult*a$g,ca$idx,reorder=FALSE)[match(seq_len(nrow(ca$X)),unique(ca$idx)),1]);ok<-count>0
 fit<-glm.fit(ca$X[ok,,drop=FALSE],cbind(cases[ok],count[ok]-cases[ok]),family=binomial());stopifnot(fit$converged,all(is.finite(fit$coefficients)))
 ps<-plogis(drop(ca$X%*%fit$coefficients))[ca$idx];stopifnot(all(ps>0 & ps<1));w<-ifelse(a$g==1,1/ps,1/(1-ps))*mult
 res<-numeric(2)
 for(g in 0:1){
  gr<-ca$grp[[g+1]];j<-gr$j;ends<-gr$ends;cum<-cumsum(w[j]);risk<-sum(w[j])-c(0,head(cum[ends],-1));death<-diff(c(0,cumsum(w[j]*a$event[j])[ends]));valid<-gr$times<=(ca$L+3)*365.25 & risk>1e-10
  tt<-gr$times[valid]/365.25-ca$L;ss<-cumprod(1-pmin(1,death[valid]/risk[valid]));res[g+1]<-sum(diff(c(0,tt,3))*c(1,ss))*12
  if(check){sf<-survfit(Surv(days/365.25-ca$L,event)~1,a[j],weights=w[j]);keep<-sf$time<=3;ref<-sum(diff(c(0,sf$time[keep],3))*c(1,sf$surv[keep]))*12;stopifnot(abs(ref-res[g+1])<1e-8)}
 };c(res,res[2]-res[1])
}
ar<-list();diag<-list();bal<-list();rr<-list();warn<-character()
withCallingHandlers({
 for(st in c('Localized','Regional')){
  for(ag in c('<70','70+')){
   z<-d[stage3==st & if(ag=='<70')age<70 else age>=70];n<-nrow(z);cs<-lapply(LM,function(L)make_cache(z,L));base<-lapply(cs,function(ca)calcAJ(ca,rep(1,n),TRUE));bo<-array(NA_real_,c(B,5,2,4));ids<-lapply(0:1,function(g)which(z$g==g))
   for(j in 1:5){v<-calc(cs[[j]],rep(1,n),TRUE);diag[[length(diag)+1]]<-copy(v$diag)[method=='IPTW'][,`:=`(stage=st,age_group=ag,landmark=LM[j])];bal[[length(bal)+1]]<-copy(v$balance)[method=='IPTW'][,`:=`(stage=st,age_group=ag,landmark=LM[j])];a<-cs[[j]]$a
    for(g in 0:1)diag[[length(diag)]] [group==c('SCLC','LCNEC')[g+1],`:=`(deaths_next3=sum(a$event[a$g==g]==1 & a$days[a$g==g]<=(LM[j]+3)*365.25),endpoint_support=sum(a$g==g & a$days>=(LM[j]+3)*365.25))]
   }
   for(b in 1:B){m<-integer(n);for(ii in ids)m[ii]<-tabulate(sample.int(length(ii),length(ii),replace=TRUE),length(ii));for(j in 1:5)bo[b,j,,]<-calcAJ(cs[[j]],m)[2,,];if(b%%250==0)cat('AGE',st,ag,b,'\n')}
   for(j in 1:5)for(o in 1:4){v<-base[[j]][2,,o];dd<-bo[,j,2,o]-bo[,j,1,o];ar[[length(ar)+1]]<-data.table(stage=st,age_group=ag,landmark=LM[j],outcome=c('OS','cancer_CIF','other_CIF','unknown_CIF')[o],SCLC=v[1],LCNEC=v[2],difference=v[2]-v[1],lower=quantile(dd,.025),upper=quantile(dd,.975),B=B)}
   fwrite(rbindlist(ar),file.path(out,'age_estimates.csv'));fwrite(rbindlist(diag),file.path(out,'age_support.csv'));fwrite(rbindlist(bal),file.path(out,'age_balance.csv'))
  }
  z<-d[stage3==st];n<-nrow(z);cs<-lapply(LM,function(L)make_cache(z,L));bs<-lapply(cs,function(ca)rmst(ca,rep(1,n),TRUE));bo<-array(NA_real_,c(B,5,3));ids<-lapply(0:1,function(g)which(z$g==g))
  for(b in 1:B){m<-integer(n);for(ii in ids)m[ii]<-tabulate(sample.int(length(ii),length(ii),replace=TRUE),length(ii));for(j in 1:5)bo[b,j,]<-rmst(cs[[j]],m);if(b%%250==0)cat('RMST',st,b,'\n')}
  for(j in 1:5)rr[[length(rr)+1]]<-data.table(stage=st,landmark=LM[j],horizon_years=3,unit='months',SCLC=bs[[j]][1],LCNEC=bs[[j]][2],difference=bs[[j]][3],lower=quantile(bo[,j,3],.025),upper=quantile(bo[,j,3],.975),B=B)
  fwrite(rbindlist(rr),file.path(out,'RMST_estimates.csv'))
 }
},warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')})
fwrite(data.table(message=warn),file.path(out,'fit_warnings.csv'))
writeLines(c('Age at diagnosis <70 versus >=70; same historical cohort and landmark definitions.','B1000 seed20261014; joint patient bootstrap stratified by stage, age group and histology for age estimates; by stage and histology for RMST.','All pointwise secondary intervals; no multiplicity-adjusted age interaction or RMST change tests.','Weighted AJ and RMST independent survival::survfit checks passed tolerance1e-8. No failed resamples omitted.',capture.output(sessionInfo())),file.path(out,'validation.txt'))
cat('COMPLETE\n')
