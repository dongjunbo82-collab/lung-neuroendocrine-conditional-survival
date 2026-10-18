code<-readLines('work/pnen_adjusted.R',warn=FALSE,encoding='UTF-8')
eval(parse(text=code[seq_len(grep('^allres<-',code)[1]-1)]))
out<-'outputs/pnen_formal_sensitivity';dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(20261007);B<-1000L
results<-list();diags<-list()
for(scenario in c('recent_2010_2015','single_primary','histology_confirmed'))for(st in c('Localized','Regional','Distant')) {
 z<-copy(d[stage3==st])
 if(scenario=='recent_2010_2015')z<-z[year>=2010]
 if(scenario=='single_primary')z<-z[total==1]
 if(scenario=='histology_confirmed')z<-z[confirmation==1]
 caches<-lapply(LM,function(L)make_cache(z,L));n<-nrow(z)
 base<-lapply(caches,function(ca)calc(ca,rep(1,n),TRUE))
 boot<-array(NA_real_,c(B,length(LM),3,3));ix0<-which(z$g==0);ix1<-which(z$g==1)
 for(b in seq_len(B)) {
  m<-integer(n);m[ix0]<-tabulate(sample.int(length(ix0),length(ix0),TRUE),length(ix0));m[ix1]<-tabulate(sample.int(length(ix1),length(ix1),TRUE),length(ix1))
  for(j in seq_along(LM))boot[b,j,,]<-calc(caches[[j]],m)
 }
 saveRDS(boot,file.path(out,paste0('bootstrap_',scenario,'_',st,'.rds')))
 for(j in seq_along(LM))for(k in seq_along(methods)) {
  ci<-quantile(boot[,j,k,3],c(.025,.975))
  results[[length(results)+1]]<-data.table(scenario=scenario,stage=st,landmark=LM[j],method=methods[k],LCNEC=base[[j]]$result[k,2],SCLC=base[[j]]$result[k,1],difference=base[[j]]$result[k,3],lower=ci[1],upper=ci[2],B=B)
 }
 for(j in seq_along(LM))diags[[length(diags)+1]]<-copy(base[[j]]$diag)[,`:=`(scenario=scenario,stage=st,landmark=LM[j])]
 cat(scenario,st,'complete\n')
 fwrite(rbindlist(results),file.path(out,'adjusted_sensitivity.csv'))
 fwrite(rbindlist(diags),file.path(out,'sensitivity_diagnostics.csv'))
}
cat('FORMAL_SENSITIVITY_COMPLETE\n')
