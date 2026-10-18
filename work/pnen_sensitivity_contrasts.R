library(data.table)
stages<-c('Localized','Regional','Distant');methods<-c('crude','IPTW','overlap')
result<-list()
for(sc in c('main','recent_2010_2015','single_primary','histology_confirmed')) {
 folder<-if(sc=='main')'outputs/pnen_formal' else 'outputs/pnen_formal_sensitivity'
 estimates<-fread(file.path(folder,if(sc=='main')'adjusted_OS_differences.csv' else 'adjusted_sensitivity.csv'))
 if(sc!='main')estimates<-estimates[scenario==sc]
 boots<-lapply(stages,function(st)readRDS(file.path(folder,if(sc=='main')paste0('bootstrap_',st,'.rds') else paste0('bootstrap_',sc,'_',st,'.rds'))))
 for(k in 1:3) {
  pair<-combn(1:3,2); B<-dim(boots[[1]])[1];draws<-matrix(NA_real_,B,3);est<-numeric(3)
  change<-function(st) {
   end<-estimates[stage==st & method==methods[k] & landmark==5,difference]
   start<-estimates[stage==st & method==methods[k] & landmark==0,difference]
   stopifnot(length(end)==1L,length(start)==1L,is.finite(end),is.finite(start))
   end-start
  }
  for(i in 1:3) {
   a<-pair[1,i];b<-pair[2,i]
   draws[,i]<-(boots[[a]][,5,k,3]-boots[[a]][,1,k,3])-(boots[[b]][,5,k,3]-boots[[b]][,1,k,3])
   est[i]<-change(stages[a])-change(stages[b])
  }
  se<-apply(draws,2,sd);crit<-quantile(apply(abs(sweep(sweep(draws,2,est,'-'),2,se,'/')),1,max),.95)
  for(i in 1:3)result[[length(result)+1]]<-data.table(scenario=sc,method=methods[k],stage_A=stages[pair[1,i]],stage_B=stages[pair[2,i]],estimate=est[i],lower=quantile(draws[,i],.025),upper=quantile(draws[,i],.975),sim_lower=est[i]-crit*se[i],sim_upper=est[i]+crit*se[i],B=B)
 }
}
r<-rbindlist(result)
numeric_fields<-c('estimate','lower','upper','sim_lower','sim_upper')
stopifnot(nrow(r)==36,all(is.finite(as.matrix(r[,..numeric_fields]))),all(r$lower<=r$upper),all(r$sim_lower<=r$sim_upper))
previous<-fread('outputs/pnen_formal/cross_stage_change_contrasts.csv')
current<-r[scenario=='main' & method=='IPTW']
stopifnot(identical(current$stage_A,previous$stage_A),identical(current$stage_B,previous$stage_B),max(abs(as.matrix(current[,..numeric_fields])-as.matrix(previous[,..numeric_fields])))<1e-10)
fwrite(r,'outputs/pnen_formal_sensitivity/change_contrasts.csv')
print(r[method=='IPTW' & stage_A=='Localized' & stage_B=='Regional'])

