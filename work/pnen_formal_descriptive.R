code<-readLines('work/pnen_adjusted.R',warn=FALSE,encoding='UTF-8')
eval(parse(text=code[seq_len(grep('^allres<-',code)[1]-1)]))
out<-'outputs/pnen_formal'
flow<-data.table(step=c('Eligible LCNEC/SCLC, diagnosis2004-2015, known/unknown stage','Known initial stage'),n=c(nrow(d),nrow(d[stage3!='Unknown'])))
fwrite(flow,file.path(out,'primary_cohort_flow.csv'))
fwrite(d[,.N,by=.(hist,stage3)][order(hist,stage3)],file.path(out,'primary_stage_counts.csv'))
z<-d[stage3!='Unknown']
fwrite(z[,.(n=.N,age_median=as.numeric(median(age)),age_q1=quantile(age,.25),age_q3=quantile(age,.75),female=sum(sex==2),year_median=as.numeric(median(year)),single_primary=sum(total==1),histologically_confirmed=sum(confirmation==1),chemo_yes=sum(chemo==1),deaths=sum(event==1),unknown_cause_deaths=sum(event==1 & cancer!=1 & other!=1)),by=.(hist,stage3)],file.path(out,'baseline_by_stage_histology.csv'))
fwrite(z[,.(n=.N,age_median=as.numeric(median(age)),age_q1=quantile(age,.25),age_q3=quantile(age,.75),female=sum(sex==2),year_median=as.numeric(median(year)),single_primary=sum(total==1),histologically_confirmed=sum(confirmation==1),chemo_yes=sum(chemo==1),deaths=sum(event==1)),by=hist],file.path(out,'baseline_by_histology.csv'))
for(v in c('race','confirmation','source','stage','cancer','other','chemo','surgery_old','radiation'))fwrite(z[,.N,by=c('hist',v)],file.path(out,paste0('audit_',v,'.csv')))
fwrite(rbindlist(lapply(LM,function(L)z[days>L*365.25 | (L==0 & days==0),.(n=.N,age_median=as.numeric(median(age)),female_pct=100*mean(sex==2),year_median=as.numeric(median(year)),risk_end=sum(days>=(L+3)*365.25),events_horizon=sum(event==1 & days<=(L+3)*365.25)),by=.(hist,stage3)][,landmark:=L])),file.path(out,'risk_sets.csv'))

