library(data.table)
library(survival)
setDTthreads(1)
out <- 'outputs/pnen_extension'; dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- fread(Sys.getenv('SEER_INPUT'),sep='\t')
setnames(d,names(d),c('id','record','sequence','total','sex','race','year','marital','site','laterality','hist','behavior','confirmation','source','grade_old','grade_new','stage','surgery_old','surgery_new','no_surgery','radiation','chemo','months','days','flag','vital','cancer','other','cod','age'))
d[,`:=`(year=year+1800L,hist=hist+8000L)]
flow <- data.table(step='Exported records',n=nrow(d))
setorder(d,id,year,sequence,record); d<-d[!duplicated(id)]
flow<-rbind(flow,data.table(step='Unique patients with earliest eligible record',n=nrow(d)))
d<-d[days!=32765 & !is.na(days) & !source %in% c(6,7) & flag!=8 & vital %in% c(0,1)]
flow<-rbind(flow,data.table(step='Valid survival across exported histologies and years',n=nrow(d)))
d<-d[hist %in% c(8013,8041)]
d[,`:=`(g=as.integer(hist==8013),event=as.integer(vital==0),histology=ifelse(hist==8013,'LCNEC','SCLC'),stage3=fcase(stage==1,'Localized',stage %in% 2:4,'Regional',stage==7,'Distant',default='Unknown'))]
# All inferential extension cohorts end by 2018; 2023 surgery field is not used.
d[,surgery:=fcase(surgery_old==0,'None',surgery_old %in% 10:90,'Recorded',default='Unknown')]
d[,radiotherapy:=fcase(radiation %in% 1:5,'Recorded',radiation==7,'Refused',radiation==8,'Recommended unknown',default='None unknown')]
d[,rt_yes:=as.integer(radiation %in% 1:5)]
fwrite(d[year<=2022,.N,by=.(year,surgery_old,radiation,chemo)],file.path(out,'treatment_code_audit.csv'))
h<-d[year<=2015]
flow<-rbind(flow,data.table(step=c('Historical LCNEC or SCLC','Unknown stage excluded','Historical known stage'),n=c(nrow(h),nrow(h[stage3=='Unknown']),nrow(h[stage3!='Unknown']))))
fwrite(flow,file.path(out,'cohort_flow.csv'))
fwrite(h[,.(n=.N,age_median=median(age),age_q1=quantile(age,.25),age_q3=quantile(age,.75),female_n=sum(sex==2),LCNEC_n=sum(g),surgery_n=sum(surgery=='Recorded'),surgery_unknown_n=sum(surgery=='Unknown'),radiotherapy_n=sum(rt_yes),chemotherapy_n=sum(chemo==1),year_median=median(year),hist_confirmed_n=sum(confirmation==1),death_n=sum(event)),by=.(stage_known=stage3!='Unknown')],file.path(out,'unknown_stage_comparison.csv'))
features<-c('age','year','sex','g','surgery','rt_yes','chemo','confirmation')
bal<-rbindlist(lapply(features,function(v){
 x<-h[[v]]; if(v %in% c('sex','confirmation')) x<-as.integer(x==ifelse(v=='sex',2,1))
 if(v=='surgery') x<-as.integer(x=='Recorded')
 a<-x[h$stage3=='Unknown']; b<-x[h$stage3!='Unknown']
 data.table(variable=v,unknown_mean=mean(a),known_mean=mean(b),SMD=(mean(a)-mean(b))/sqrt((var(a)+var(b))/2))
}))
fwrite(bal,file.path(out,'unknown_stage_SMD.csv'))
risk<-rbindlist(lapply(c(0,1,2,3,5),function(L)h[stage3!='Unknown' & (days>L*365.25 | (L==0 & days==0)),.(n=.N,deaths_next3=sum(event==1 & days<=(L+3)*365.25),endpoint_support=sum(days>=(L+3)*365.25)),by=.(stage3,histology,surgery)][,landmark:=L]))
fwrite(risk,file.path(out,'surgery_risk_sets.csv'))
era<-rbindlist(lapply(list(c(2010,2015),c(2016,2018)),function(y)rbindlist(lapply(c(0,1,2),function(L)d[year>=y[1] & year<=y[2] & stage3!='Unknown' & (days>L*365.25 | (L==0 & days==0)),.(n=.N,deaths_next3=sum(event==1 & days<=(L+3)*365.25),endpoint_support=sum(days>=(L+3)*365.25),early_censor=sum(event==0 & days<(L+3)*365.25)),by=.(stage3,histology)][,`:=`(landmark=L,era=paste(y,collapse='-'))]))))
fwrite(era,file.path(out,'era_support.csv'))
saveRDS(d,'work/pnen_extension_data.rds')
print(risk[stage3 %in% c('Localized','Regional') & landmark %in% c(0,5)])
print(era)
print(bal)
