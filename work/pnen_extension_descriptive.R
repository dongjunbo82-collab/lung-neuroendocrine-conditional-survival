library(data.table)
d<-readRDS('work/pnen_extension_data.rds');out<-'outputs/pnen_extension'
h<-d[year<=2015 & stage3!='Unknown']
fwrite(h[,.(n=.N,surgery_recorded=sum(surgery=='Recorded'),surgery_none=sum(surgery=='None'),surgery_unknown=sum(surgery=='Unknown'),radiation_recorded=sum(rt_yes),radiation_none_unknown=sum(radiation==0),radiation_refused=sum(radiation==7),radiation_recommended_unknown=sum(radiation==8),chemo_recorded=sum(chemo==1)),by=histology],file.path(out,'baseline_treatment.csv'))
fwrite(h[,.(n=.N,surgery_recorded=sum(surgery=='Recorded'),surgery_percent=100*mean(surgery=='Recorded')),by=.(stage3,histology)],file.path(out,'stage_surgery_percent.csv'))
fwrite(h[stage3 %in% c('Localized','Regional'),.(n=.N,excluded_unknown_surgery=sum(surgery=='Unknown')),by=.(stage3,histology)],file.path(out,'treatment_exclusions.csv'))
fwrite(d[year %in% 2010:2018,.N,by=.(era=ifelse(year<=2015,'2010-2015','2016-2018'),stage3,histology)],file.path(out,'era_cohort_counts.csv'))
fwrite(h[,.(n=.N),by=.(histology,surgery_old)],file.path(out,'surgery_codes_historical.csv'))
