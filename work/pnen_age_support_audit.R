library(data.table)
d<-readRDS('work/pnen_extension_data.rds')[year<=2015 & stage3 %in% c('Localized','Regional')]
d[,age_group:=fcase(age<65,'<65',age<75,'65-74',default='75+')]
r<-rbindlist(lapply(c(0,3,5),function(L)d[days>L*365.25 | (L==0 & days==0),.(n=.N,deaths_next3=sum(event==1 & days<=(L+3)*365.25),other_deaths_next3=sum(event==1 & cancer!=1 & other==1 & days<=(L+3)*365.25),endpoint_support=sum(days>=(L+3)*365.25)),by=.(stage3,histology,age_group)][,landmark:=L]))
fwrite(r,'outputs/pnen_extension/age_support_audit.csv')
print(r[histology=='LCNEC' & landmark==5])
