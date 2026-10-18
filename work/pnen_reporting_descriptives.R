code<-readLines('work/pnen_formal_OS.R',warn=FALSE,encoding='UTF-8')
eval(parse(text=code[seq_len(grep('^LM<-',code)[1]-1)]))
z<-d[stage3!='Unknown']
stopifnot(nrow(z)==67305)
summarize<-function(a) {
 sf<-survfit(Surv(days/365.25,1-event)~1,data=a)
 list(n=nrow(a),deaths=sum(a$event),person_years=sum(a$days)/365.25,
      reverse_KM_median_years=as.numeric(summary(sf)$table['median']),
      missing_age=sum(is.na(a$age)),missing_sex=sum(is.na(a$sex)),missing_year=sum(is.na(a$year)))
}
r<-rbindlist(list(as.data.table(summarize(z))[,hist:='All'],z[,summarize(.SD),by=hist]),use.names=TRUE)
fwrite(r,'outputs/pnen_formal/reporting_descriptives.csv')
print(r)
