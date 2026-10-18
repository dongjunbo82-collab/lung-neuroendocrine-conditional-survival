library(data.table)
f<-Sys.getenv('SEER_INPUT')
d<-fread(f,sep='\t');setnames(d,c('id','record','sequence','total','sex','race','year','marital','site','laterality','hist','behavior','confirmation','source','grade_old','grade_new','stage','surgery_old','surgery_new','no_surgery','radiation','chemo','months','days','flag','vital','cancer','other','cod','age'))
d[,`:=`(year=year+1800L,hist=hist+8000L)];setorder(d,id,year,sequence,record);d<-d[!duplicated(id)];d<-d[days!=32765 & !is.na(days) & !source %in% c(6,7) & flag!=8 & vital %in% c(0,1) & year<=2015 & hist %in% c(8013,8041)]
d[,stage3:=fcase(stage==1,'Localized',stage %in% 2:4,'Regional',default='Other')]
a<-d[stage3!='Other' & days>5*365.25,.(n=.N,deaths1=sum(vital==0 & days<=6*365.25),deaths3=sum(vital==0 & days<=8*365.25)),by=.(stage3,hist)]
print(a)
