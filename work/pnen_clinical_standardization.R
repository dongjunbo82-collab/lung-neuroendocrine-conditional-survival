library(data.table);library(survival);library(splines)
setDTthreads(1);set.seed(20261017);B<-1000L
base<-Sys.getenv('PNEN_PROJECT', getwd())
setwd(base)
out<-'outputs/pnen_clinical';dir.create(out,recursive=TRUE,showWarnings=FALSE)
rawpath<-Sys.getenv('SEER_INPUT')
newpath<-Sys.getenv('SEER_EXTENDED')
src<-readLines('work/pnen_formal_OS.R')
init<-src[seq_len(grep('^allres<-list',src)[1]-1)]
init[grep('^d <- fread',init)]<-"d <- fread(rawpath,sep='\\t')"
eval(parse(text=init))
set.seed(20261017);B<-1000L;out<-'outputs/pnen_clinical'
extra<-fread(newpath,sep='\t',na.strings=NULL,colClasses='character')
stopifnot(ncol(extra)==37L,nrow(extra)==141383L)
setnames(extra,c('id','record',paste0('orig',3:30),'size_raw','cs_size_raw','ajcc6','t6','n6','m6','schema6'))
extra[,`:=`(id=as.numeric(id),record=as.integer(record))]
stopifnot(!anyDuplicated(extra[,.(id,record)]))
d<-merge(d,extra[,c('id','record','size_raw','cs_size_raw','ajcc6','t6','n6','m6','schema6'),with=FALSE],by=c('id','record'),all.x=TRUE,sort=FALSE)
d<-d[stage3!='Unknown'];stopifnot(nrow(d)==67305L,all(d$schema6=='Lung'))
# Numeric millimeter recode values 001--400 (including SEER-assigned range midpoints).
d[,size_mm:=suppressWarnings(as.numeric(fifelse(grepl('^[0-9]{3}$',size_raw),size_raw,NA_character_)))]
d[size_mm<1 | size_mm>400,size_mm:=NA_real_]
d[,size_group:=factor(fcase(is.na(size_mm),'Unknown/special',size_mm<=30,'1-30 mm',size_mm<=50,'31-50 mm',default='51-400 mm'),levels=c('1-30 mm','31-50 mm','51-400 mm','Unknown/special'))]
d[,extent_group:=fcase(stage3=='Localized' & ajcc6=='IA','IA',stage3=='Localized' & ajcc6=='IB','IB',stage3=='Regional' & ajcc6 %in% c('IA','IB','IIA','IIB'),'I/II',stage3=='Regional' & ajcc6=='IIIA','IIIA',stage3=='Regional' & ajcc6=='IIIB','IIIB',default='Other/unknown')]
fwrite(d[,.N,by=.(stage3,g,year,ajcc6,t6,n6,m6)],file.path(out,'stage_year_distribution.csv'))
fwrite(d[,.N,by=.(stage3,g,year,size_raw,size_group)],file.path(out,'size_year_distribution.csv'))
fwrite(d[,.N,by=.(stage3,g,ajcc6)],file.path(out,'summary_stage_AJCC_crosswalk.csv'))
fwrite(d[,.(n=.N,size_unknown_special=sum(is.na(size_mm)),size_unknown_special_percent=100*mean(is.na(size_mm)),AJCC_unknown=sum(ajcc6 %in% c('NA','Blank(s)','','UNK Stage')),AJCC_unknown_percent=100*mean(ajcc6 %in% c('NA','Blank(s)','','UNK Stage'))),by=.(stage3,g,year)],file.path(out,'missingness_by_year.csv'))
support<-rbindlist(lapply(c('Localized','Regional'),function(st)rbindlist(lapply(LM,function(L)d[stage3==st & (days>L*365.25 | (L==0 & days==0)),.N,by=.(g,extent_group,size_group)][,`:=`(stage=st,landmark=L)]))))
fwrite(support,file.path(out,'clinical_cell_support.csv'))
# Reuse the validated KM engine, generalize the design and remove replicate-empty columns.
code<-paste(src[grep('^make_cache<-function',src):(grep('^allres<-list',src)[1]-1)],collapse='\n')
code<-sub('cv<-unique(a[,.(age,sex,year)])','cv<-unique(a[,..covars])',code,fixed=TRUE)
code<-sub('idx<-cv[a,on=.(age,sex,year),which=TRUE]','idx<-cv[a,on=covars,which=TRUE]',code,fixed=TRUE)
old<-"~ns(age,knots=c(50,70),Boundary.knots=c(18,90))+factor(sex)+ns(year,knots=c(2008,2012),Boundary.knots=c(2004,2015))"
code<-sub(old,'model_formula',code,fixed=TRUE)
code<-sub('fit<-glm.fit(ca$X[ok,,drop=FALSE],cbind(cases[ok],count[ok]-cases[ok]),family=binomial())',"X<-ca$X; q<-qr(X[ok,,drop=FALSE]); X<-X[,sort(q$pivot[seq_len(q$rank)]),drop=FALSE]; fit<-glm.fit(X[ok,,drop=FALSE],cbind(cases[ok],count[ok]-cases[ok]),family=binomial(),control=glm.control(maxit=100))",code,fixed=TRUE)
code<-sub('ca$X%*%fit$coefficients','X%*%fit$coefficients',code,fixed=TRUE)
code<-sub('if(any(ps<=0|ps>=1))','if(any((ps<=0|ps>=1) & mult>0))',code,fixed=TRUE)
eval(parse(text=code));methods<-c('crude','IPTW','overlap')
forms<-list(C0=as.formula(old),C1=update(as.formula(old),~.+extent_group),C2=update(as.formula(old),~.+extent_group+size_group))
estimates<-list();diagnostics<-list();balances<-list();warnings<-list();fails<-list();boots<-list();failure_messages<-list();cellstage<-'';cellmodel<-'';cellL<-0;cellB<-0
impl<-calc
calc<-function(...)withCallingHandlers(impl(...),warning=function(w){warnings[[length(warnings)+1L]]<<-data.table(stage=cellstage,model=cellmodel,landmark=cellL,replicate=cellB,message=conditionMessage(w));invokeRestart('muffleWarning')})
for(st in c('Localized','Regional')) {
 cellstage<-st;z<-copy(d[stage3==st]);z[,extent_group:=factor(extent_group)];z[,size_group:=droplevels(size_group)]
 n<-nrow(z);caches<-list();points<-list()
 for(mm in names(forms)) {
  cellmodel<-mm;cellB<-0;model_formula<-forms[[mm]];covars<-all.vars(model_formula)
  caches[[mm]]<-lapply(LM,function(L){cellL<<-L;make_cache(z,L)})
  points[[mm]]<-lapply(seq_along(LM),function(j){cellL<<-LM[j];calc(caches[[mm]][[j]],rep(1,n),TRUE)})
  for(j in seq_along(LM)) {
   diagnostics[[length(diagnostics)+1L]]<-copy(points[[mm]][[j]]$diag)[,`:=`(stage=st,model=mm,landmark=LM[j])]
   balances[[length(balances)+1L]]<-copy(points[[mm]][[j]]$balance)[,`:=`(stage=st,model=mm,landmark=LM[j])]
  }
 }
 bb<-array(NA_real_,c(B,3,5,3,3),dimnames=list(NULL,names(forms),as.character(LM),methods,c('SCLC','LCNEC','difference')))
 i0<-which(z$g==0);i1<-which(z$g==1)
 for(b in seq_len(B)) {
  cellB<-b;m<-integer(n);m[i0]<-tabulate(sample.int(length(i0),length(i0),TRUE),length(i0));m[i1]<-tabulate(sample.int(length(i1),length(i1),TRUE),length(i1))
  for(mm in names(forms))for(j in seq_along(LM)) {
   cellmodel<-mm;cellL<-LM[j];ans<-tryCatch(calc(caches[[mm]][[j]],m),error=function(e)e)
   if(inherits(ans,'error'))failure_messages[[length(failure_messages)+1L]]<-data.table(stage=st,model=mm,landmark=LM[j],replicate=b,error=conditionMessage(ans)) else bb[b,mm,j,,]<-ans
  }
  if(b%%100==0)cat(st,b,'/',B,'\n')
 }
 for(mm in names(forms))for(j in seq_along(LM))for(met in methods) {
  xx<-bb[,mm,j,met,3];valid<-is.finite(xx);ci<-quantile(xx[valid],c(.025,.975))
  rr<-points[[mm]][[j]]$result[met,]
  estimates[[length(estimates)+1L]]<-data.table(stage=st,model=mm,landmark=LM[j],method=met,SCLC=rr[1],LCNEC=rr[2],difference=rr[3],lower=ci[1],upper=ci[2],valid_B=sum(valid))
 }
 boots[[st]]<-bb;saveRDS(bb,file.path(out,paste0('bootstrap_',st,'.rds')))
}
rr<-rbindlist(estimates);cc<-list()
for(met in methods){
 ee<-sapply(names(forms),function(mm){x<-rr[model==mm & method==met];x[stage=='Localized' & landmark==5,difference]-x[stage=='Localized' & landmark==0,difference]-x[stage=='Regional' & landmark==5,difference]+x[stage=='Regional' & landmark==0,difference]})
 bb<-sapply(names(forms),function(mm)boots$Localized[,mm,'5',met,3]-boots$Localized[,mm,'0',met,3]-boots$Regional[,mm,'5',met,3]+boots$Regional[,mm,'0',met,3])
 ok<-apply(bb,1,function(x)all(is.finite(x)));bf<-bb[ok,,drop=FALSE];se<-apply(bf,2,sd);crit<-quantile(apply(abs(sweep(sweep(bf,2,ee,'-'),2,se,'/')),1,max),.95)
 for(k in 1:3)cc[[length(cc)+1L]]<-data.table(method=met,model=names(forms)[k],contrast='Localized minus Regional change 0 to 5',estimate=ee[k],lower=quantile(bf[,k],.025),upper=quantile(bf[,k],.975),sim_lower=ee[k]-crit*se[k],sim_upper=ee[k]+crit*se[k],valid_B=sum(ok))
 for(k in 2:3){v<-bf[,k]-bf[,1];cc[[length(cc)+1L]]<-data.table(method=met,model=paste0(names(forms)[k],'-C0'),contrast='Paired change in stage-change contrast',estimate=ee[k]-ee[1],lower=quantile(v,.025),upper=quantile(v,.975),sim_lower=NA_real_,sim_upper=NA_real_,valid_B=sum(ok))}
}
fwrite(rr,file.path(out,'clinical_estimates.csv'));fwrite(rbindlist(cc),file.path(out,'clinical_change_contrasts.csv'))
fwrite(rbindlist(diagnostics),file.path(out,'clinical_weight_diagnostics.csv'));fwrite(rbindlist(balances),file.path(out,'clinical_balance.csv'))
fwrite(rbindlist(warnings,fill=TRUE),file.path(out,'clinical_fit_warnings.csv'))
if(length(failure_messages))fwrite(rbindlist(failure_messages),file.path(out,'clinical_failures.csv')) else writeLines('No failed bootstrap cells.',file.path(out,'clinical_failures.txt'))
primary<-fread('outputs/pnen_formal/adjusted_OS_differences.csv')[stage %in% c('Localized','Regional')]
check<-merge(rr[model=='C0'],primary,by=c('stage','landmark','method'),suffixes=c('_new','_original'))
stopifnot(max(abs(check$difference_new-check$difference_original))<1e-10)
writeLines(c('Exploratory clinical standardization. Seed 20261017; 1000 joint bootstrap replicates.', 'Original main cohort preserved; full cohort retained using unknown/special categories. No outcome-driven trimming.', 'C0 demographic; C1 adds stage-specific AJCC6 grouping; C2 adds tumor-size category.', 'Localized extent: IA, IB, other/unknown. Regional extent: I/II, IIIA, IIIB, other/unknown.', 'Size groups: numeric recode 001-030, 031-050, 051-400 mm; all special/unknown codes a separate category.', 'Common patient resamples across landmarks and models; model columns absent from a replicate are removed via QR.', 'Separate supporting simultaneous family across C0/C1/C2 focused stage-change contrasts within each weighting method; paired model differences pointwise.', 'All point-estimate KM results independently verified against survival::survfit; C0 point estimates reproduce original results within 1e-10.',capture.output(sessionInfo())),file.path(out,'clinical_validation_and_session.txt'))
print(rbindlist(cc)[method=='IPTW']);print(rbindlist(diagnostics)[method=='IPTW',.(max_weight=max(max_normalized_weight),min_ess=min(ess)),by=model])
