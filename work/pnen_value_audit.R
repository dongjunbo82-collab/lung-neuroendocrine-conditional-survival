library(data.table)
library(survival)
library(splines)
setDTthreads(1)
base <- normalizePath(getwd(), winslash='/')
out <- file.path(base,'outputs/pnen_value_audit')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src <- readLines(file.path(base,'work/pnen_formal_OS.R'))
eval(parse(text=src[seq_len(which(grepl('^allres<-',src))[1]-1)]))
out <- file.path(base,'outputs/pnen_value_audit')
d[,cause:=fcase(event==0,0L,cancer==1,1L,other==1,2L,default=3L)]
stopifnot(nrow(d[stage3!='Unknown'])==67305L)
fwrite(d[stage3!='Unknown',.N,by=.(event,cancer,other,cod)],file.path(out,'cause_classification_audit.csv'))
counts<-list();support<-list();probs<-list()
for(st in c('Localized','Regional')) {
 z<-d[stage3==st]
 ref<-z[days>5*365.25,.(ref_n=.N),by=.(age,sex,year)]
 for(L in LM) {
  a<-z[days>L*365.25 | (L==0 & days==0)]
  for(gval in 0:1) {
   zz<-a[g==gval]
   cc<-zz[,.(source_n=.N),by=.(age,sex,year)]
   match<-merge(ref,cc,by=c('age','sex','year'),all.x=TRUE)
   support[[length(support)+1]]<-data.table(stage=st,landmark=L,histology=c('SCLC','LCNEC')[gval+1],reference='Pooled 5-year survivors within stage',reference_n=sum(ref$ref_n),reference_cells=nrow(ref),missing_cells=sum(is.na(match$source_n)),reference_mass_unrepresented=sum(match$ref_n[is.na(match$source_n)])/sum(ref$ref_n),min_nonzero_cell=min(cc$source_n))
  }
  b<-readRDS(file.path(base,'outputs/pnen_formal',paste0('bootstrap_',st,'.rds')))
  rr<-fread(file.path(base,'outputs/pnen_formal/adjusted_OS_differences.csv'))[stage==st & landmark==L & method=='IPTW']
  j<-match(L,LM)
  probs[[length(probs)+1]]<-data.table(stage=st,landmark=L,LCNEC=rr$LCNEC,LCNEC_lower=quantile(b[,j,2,2],.025),LCNEC_upper=quantile(b[,j,2,2],.975),SCLC=rr$SCLC,SCLC_lower=quantile(b[,j,2,1],.025),SCLC_upper=quantile(b[,j,2,1],.975),difference=rr$difference,lower=rr$lower,upper=rr$upper,n_LCNEC=sum(a$g==1),n_SCLC=sum(a$g==0),B=dim(b)[1])
  if(L==5) {
   a[,cod_group:=fcase(cause==1,'Cancer-related',cause==3,'Unknown',cause==0,'Censored/alive',cod %in% 106:111,'Circulatory',cod %in% 112:113,'Respiratory',cod %in% 1:79,'Other neoplasm',default='Other specified/residual')]
   counts[[length(counts)+1]]<-a[event==1 & days<=(L+3)*365.25,.N,by=.(g,cause,cod_group)][,stage:=st]
  }
 }
}
fwrite(rbindlist(support),file.path(out,'fixed_reference_support.csv'))
fwrite(rbindlist(counts),file.path(out,'year5_cause_group_counts.csv'))
fwrite(rbindlist(probs),file.path(out,'consultation_table_probabilities.csv'))
print(rbindlist(counts))
print(rbindlist(support)[landmark %in% c(0,5)])
writeLines(c('Post-exploration descriptive audit, 2026-10-15. No primary model or comparison family changed.','Probability confidence intervals derived from the saved 2,000 primary bootstrap replicates; no new resampling.','Fixed-reference check uses exact diagnosis age x sex x diagnosis year cells, pooled 5-year survivors within each stage. Empty cells diagnose lack of exact empirical support, not failure of every possible model-based standardization.',capture.output(sessionInfo())),file.path(out,'validation_and_session.txt'))
