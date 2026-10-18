# R figure contract: absolute subsequent mortality components contextualize
# the stage-specific survivor OS gap; secondary surgery curves remain descriptive.
# 183mm width; SVG, 600dpi TIFF and PNG; no PDF outputs.
source('work/pnen_figures.R')
aj<-fread('outputs/pnen_formal_competing/conditional_competing_differences.csv')[method=='IPTW' & landmark==5 & stage %in% c('Localized','Regional')]
aj<-melt(aj,id.vars=c('stage','outcome'),measure.vars=c('LCNEC','SCLC'),variable.name='Histology',value.name='probability')
aj[,outcome:=factor(outcome,levels=c('unknown_CIF','other_CIF','cancer_CIF','OS'),labels=c('Unknown-cause death','Other-cause death','Cancer-related death','Alive'))]
pd<-ggplot(aj,aes(Histology,probability*100,fill=outcome))+geom_col(width=.6)+facet_wrap(~stage,nrow=1)+scale_fill_manual(values=c('#A8A8A8','#C49850','#BA6844','#83AFC2'))+scale_y_continuous(limits=c(0,100),breaks=c(0,50,100))+labs(x=NULL,y='Probability (%)',fill=NULL,title='d  Subsequent 3-year outcomes among 5-year survivors')+theme(legend.position='bottom',legend.text=element_text(size=7),legend.key.size=grid::unit(3,'mm'),strip.text=element_blank())+guides(fill=guide_legend(nrow=1,reverse=TRUE))
saveplot((p1/p2/pr/pd)+plot_layout(heights=c(1,1,.3,.95))+plot_annotation(caption='Mortality components are absolute probabilities; unknown cause is retained.'),'Figure1_conditional_OS',225)
fwrite(aj,file.path(out,'Figure1_mortality_source_data.csv'))
ss<-fread(file.path(out,'Supplement_extension_extension_estimates.csv'))[job %in% c('Surgery_Recorded_Localized','Surgery_None_Localized') & method=='IPTW']
ss[,surgery:=ifelse(grepl('None',job),'No recorded surgery','Recorded surgery')]
pp<-ggplot(ss,aes(landmark,difference*100))+geom_hline(yintercept=0,linetype=2,color='grey60')+geom_line(color='#216B8A')+geom_point(color='#216B8A')+geom_errorbar(aes(ymin=lower*100,ymax=upper*100),width=.1,color='#216B8A')+facet_wrap(~surgery)+scale_x_continuous(breaks=c(1,3,5))+labs(x='Years already survived',y='LCNEC minus SCLC (percentage points)')
saveplot(pp+plot_annotation(caption='Localized disease; future 3-year OS. No-surgery year 5 suppressed (15 LCNEC survivors).\nTreatment dates unavailable; selection and immortal-time bias may affect comparisons.'),'Supplementary_Figure_S3_surgery',110)
fwrite(ss,file.path(out,'Supplementary_Figure_S3_source_data.csv'))
