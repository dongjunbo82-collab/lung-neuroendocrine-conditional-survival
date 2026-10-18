library(data.table);library(ggplot2);library(patchwork)
base<-Sys.getenv('PNEN_PROJECT',getwd());setwd(base)
out<-'outputs/pnen_manuscript'
os<-fread('outputs/pnen_value_audit/consultation_table_probabilities.csv')
curves<-rbindlist(lapply(c('LCNEC','SCLC'),function(g)os[,.(stage,landmark,Histology=g,probability=get(g),lower=get(paste0(g,'_lower')),upper=get(paste0(g,'_upper')))]))
risk<-fread('outputs/pnen_formal/weight_diagnostics.csv')[method=='IPTW' & stage %in% c('Localized','Regional')]
aj<-fread('outputs/pnen_formal_competing/conditional_competing_differences.csv')[method=='IPTW' & landmark==5 & stage %in% c('Localized','Regional')]
cp<-list()
for(st in c('Localized','Regional')) {
 bb<-readRDS(file.path('outputs/pnen_formal_competing',paste0('bootstrap_',st,'.rds')))
 stopifnot(identical(as.integer(dim(bb)),c(2000L,5L,3L,2L,4L)))
 for(g in c('SCLC','LCNEC'))for(oc in c('cancer_CIF','other_CIF')) {
  j<-match(oc,c('OS','cancer_CIF','other_CIF','unknown_CIF'));gi<-match(g,c('SCLC','LCNEC'));q<-quantile(bb[,5,2,gi,j],c(.025,.975))
  cp[[length(cp)+1L]]<-data.table(stage=st,Histology=g,outcome=oc,probability=aj[stage==st & outcome==oc,get(g)],lower=q[1],upper=q[2])
 }
}
cp<-rbindlist(cp);pal<-c(LCNEC='#216B8A',SCLC='#BA6844')
theme_set(theme_classic(base_family='Arial',base_size=8)+theme(strip.background=element_blank(),strip.text=element_text(face='bold'),plot.title=element_text(size=9,face='bold'),legend.position='bottom',plot.caption=element_text(size=7,hjust=0)))
curves[,stage:=factor(stage,levels=c('Localized','Regional'))];cp[,stage:=factor(stage,levels=c('Localized','Regional'))];risk[,stage:=factor(stage,levels=c('Localized','Regional'))]
p1<-ggplot(curves,aes(landmark,probability*100,color=Histology,shape=Histology,group=Histology))+geom_errorbar(aes(ymin=lower*100,ymax=upper*100),position=position_dodge(.12),width=.1,linewidth=.35)+geom_line(linewidth=.45)+geom_point(position=position_dodge(.12),size=1.6)+facet_wrap(~stage,nrow=1)+scale_color_manual(values=pal)+scale_x_continuous(breaks=c(0,1,2,3,5))+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20))+labs(x='Years already survived',y='Probability of surviving another 3 years (%)',title='a  Historical conditional overall survival',color=NULL,shape=NULL)
pr<-ggplot(risk,aes(landmark,group,label=format(n,big.mark=','),color=group))+geom_text(size=2.5,show.legend=FALSE)+facet_wrap(~stage,nrow=1)+scale_color_manual(values=pal)+scale_x_continuous(breaks=c(0,1,2,3,5),limits=c(-.25,5.3))+labs(x=NULL,y=NULL,title='b  Patients alive and observed at each landmark')+theme_void(base_size=8)+theme(strip.text=element_blank(),axis.text.y=element_text(size=8),plot.title=element_text(size=9,face='bold'))
cp[,outcome:=factor(outcome,levels=c('cancer_CIF','other_CIF'),labels=c('Cancer-related death','Other-cause death'))]
p2<-ggplot(cp,aes(outcome,probability*100,color=Histology,shape=Histology))+geom_errorbar(aes(ymin=lower*100,ymax=upper*100),position=position_dodge(.4),width=.15,linewidth=.4)+geom_point(position=position_dodge(.4),size=2)+geom_text(aes(y=upper*100+1.5,label=sprintf('%.1f',probability*100)),position=position_dodge(.4),size=2.3,show.legend=FALSE)+facet_wrap(~stage,nrow=1)+scale_color_manual(values=pal)+scale_y_continuous(limits=c(0,35),breaks=seq(0,35,5))+labs(x=NULL,y='Death probability over the next 3 years (%)',title='c  Mortality components among 5-year survivors',color=NULL,shape=NULL)+guides(color='none',shape='none')
fig<-(p1/pr/p2)+plot_layout(heights=c(1,.3,.95),guides='collect')+plot_annotation(caption='Historical SEER diagnoses, 2004-2015. C0 demographic standardization; pointwise 95% intervals.\nGroup-level reference estimates; unknown-cause deaths remain a separate event in the analysis.')
svglite::svglite(file.path(out,'Figure1_conditional_OS.svg'),width=183/25.4,height=200/25.4);print(fig);dev.off()
ragg::agg_tiff(file.path(out,'Figure1_conditional_OS.tiff'),width=183,height=200,units='mm',res=600,compression='lzw');print(fig);dev.off()
ragg::agg_png(file.path(out,'Figure1_conditional_OS.png'),width=183,height=200,units='mm',res=180);print(fig);dev.off()
fwrite(curves,file.path(out,'Figure1_probability_interval_source_data.csv'));fwrite(cp,file.path(out,'Figure1_mortality_interval_source_data.csv'))
