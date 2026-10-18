# Figure contract: stage-specific differences in conditional OS narrow more in
# localized than regional disease; distant late-survivor uncertainty stays visible.
# Quantitative grid; R backend; 183mm width; SVG and600dpiTIFF with PNG previews.
library(data.table);library(ggplot2);library(patchwork)
out<-'outputs/pnen_manuscript';dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-fread('outputs/pnen_formal/adjusted_OS_differences.csv')[method=='IPTW']
d_all<-copy(d)
d<-d[stage!='Distant']
d[,stage:=factor(stage,levels=c('Localized','Regional','Distant'))]
theme_set(theme_classic(base_size=9,base_family='Arial')+theme(strip.background=element_blank(),strip.text=element_text(face='bold'),legend.position='bottom',plot.title=element_text(size=10,face='bold')))
pal<-c('LCNEC'='#216B8A','SCLC'='#BA6844')
long<-melt(d,id.vars=c('stage','landmark'),measure.vars=c('LCNEC','SCLC'),variable.name='Histology',value.name='OS')
p1<-ggplot(long,aes(landmark,OS*100,color=Histology,shape=Histology))+geom_line(linewidth=.5)+geom_point(size=2)+facet_wrap(~stage,nrow=1)+scale_color_manual(values=pal)+scale_x_continuous(breaks=c(0,1,2,3,5))+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20))+labs(x='Years already survived',y='Conditional 3-year overall survival (%)',title='a  Standardized survival probabilities')
p2<-ggplot(d,aes(landmark,difference*100))+geom_hline(yintercept=0,linetype=2,color='grey60')+geom_line(color='#216B8A',linewidth=.5)+geom_errorbar(aes(ymin=lower*100,ymax=upper*100),width=.12,linewidth=.4)+geom_point(color='#216B8A',size=2)+facet_wrap(~stage,nrow=1)+scale_x_continuous(breaks=c(0,1,2,3,5))+scale_y_continuous(limits=c(-25,35),breaks=seq(-20,30,10))+labs(x='Years already survived',y='LCNEC minus SCLC (percentage points)',title='b  Absolute differences with pointwise 95% intervals')
saveplot<-function(p,name,h) {
 svglite::svglite(file.path(out,paste0(name,'.svg')),width=183/25.4,height=h/25.4);print(p);dev.off()
 ragg::agg_tiff(file.path(out,paste0(name,'.tiff')),width=183,height=h,units='mm',res=600,compression='lzw');print(p);dev.off()
 ragg::agg_png(file.path(out,paste0(name,'.png')),width=183,height=h,units='mm',res=180);print(p);dev.off()
}
risk<-fread('outputs/pnen_formal/weight_diagnostics.csv')[method=='IPTW']
risk_all<-copy(risk)
risk<-risk[stage!='Distant']
risk[,stage:=factor(stage,levels=c('Localized','Regional','Distant'))]
risk[,group:=factor(group,levels=c('SCLC','LCNEC'))]
pr<-ggplot(risk,aes(landmark,group,label=n,color=group))+geom_text(size=2.5,show.legend=FALSE)+facet_wrap(~stage,nrow=1)+scale_color_manual(values=pal)+scale_x_continuous(breaks=c(0,1,2,3,5),limits=c(-.25,5.3))+labs(x=NULL,y=NULL,title='c  Patients alive and under observation at each landmark')+theme_void(base_size=9)+theme(strip.text=element_blank(),axis.text.y=element_text(size=8),plot.title=element_text(size=10,face='bold'))
saveplot((p1/p2/pr)+plot_layout(heights=c(1,1,.3))+plot_annotation(caption='Complete three-stage estimates are provided in Supplementary Figure S1.'),'Figure1_conditional_OS',182)
c<-fread('outputs/pnen_formal/cross_stage_change_contrasts.csv')
c[,contrast:=factor(paste(stage_A,'vs',stage_B),levels=rev(paste(stage_A,'vs',stage_B)))]
p3<-ggplot(c,aes(estimate*100,contrast))+geom_vline(xintercept=0,linetype=2,color='grey60')+geom_errorbar(aes(xmin=sim_lower*100,xmax=sim_upper*100),orientation='y',width=.12,linewidth=.7,color='#216B8A')+geom_point(size=2.5,color='#216B8A')+labs(x='Between-stage contrast in change from landmark 0 to 5 (percentage points)',y=NULL)+theme(axis.text.y=element_text(size=10))
sens<-fread('outputs/pnen_formal_sensitivity/change_contrasts.csv')[stage_A=='Localized' & stage_B=='Regional' & ((method=='IPTW') | (scenario=='main' & method=='overlap'))]
sens[,label:=fcase(scenario=='main' & method=='IPTW','Main cohort',scenario=='main' & method=='overlap','Overlap weighting',grepl('2010|recent',scenario),'Diagnosed 2010-2015',grepl('single',scenario),'Single recorded primary',default='Histologically confirmed')]
sens[,label:=factor(label,levels=rev(unique(label)))]
p4<-ggplot(sens,aes(estimate*100,label))+geom_vline(xintercept=0,linetype=2,color='grey60')+geom_errorbar(aes(xmin=sim_lower*100,xmax=sim_upper*100),orientation='y',width=.12,linewidth=.7,color='#216B8A')+geom_point(size=2.5,color='#216B8A')+labs(x='Localized minus regional change (percentage points)',y=NULL,title='b  Sensitivity analyses of the focused contrast')
saveplot((p3+labs(title='a  Comparisons between initial stages'))/p4,'Figure2_stage_change_contrasts',160)
fwrite(d,file.path(out,'Figure1_source_data.csv'));fwrite(risk,file.path(out,'Figure1_risk_sets.csv'));fwrite(c,file.path(out,'Figure2_source_data.csv'));fwrite(sens,file.path(out,'Figure2_sensitivity_source_data.csv'))

# Full three-stage display retained without altering estimates or intervals.
d_all[,stage:=factor(stage,levels=c('Localized','Regional','Distant'))]
risk_all[,stage:=factor(stage,levels=c('Localized','Regional','Distant'))]
risk_all[,group:=factor(group,levels=c('SCLC','LCNEC'))]
long_all<-melt(d_all,id.vars=c('stage','landmark'),measure.vars=c('LCNEC','SCLC'),variable.name='Histology',value.name='OS')
saveplot(((p1 %+% long_all)/(p2 %+% d_all)/(pr %+% risk_all))+plot_layout(heights=c(1,1,.3))+plot_annotation(caption='Distant-stage LCNEC: only 65 patients at year 5; late estimates are imprecise.'),'Supplementary_Figure_S1_all_stages',182)
fwrite(d_all,file.path(out,'Supplementary_Figure_S1_source_data.csv'))
fwrite(risk_all,file.path(out,'Supplementary_Figure_S1_risk_sets.csv'))
