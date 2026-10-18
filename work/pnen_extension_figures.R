# Figure contract: make selection and inferential limits visible; treatment
# standardization and era contrasts are descriptive extensions, not mechanisms.
# R quantitative grid and cohort-flow diagram, 183 mm, PNG/SVG/600 dpi TIFF.
library(data.table); library(ggplot2);library(patchwork)
out<-'outputs/pnen_manuscript'
theme_set(theme_classic(base_size=9,base_family='Arial')+theme(plot.title=element_text(face='bold',size=10)))
saveplot<-function(p,name,h){
 svglite::svglite(file.path(out,paste0(name,'.svg')),width=183/25.4,height=h/25.4);print(p);dev.off()
 ragg::agg_tiff(file.path(out,paste0(name,'.tiff')),width=183,height=h,units='mm',res=600,compression='lzw');print(p);dev.off()
 ragg::agg_png(file.path(out,paste0(name,'.png')),width=183,height=h,units='mm',res=180);print(p);dev.off()
}
nodes<-data.table(x=c(0,0,0,0,0,0,0),y=7:1,label=c('Pulmonary neuroendocrine export, 2004-2023\n141,383 records','Earliest eligible record per patient\n140,776 patients','Valid survival and reporting information\n140,395 patients','LCNEC or SCLC diagnosed 2004-2015\n68,946 patients','Known initial Summary Stage\n67,305 patients','Localized: 3,652   |   Regional: 14,242\nDistant: 49,411','Treatment analysis, known surgery\nLocalized: 3,613   |   Regional: 14,152'))
ar<-data.table(x=0,xend=0,y=seq(6.68,1.68,-1),yend=seq(6.32,1.32,-1))
flow<-ggplot()+geom_segment(data=ar,aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(length=unit(1.5,'mm')),linewidth=.35)+geom_label(data=nodes,aes(x=x,y=y,label=label),size=3.1,lineheight=1.15,label.padding=unit(.17,'lines'),fill='white',linewidth=.3)+annotate('text',x=3.1,y=c(6.5,5.5,4.5,3.5,1.5),label=c('607 additional records','381 survival/reporting exclusions','71,449 other histologies\nor diagnosis years','1,641 unknown stage','129 unknown surgery\n(localized/regional)'),size=2.8,hjust=.5)+xlim(-2.9,4.7)+ylim(.55,7.5)+theme_void()+labs(caption='The diagnosis-era analysis uses a separate 2016-2018 cohort from the valid export.\nTreatment analyses include localized and regional disease only.')
saveplot(flow,'Figure3_cohort_flow',160);fwrite(nodes,file.path(out,'Figure3_source_data.csv'))
c<-fread('outputs/pnen_extension/extension_contrasts.csv')[method=='IPTW']
a<-c[family=='Treatment stage change 0 to 5']
a[,label:=factor(contrast,levels=c('M2','M1','M0'),labels=c('M2: surgery, RT and chemotherapy','M1: surgery added','M0: age, sex and year'))]
b<-c[family=='Era stage change 0 to 2']
b[,label:=factor(contrast,levels=c('Later minus earlier','2016-2018','2010-2015'),labels=c('Between-era contrast','Diagnosed 2016-2018','Diagnosed 2010-2015'))]
forest<-function(z,title){ggplot(z,aes(estimate*100,label))+geom_vline(xintercept=0,linetype=2,color='grey60')+geom_errorbar(aes(xmin=sim_lower*100,xmax=sim_upper*100),orientation='y',width=.15,color='#216B8A')+geom_point(size=2.5,color='#216B8A')+labs(x='Contrast (percentage points)',y=NULL,title=title)}
p<-forest(a,'a  Treatment standardization: change from year 0 to 5')/forest(b,'b  Diagnosis era: change from year 0 to 2')+plot_annotation(caption='Stage contrasts are localized minus regional; models use the same surgery-known cohort.\nBars: simultaneous 95% intervals within each panel (1,000 bootstrap samples).\nM1 and M2 retain age, sex and year. RT = radiotherapy. No causal treatment effect is estimated.')
saveplot(p,'Figure4_treatment_and_era',145);fwrite(c,file.path(out,'Figure4_source_data.csv'))
