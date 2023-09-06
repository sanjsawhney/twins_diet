# twins_diet
Code used for data visualization and analysis in Sawhney et al (Weaning accelerates and transforms within-host adaptation in the infant gut microbiome).

# Figure 3A

Made in Biorender.com

# Figure 3B

<pre>
  #--------------------------220330 MQ Winner vs. HQ Winner vs. Ref dRep 0.98  MAG Completeness-Contamination Scatter Plot and Histograms--------------------------
#Read in CSV
df_All_MAG_metrics<-read.csv('220330_All_MAGs_LQ-vs-MQ-vs-HQ_dRep_0.98_metrics.csv',
                              sep=",",
                              header = T)
df_All_MAG_metrics$Group_ScatterPlot<-factor(df_All_MAG_metrics$Group_ScatterPlot, levels=c("HQ","MQ", "LQ"))
df_All_MAG_metrics$Completeness<-as.numeric(df_All_MAG_metrics$Completeness)
df_All_MAG_metrics$Contamination<-as.numeric(df_All_MAG_metrics$Contamination)


library(ggplot2)

#Completeness vs. Contamination Scatter plot
scatterplot_All_MAG_completeness_contamination<-ggplot(df_All_MAG_metrics, aes(x=Completeness, y=Contamination))+
  geom_point(size=0.5)+
  scale_x_continuous(limits=c(0,100),expand = c(0,1))+
  scale_y_continuous(limits=c(0,100),expand = c(0,1))+
  scale_fill_manual(values=c("#3399FF", "#33FF99", "#FF9900"))+
  #Remove fill= if not showing grouped boxplot
  labs(x="Completeness (%)", y="Contamination (%)")+
  theme_bw()+
  theme(legend.position='none', plot.title = element_text(face="bold",size=20,hjust=0.5),axis.title = element_text(face="bold", size=21),axis.text = element_text(size = 14, face="bold"),plot.margin = margin(5, 10, 5, 5))

scatterplot_All_MAG_completeness_contamination

histogram_All_MAG_completeness <- ggplot(df_All_MAG_metrics, aes(x=Completeness,fill=Group_ScatterPlot, color=Group_ScatterPlot))+
  geom_histogram(color="black")+
  scale_fill_manual(values=c("#33FF99","#3399FF","#FF6666"))+
  labs(x="Completeness (%)",y="Count")+
  theme_classic()+
  theme(legend.position="none",legend.title=element_blank(),axis.title.y = element_text(face="bold", size=19),legend.text = element_text(size=14, face="bold"), axis.title.x = element_blank(),axis.text.y = element_text(size=12,face="bold"),axis.text.x=element_blank())

histogram_All_MAG_completeness

histogram_All_MAG_contamination <- ggplot(df_All_MAG_metrics, aes(y=Contamination,fill=Group_ScatterPlot))+
  geom_histogram(color="black")+
  #  scale_fill_discrete(breaks=c('LQ', 'MQ', 'HQ'))+
  scale_fill_manual(values=c("#33FF99","#3399FF","#FF6666"))+
  labs(y="Contamination (%)",x="Count")+
  theme_classic()+
  theme(legend.position="none",legend.title=element_blank(),axis.title.x = element_blank(),legend.text = element_text(size=14, face="bold"), axis.title.y = element_blank(),axis.text.x = element_text(size=16,face="bold",angle=270,hjust=0,vjust=.5),axis.text.y=element_blank())

histogram_All_MAG_contamination

</pre>

# Figure 3C

<pre>

#--------------------------220329 MQ Winner vs. HQ Winner vs. Ref dRep 0.98  MAG Completeness--------------------------

  #Read in CSV
df_dRep_MAG_metrics<-read.csv('220329_WinnersMQ-vs-WinnersHQ-vs-Refs_dRep_0.98_metrics.csv',
                              sep=",",
                              header = T)
df_dRep_MAG_metrics$Group<-factor(df_dRep_MAG_metrics$Group, levels=c("MQ Winners", "HQ Winners", "Reference Genomes"))
df_dRep_MAG_metrics$Completeness<-as.numeric(df_dRep_MAG_metrics$Completeness)
df_dRep_MAG_metrics$Contamination<-as.numeric(df_dRep_MAG_metrics$Contamination)
df_dRep_MAG_metrics$Strain_Heterogeneity<-as.numeric(df_dRep_MAG_metrics$Strain_Heterogeneity)
df_dRep_MAG_metrics$N50<-as.numeric(df_dRep_MAG_metrics$N50)
df_dRep_MAG_metrics$Contig_count<-as.numeric(df_dRep_MAG_metrics$Contig_count)

library(ggplot2)
library(ggpubr)
library(rstatix)

#Check for normality - all are non-parametric distribution
shapiro.test(df_dRep_MAG_metrics$Completeness)

ggdensity(df_dRep_MAG_metrics$Completeness)
ggqqplot(df_dRep_MAG_metrics$Completeness)

ggdensity(df_dRep_MAG_metrics$Contig_count)
ggqqplot(df_dRep_MAG_metrics$Contig_count)

ggdensity(df_dRep_MAG_metrics$N50)
ggqqplot(df_dRep_MAG_metrics$N50)

ggdensity(df_dRep_MAG_metrics$Contamination)
ggqqplot(df_dRep_MAG_metrics$Contamination)

ggdensity(df_dRep_MAG_metrics$Strain_Heterogeneity)
ggqqplot(df_dRep_MAG_metrics$Strain_Heterogeneity)

#Significance comparisons
my_comparisons <- list( c("MQ Winners", "HQ Winners"), c("HQ Winners", "Reference Genomes"), c("MQ Winners", "Reference Genomes"))
#my_comparisons <- list(c("HQ Winners", "Reference Genomes"), c("MQ Winners", "Reference Genomes"))

#Completeness Boxplot (change or remove fill= to bin by SR vs LR, Plate vs Stool, no bin)
#boxplot_dRep_MAG_completeness<-ggplot(df_dRep_MAG_metrics, aes(x=Group, y=Completeness, fill=Assembly_Type))+
boxplot_dRep_MAG_completeness<-ggplot(df_dRep_MAG_metrics, aes(x=Group, y=Completeness, fill=Group))+  
  #remove position=position_dodge() if not showing grouped boxplot
  #geom_boxplot(width=0.5, outlier.size=0.3, position = position_dodge(0.65))+
  geom_boxplot(width=0.5, outlier.size=0.3, lwd=0.7)+
  scale_x_discrete(labels = c("MQ","HQ","Ref."))+
  scale_y_continuous(limits=c(30,107), expand = c(0,0), breaks=c(40,60,80,100))+
  scale_fill_manual(values=c("#3399FF", "#33FF99", "#FF9900"))+
  theme(axis.title.y=element_blank())+
  #Remove fill= if not showing grouped boxplot
  #labs(y="Completeness (%)", fill="Assembly Type")+
  labs(y="Completeness (%)")+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE, label = "p.signif", test = "kruskal.test", tip.length = 0.005, label.y = c(98,99.5,101.5), vjust=.75, size =5)+
  theme(legend.position='none',panel.border = element_rect(color = "black", size=1),axis.ticks = element_line(colour = "black", size = 0.7), axis.title = element_blank(),axis.text = element_text(size = 24, face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

boxplot_dRep_MAG_completeness

#--------------------------220329 MQ Winner vs. HQ Winner vs. Ref dRep 0.98  MAG Contamination--------------------------
#Read in CSV
df_dRep_MAG_metrics<-read.csv('220329_WinnersMQ-vs-WinnersHQ-vs-Refs_dRep_0.98_metrics.csv',
                              sep=",",
                              header = T)
df_dRep_MAG_metrics$Group<-factor(df_dRep_MAG_metrics$Group, levels=c("MQ Winners", "HQ Winners", "Reference Genomes"))
df_dRep_MAG_metrics$Completeness<-as.numeric(df_dRep_MAG_metrics$Completeness)
df_dRep_MAG_metrics$Contamination<-as.numeric(df_dRep_MAG_metrics$Contamination)
df_dRep_MAG_metrics$Strain_Heterogeneity<-as.numeric(df_dRep_MAG_metrics$Strain_Heterogeneity)
df_dRep_MAG_metrics$N50<-as.numeric(df_dRep_MAG_metrics$N50)
df_dRep_MAG_metrics$Contig_count<-as.numeric(df_dRep_MAG_metrics$Contig_count)

library(ggplot2)
library(ggpubr)

#Contamination Boxplot (change or remove fill= to bin by SR vs LR, Plate vs Stool, no bin)
boxplot_dRep_MAG_contamination<-ggplot(df_dRep_MAG_metrics, aes(x=Group, y=Contamination, fill=Group))+  
  #remove position=position_dodge() if not showing grouped boxplot
  #geom_boxplot(width=0.5, outlier.size=0.3, position = position_dodge(0.65))+
  geom_boxplot(width=0.5, outlier.size=0.3, lwd=0.7)+
  scale_x_discrete(labels = c("MQ","HQ","Ref."))+
  scale_y_continuous(limits=c(0,12.55), expand = c(0,0.2), breaks=c(0.0,2.5,5.0,7.5,10.0,12.5))+
  scale_fill_manual(values=c("#3399FF", "#33FF99", "#FF9900"))+
  theme(axis.title.y=element_blank())+
  #Remove fill= if not showing grouped boxplot
  #labs(y="Contamination (%)", fill="Assembly Type")+
  labs(y="Contamination (%)")+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.signif", test = "kruskal.test", tip.length = 0.005, label.y = c(11,11.4,11.8), vjust=.75, size =5)+
  theme(legend.position='none',panel.border = element_rect(color = "black", size=1),axis.ticks = element_line(colour = "black", size = 0.7), axis.title = element_blank(),axis.text = element_text(size = 24, face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

boxplot_dRep_MAG_contamination

#--------------------------220329 MQ Winner vs. HQ Winner vs. Ref dRep 0.98  MAG Strain Heterogeneity--------------------------
#Read in CSV
df_dRep_MAG_metrics<-read.csv('220329_WinnersMQ-vs-WinnersHQ-vs-Refs_dRep_0.98_metrics.csv',
                              sep=",",
                              header = T)
df_dRep_MAG_metrics$Group<-factor(df_dRep_MAG_metrics$Group, levels=c("MQ Winners", "HQ Winners", "Reference Genomes"))
df_dRep_MAG_metrics$Completeness<-as.numeric(df_dRep_MAG_metrics$Completeness)
df_dRep_MAG_metrics$Contamination<-as.numeric(df_dRep_MAG_metrics$Contamination)
df_dRep_MAG_metrics$Strain_Heterogeneity<-as.numeric(df_dRep_MAG_metrics$Strain_Heterogeneity)
df_dRep_MAG_metrics$N50<-as.numeric(df_dRep_MAG_metrics$N50)
df_dRep_MAG_metrics$Contig_count<-as.numeric(df_dRep_MAG_metrics$Contig_count)

library(ggplot2)
library(ggpubr)

#Strain Heterogeneity Boxplot
boxplot_dRep_MAG_het<-ggplot(df_dRep_MAG_metrics, aes(x=Group, y=Strain_Heterogeneity,fill=Group))+
  #remove position=position_dodge() if not showing grouped boxplot
  #geom_boxplot(width=0.5, outlier.size=0.3, position = position_dodge(0.65))+
  geom_boxplot(width=0.5, outlier.size=0.3, lwd=0.7)+
  scale_x_discrete(labels = c("MQ","HQ","Ref."))+
  scale_y_continuous(limits=c(0,100),expand = c(0,1))+
  scale_y_continuous(limits=c(0,110), expand = c(0,1), breaks=c(0,25,50,75,100))+
  scale_fill_manual(values=c("#3399FF", "#33FF99", "#FF9900"))+
  theme(axis.title.y=element_blank())+
  #Remove fill= if not showing grouped boxplot
  #labs(y="Strain Heterogeneity (%)", fill="Assembly Type")+
  labs(y="Strain Het. (%)")+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE, label = "p.signif", test = "kruskal.test", tip.length = 0.005, label.y = c(97,100,103), vjust=.75, size =5)+
  theme(legend.position='none',panel.border = element_rect(color = "black", size=1),axis.ticks = element_line(colour = "black", size = 0.7), axis.title = element_blank(),axis.text = element_text(size = 24, face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

boxplot_dRep_MAG_het
  
#--------------------------220329 MQ Winner vs. HQ Winner vs. Ref dRep 0.98  MAG N50--------------------------
#Read in CSV
df_dRep_MAG_metrics<-read.csv('220329_WinnersMQ-vs-WinnersHQ-vs-Refs_dRep_0.98_metrics.csv',
                              sep=",",
                              header = T)
df_dRep_MAG_metrics$Group<-factor(df_dRep_MAG_metrics$Group, levels=c("MQ Winners", "HQ Winners", "Reference Genomes"))
df_dRep_MAG_metrics$Completeness<-as.numeric(df_dRep_MAG_metrics$Completeness)
df_dRep_MAG_metrics$Contamination<-as.numeric(df_dRep_MAG_metrics$Contamination)
df_dRep_MAG_metrics$Strain_Heterogeneity<-as.numeric(df_dRep_MAG_metrics$Strain_Heterogeneity)
df_dRep_MAG_metrics$N50<-as.numeric(df_dRep_MAG_metrics$N50)
df_dRep_MAG_metrics$Contig_count<-as.numeric(df_dRep_MAG_metrics$Contig_count)

library(ggplot2)
library(ggpubr)

#N50 boxplot
library(scales)

boxplot_dRep_MAG_N50<-ggplot(df_dRep_MAG_metrics, aes(x=Group, y=log10(N50), fill=Group))+
  geom_boxplot(width=0.5, outlier.size=0.3, lwd=0.7)+
  scale_x_discrete(labels = c("MQ","HQ","Ref."))+
  scale_fill_manual(values=c("#3399FF", "#33FF99", "#FF9900"))+
  scale_y_continuous(breaks = c(3,4,5,6,7),
                     limits=c(3,7.15), expand=c(0,0.09),
                     labels = function(lab) {
                       do.call( expression,lapply(paste(lab), function(x) bquote(bold("10"^.(x)))))
                     })+
 # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
 #               labels = trans_format("log10", math_format(10^.x)),
 #               limits=c(1e3,1e7),expand = c(0,0))+
  annotation_logticks(sides="l")+
  theme_bw()+
  labs(y="N50 (log)")+
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.signif", test = "kruskal.test", tip.length = 0.005, label.y = c(6.7,6.82,6.94), vjust=.75, size =5)+
  theme(legend.position='none',panel.border = element_rect(color = "black", size=1),axis.ticks = element_line(colour = "black", size = 0.7), axis.title = element_blank(),axis.text = element_text(size = 24, face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

boxplot_dRep_MAG_N50
  
#--------------------------220329 MQ Winner vs. HQ Winner vs. Ref dRep 0.98  MAG Contig Count--------------------------
#Read in CSV
df_dRep_MAG_metrics<-read.csv('220329_WinnersMQ-vs-WinnersHQ-vs-Refs_dRep_0.98_metrics.csv',
                              sep=",",
                              header = T)
df_dRep_MAG_metrics$Group<-factor(df_dRep_MAG_metrics$Group, levels=c("MQ Winners", "HQ Winners", "Reference Genomes"))
df_dRep_MAG_metrics$Completeness<-as.numeric(df_dRep_MAG_metrics$Completeness)
df_dRep_MAG_metrics$Contamination<-as.numeric(df_dRep_MAG_metrics$Contamination)
df_dRep_MAG_metrics$Strain_Heterogeneity<-as.numeric(df_dRep_MAG_metrics$Strain_Heterogeneity)
df_dRep_MAG_metrics$N50<-as.numeric(df_dRep_MAG_metrics$N50)
df_dRep_MAG_metrics$Contig_count<-as.numeric(df_dRep_MAG_metrics$Contig_count)

library(ggplot2)
library(ggpubr)
library(scales)

#Contig Count Boxplot
boxplot_dRep_MAG_contig_count<-ggplot(df_dRep_MAG_metrics, aes(x=Group, y=log10(Contig_count), fill=Group))+  
  geom_boxplot(width=0.5, outlier.size=0.3,lwd=0.7)+
  scale_x_discrete(labels = c("MQ","HQ","Ref."))+
  scale_y_continuous(breaks = c(0,1,2,3),
                     limits=c(0,3.48), expand=c(0,0.09),
                     labels = function(lab) {
                       do.call( expression,lapply(paste(lab), function(x) bquote(bold("10"^.(x)))))
                     })+
  annotation_logticks(sides="l")+
  scale_fill_manual(values=c("#3399FF", "#33FF99", "#FF9900"))+
  theme(axis.title.y=element_blank())+
  #Remove fill= if not showing grouped boxplot
  #labs(y="Contamination (%)", fill="Assembly Type")+
  labs(y="# contigs (log)")+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.signif", test = "kruskal.test", tip.length = 0.005, label.y = c(3.1,3.2,3.3), vjust=.75, size =5)+
  theme(legend.position='none',panel.border = element_rect(color = "black", size=1),axis.ticks = element_line(colour = "black", size = 0.7), axis.title = element_blank(),axis.text = element_text(size = 24, face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

boxplot_dRep_MAG_contig_count

</pre>

# Figure 4A

Made in https://itol.embl.de/

# Figure 4B

<pre>
#------------------------Figure 4B: Infant vs. Mother: Stool, RG, Persisting RG------------------------
library(ggplot2)
library(ggprism)

df_infant_vs_mother<-read.csv('230412_Infant_vs_Mother_Stool_RG_Persisting.csv',
                                        sep=",",
                                        header = T)


df_infant_vs_mother$Group <- factor(df_infant_vs_mother$Group, levels = c("Mother","Infant"))
df_infant_vs_mother$Comparison <- factor(df_infant_vs_mother$Comparison, levels = c("Stool","RGs","Persisting RGs"))

#Plot
ggplot(data=df_infant_vs_mother, aes(x=Comparison, y=Percentage, fill=Group)) +
  geom_bar(stat="identity",width=0.45, color="black")+
  theme_bw()+
 # facet_wrap( ~ Comparison)+
  scale_fill_manual(values=c("black", "#F9DA78"))+
  scale_x_discrete(labels=c("Stool" = "Stool\n(n=214)", "RGs" = "RGs\n(n=3,995)","Persisting RGs" = "Persisting\nRGs (n=1,093)"))+
  scale_y_continuous(limits=c(0,100), expand=c(0,0), guide = guide_prism_minor())+
  ylab("Percentage")+
  theme(axis.text = element_text(size=11, face="bold"), plot.title = element_text(size=13, face="bold", hjust=0.5), axis.title.y = element_text(size=13, face="bold"), axis.title.x = element_blank(), legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"), legend.position="bottom", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

</pre>

# Supplementary Figure 1

<pre>

#------------------------Suppl Figure 1: Sampling schematic------------------------

library(ggplot2)

df_sampling<-read.csv('230802_sampling_timelines.csv',
                      sep=",",
                      header = T)

df_sampling$MOL <- as.numeric(df_sampling$MOL)
df_sampling$Individual <- as.factor(df_sampling$Individual)
df_sampling$Individual <- factor(df_sampling$Individual, levels = c("48-2","48-1","47-2","47-1","45-2","45-1","44-2","44-1","43-2","43-1","40-2","40-1","39-2","39-1","37-2","37-1","30-2","30-1","29-2","29-1","28-2","28-1","27-2","27-1","25-2","25-1","24-2","24-1","21-2","21-1","20-2","20-1","19-2","19-1","18-2","18-1","17-2","17-1","16-2","16-1","14-2","14-1","13-2","13-1","12-2","12-1","10-2","10-1","08-2","08-1","06-2","06-1","C048","C047","C045","C044","C043","C040","C039","C037","C030","C029","C028","C027","C025","C024","C021","C020","C019","C018","C017","C016","C014","C013","C012","C010","C008","C006"))
df_sampling$Cohort_sequencing <- factor(df_sampling$Cohort_sequencing, levels = c("shallow shotgun","extended"))

#Plot
ggplot(data=df_sampling, aes(x=MOL, y=Individual)) +
  geom_line()+
  geom_point(shape=21, aes(color=Cohort_sequencing, fill=Cohort_sequencing))+
  scale_color_manual(values=c("black", "black"))+
  scale_fill_manual(values=c("gray", "black"))+
  xlab("Months post childbirth")+
  guides(fill=guide_legend(title="Sequencing type"))+
  theme_bw()+
  theme(axis.text.x = element_text(size=14, face="bold"), axis.text.y = element_text(size=8, face="bold"), plot.title=element_blank(), axis.title=element_text(size=14, face="bold"), panel.border = element_rect(colour = "black", size=0.75))

</pre>
  
# Supplementary Figure 4
<pre>
library(ggplot2)

twins_results <- read.delim("210722_MAG_evaluation_nopilon.txt", header = TRUE,check.names = FALSE, dec=".")
twins_results$method <- factor(twins_results$method, levels = c("alignedFracs", "highcov", "cmb", "rawMSpades", "mSpadesNanopore", "mSpades_alignedFracs", "mSpades_highcov", "mSpades_cmb", "mSpadesNanopore_alignedFracs", "mSpadesNanopore_highcov", "mSpadesNanopore_cmb", "mSpades_mFlye", "mSpadesNanopore_mFlye","OPERA-MS"))


HQ_MAGs<-ggplot(twins_results, aes(x=method, y=HQ, fill=method)) +
  theme_bw()+
  geom_boxplot()+
  scale_fill_manual(values=c("#0033CC","#6699FF","#33CCFF","#669933","#339900","#CC0000","#FF3333","#FF6666","#FF99FF","#CC99FF","#9933FF","#999999","#666666","#FF9933"))+
  theme(legend.text=element_blank(),legend.title=element_blank(),legend.position="none",
        plot.margin = margin(10, 40, 10, 10),
        text = element_text(size=15,face = "bold"),axis.text.x = element_text(angle = 320, hjust = 0),axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black")) +
  ylab("Count - High Quality (HQ) MAGs")
HQ_MAGs

HQMQ_MAGs<-ggplot(twins_results, aes(x=method, y=`# >MQ MAGS`, fill=method)) +
  theme_bw()+
  geom_boxplot()+
  scale_fill_manual(values=c("#0033CC","#6699FF","#33CCFF","#669933","#339900","#CC0000","#FF3333","#FF6666","#FF99FF","#CC99FF","#9933FF","#999999","#666666","#FF9933"))+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"), 
        plot.margin = margin(10, 40, 10, 10),
        text = element_text(size=15,face = "bold"),axis.text.x = element_text(angle = 320, hjust = 0),axis.title.x = element_blank())+
        ylab("Count - High Quality + Medium Quality (MQ) MAGs")
HQMQ_MAGs

MQlength_MAGs<-ggplot(twins_results, aes(x=method, y=as.numeric(`>MQ length`), fill=method)) +
  theme_bw()+
  geom_boxplot()+
  scale_fill_manual(values=c("#0033CC","#6699FF","#33CCFF","#669933","#339900","#CC0000","#FF3333","#FF6666","#FF99FF","#CC99FF","#9933FF","#999999","#666666"))+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"), 
        text = element_text(size=15,face = "bold"),axis.text.x = element_text(angle = 320, hjust = 0),axis.title.x = element_blank())+
  ylab("Total Length of (HQ and MQ MAGs)")
MQlength_MAGs

MQcontigs_MAGs<-ggplot(twins_results, aes(x=method, y=as.numeric(`#contigs`), fill=method)) +
  theme_bw()+
  geom_boxplot()+
  scale_fill_manual(values=c("#0033CC","#6699FF","#33CCFF","#669933","#339900","#CC0000","#FF3333","#FF6666","#FF99FF","#CC99FF","#9933FF","#999999","#666666","#FF9933"))+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"), 
        plot.margin = margin(10, 40, 10, 10),
        text = element_text(size=15,face = "bold"),axis.text.x = element_text(angle = 320, hjust = 0),axis.title.x = element_blank())+
  ylab("Average contigs per MAG (HQ and MQ MAGs)")
MQcontigs_MAGs

MQN50_MAGs<-ggplot(twins_results, aes(x=method, y=`>MQ N50`, fill=method)) +
  theme_bw()+
  geom_boxplot()+
  scale_fill_manual(values=c("#0033CC","#6699FF","#33CCFF","#669933","#339900","#CC0000","#FF3333","#FF6666","#FF99FF","#CC99FF","#9933FF","#999999","#666666","#FF9933"))+
  theme(legend.text=element_text(size=12, face="bold"),legend.title=element_text(size=14, face="bold"),legend.position="none",
        panel.grid.minor = element_blank(),panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"), 
        plot.margin = margin(10, 40, 10, 10),
        text = element_text(size=15,face = "bold"),axis.text.x = element_text(angle = 320, hjust = 0),axis.title.x = element_blank())+
  ylab("Average N50 per MAG (HQ and MQ MAGs)")
MQN50_MAGs

</pre>

# Supplementary Figure 5A

<pre>

  #--------------------------Suppl Figure 5A: 220227 pre dRep MOL MAG Count--------------------------
library(ggplot2)
library(RColorBrewer)

#Read in csv
df_MOL_count <- read.csv('220227_pre-dRep_MOL_MAG_count.csv',
                         sep=",",
                         header = T)

#Plot
ggplot(data=df_MOL_count, aes(x=MOL,y=Count))+
  geom_point()+
  geom_smooth(span=1)+
  theme_bw()+
  ggtitle("Putative Genome Count per Sample by MOL") +
  xlab("MOL")+
  ylab("Putative Genome Count")+
  theme(plot.title = element_text(face="bold",hjust = 0.5))+
  theme(axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"))

</pre>

# Supplementary Figure 5B

<pre>
#--------------------------Suppl Figure 5B: 220227 pre dRep Timepoint MAG Count (Infant Fraction)--------------------------
#Read in CSV
df_MAG_timepoint_infants<-read.csv('220227_pre-dRep_MAG_Fraction_Timepoint_Infants.csv',
                             sep=",",
                             header = T)
df_MAG_timepoint_infants$Timepoint<-factor(df_MAG_timepoint_infants$Timepoint, levels=c("1","2","3","4","5","6"))

barplot_MAG_timepoint_infants <- ggplot(data=df_MAG_timepoint_infants, aes(x=Timepoint, y=Fraction)) +
  geom_bar(stat="identity",fill="#CCCCFF", color="black", width=0.7)+
  geom_errorbar(aes(ymin=Fraction-StdDev, ymax=Fraction+StdDev), width=.2)+
  scale_y_continuous(limits=c(0,.31),expand = c(0,0))+
  labs(title="Infant",y="Fraction of All (pre-dRep)\nMAGs Per Individual")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size=14),axis.text.x = element_text(size = 12, face="bold"),plot.title = element_text(face="bold",size=16,hjust = 0.5))

barplot_MAG_timepoint_infants

  </pre>
  
# Supplementary Figure 5C

<pre>
#--------------------------Suppl Figure 5C: 220227 pre dRep Timepoint MAG Count (Maternal Fraction)--------------------------
#Read in CSV
df_MAG_timepoint_mothers<-read.csv('220227_pre-dRep_MAG_Fraction_Timepoint_Mothers.csv',
                           sep=",",
                           header = T)
df_MAG_timepoint_mothers$Timepoint<-factor(df_MAG_timepoint_mothers$Timepoint, levels=c("1","2","3","4"))

barplot_MAG_timepoint_mothers <- ggplot(data=df_MAG_timepoint_mothers, aes(x=Timepoint, y=Fraction)) +
  geom_bar(stat="identity",fill="#CCCCFF", color="black", width=0.7)+
  geom_errorbar(aes(ymin=Fraction-StdDev, ymax=Fraction+StdDev), width=.2)+
  scale_y_continuous(limits=c(0,.4),expand = c(0,0))+
  labs(title="Maternal",y="Fraction of All (pre-dRep)\nMAGs Per Individual")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size=14),axis.text.x = element_text(size = 12, face="bold"),plot.title = element_text(face="bold",size=16,hjust = 0.5))

barplot_MAG_timepoint_mothers

</pre>
  
# Supplementary Figure 5D

<pre>
#--------------------------Suppl Figure 5D: 220301 pre-dRep MAG Count by # of Samples (All vs Filtered)--------------------------
#Read in CSV
df_pre_dRep_MAG_ct_vs<-read.csv('220301_pre-dRep_MAG_Count_Filtered.csv',
                             sep=",",
                             header = T)
df_pre_dRep_MAG_ct_vs$Timepoints<-factor(df_pre_dRep_MAG_ct_vs$Timepoints, levels=c("3", "4", "6-7","19-21"))

library(ggplot2)
library(ggpubr)

#Plot
boxplot_pre_dRep_MAG_ct_vs<-ggplot(df_pre_dRep_MAG_ct_vs, aes(x=Timepoints, y=MAG_Count, fill=Group))+
  #remove position=position_dodge() if not showing grouped boxplot
  geom_boxplot(width=0.5, outlier.shape=NA, position = position_dodge(0.65))+
  #geom_boxplot(width=0.5, outlier.size=0.3)+
  scale_x_discrete(labels = c("3\n(2 mothers)","4\n(8 mothers)","6-7\n(16 infants)","19-21\n(4 infants)"))+
  scale_fill_manual(values=c("#F8766D", "#00BA38"))+
  scale_y_continuous(limits=c(0,1500),expand = c(0,0))+
  geom_jitter(shape=16, size =0.85, position=position_jitterdodge(.2))+
  #Only include for Assembly Source grouping
  theme(axis.title.y=element_blank())+
  #Remove fill= if not showing grouped boxplot
  labs(x="Number of Samples",y="MAG count per individual")+
  theme_bw()+
  theme(legend.position='top', axis.title = element_text(face="bold", size=14),axis.text.x = element_text(size = 8, face="bold"))

boxplot_pre_dRep_MAG_ct_vs
 </pre>

# Supplementary Figure 5E

<pre>
  #--------------------------Suppl Figure 5E: 220620 Number of MAGs per dRep Winner bin--------------------------

  #Read in CSV
df_MAGs_per_dRep_Winner<-read.csv('220620_MAGs_per_dRep_Winner.csv',
                             sep=",",
                             header = T)
df_MAGs_per_dRep_Winner$MAGs_within_secondary_cluster<-factor(df_MAGs_per_dRep_Winner$MAGs_within_secondary_cluster, levels=c("1", "2", "3","4","5","6-10","11-20","21-50"))

barplot_MAGs_per_dRep_Winner <- ggplot(data=df_MAGs_per_dRep_Winner, aes(x=MAGs_within_secondary_cluster, y=Number_of_clusters)) +
  geom_bar(stat="identity", color="black", fill="#33CCCC", width=0.7)+
  scale_y_continuous(limits=c(0,1600),expand = c(0,0))+
  labs(y="Number of secondary clusters\n(Total Reconstructed Genomes)",x="Quality-filtered MAGs within secondary cluster")+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size=14),axis.text.x = element_text(size = 12, face="bold"))

barplot_MAGs_per_dRep_Winner
  </pre>

# Supplementary Figure 5F

<pre>
#--------------------------Suppl Figure 5F: 220227 pre vs. post dRep MAG Count--------------------------
library(ggplot2)
library(RColorBrewer)

#Read in csv
df_MAG_count <- read.csv('220227_pre-vs-post_dRep_MAG_count.csv',
                            sep=",",
                            header = T)

#Make lineage barplot
ggplot(data=df_MAG_count, aes(x=pre_dRep,y=post_dRep))+
  geom_point()+
  geom_smooth(span=1)+
  scale_y_continuous(limits=c(0,250),expand = c(0,0))+
  theme_bw()+
  ggtitle("Reduction in MAG Count Following Dereplication") +
  xlab("Quality-filtered Putative Genomes")+
  ylab("Reconstructed Genomes")+
  theme(plot.title = element_text(face="bold",hjust = 0.5))+
  theme(axis.title.x = element_text(face="bold"),axis.title.y = element_text(face="bold"))
</pre>

# Supplementary Figure 6A

<pre>
#--------------------------Suppl Figure 6A: 220228 Assembly Type Percent--------------------------
#Read in CSV
df_Assembly_Type<-read.csv('220228_Assembly_Type_Percent.csv',
                              sep=",",
                              header = T)
df_Assembly_Type$Group<-factor(df_Assembly_Type$Group, levels=c("All MAGs", "Quality-Filtered (HQ or MQ)", "dRep Winners"))

barplot_Assembly_Type <- ggplot(data=df_Assembly_Type, aes(x=Group, y=Percent, fill=Assembly_Type)) +
  geom_bar(stat="identity", color="black", width=0.7)+
  scale_y_continuous(limits=c(0,101),expand = c(0,0))+
  scale_x_discrete(labels = c("All Putative Genomes","QF Putative Genomes","Reconstructed Genomes"))+
  labs(y="Percent",fill="Assembly Type")+
  theme_classic()+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=14),axis.text.x = element_text(size = 12, face="bold",angle=315, hjust=0))

barplot_Assembly_Type

</pre>

# Supplementary Figure 6B

<pre>
#--------------------------Suppl Figure 6B: 220228 Assembly Source Percent--------------------------
#Read in CSV
df_Assembly_Source<-read.csv('220228_Assembly_Source_Percent.csv',
                           sep=",",
                           header = T)
df_Assembly_Source$Group<-factor(df_Assembly_Source$Group, levels=c("All MAGs", "Quality-Filtered (HQ or MQ)", "dRep Winners"))

barplot_Assembly_Source <- ggplot(data=df_Assembly_Source, aes(x=Group, y=Percent, fill=Assembly_Source)) +
  geom_bar(stat="identity", color="black", width=0.7)+
  scale_y_continuous(limits=c(0,101),expand = c(0,0))+
  scale_x_discrete(labels = c("All Putative Genomes","QF Putative Genomes","Reconstructed Genomes"))+
  scale_fill_manual(values=c("#CCFFCC", "#CC9966"))+
  labs(y="Percent",fill="Assembly Source")+
  theme_classic()+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=14),axis.text.x = element_text(size = 12, face="bold",angle=315, hjust=0))

barplot_Assembly_Source

</pre>

# Supplementary Figure 7A

<pre>
#------------------------Suppl Figure 7A: Persisting MAGs: Number of unique, persisting MAGs of a given taxa------------------------

library(ggplot2)
library(ggprism)

df_persistingMAGs_per_taxa <- data.frame(samples  = c("1","2","3","4","5","6","7","8","9","10","11","12","13-14","15-16","17-18","19-20","21-22","23-25",""," ","38"),
                                                       count = c("73","39","19","15","12","10","4","7","10","6","5","4","4","5","4","3","0","3","0","0","1")
)

df_persistingMAGs_per_taxa$samples <- factor(df_persistingMAGs_per_taxa$samples, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13-14","15-16","17-18","19-20","21-22","23-25",""," ","38"))
df_persistingMAGs_per_taxa$count <- as.numeric(df_persistingMAGs_per_taxa$count)

ggplot(data=df_persistingMAGs_per_taxa, aes(x=samples, y=count)) +
  geom_bar(stat="identity",width=0.75, fill="#CCCCFF", color="black", size=0.65)+
  theme_bw()+
  scale_y_continuous(limits=c(0,74), expand=c(0,0), guide = guide_prism_minor())+
  xlab("Persisting MAG count")+
  ylab("Taxa Count")+
  ggtitle("Number of unique, persisting MAGs of a given taxon across all individuals")+
  annotate("text",x=17.6,y=70,label="n=224 total unique taxa identifiers")+
  theme(axis.text.x = element_text(size=9.5, face="bold", angle = 270, vjust=0.5, hjust=0), axis.text.y = element_text(size=11, face="bold"), plot.title = element_text(size=13, face="bold", hjust=0.5), axis.title = element_text(size=11.5, face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

</pre>

# Supplementary Figure 7B

<pre>
#------------------------Suppl Figure 7B: Persisting MAGs: Samples within individual------------------------

library(ggplot2)
library(ggprism)

df_persistingMAGs_samples_per_individual <- data.frame(samples  = c("1","2", "3", "4","5","6","7","8","9","10-15","16-21"),
                                                       count = c("0","549","206","154","66","33","24","12","13","24","12")
)

df_persistingMAGs_samples_per_individual$samples <- factor(df_persistingMAGs_samples_per_individual$samples, levels = c("1","2", "3", "4","5","6","7","8","9","10-15","16-21"))
df_persistingMAGs_samples_per_individual$count <- as.numeric(df_persistingMAGs_samples_per_individual$count)

ggplot(data=df_persistingMAGs_samples_per_individual, aes(x=samples, y=count)) +
  geom_bar(stat="identity",width=0.75, fill="#CCCCFF", color="black", size=0.65)+
  theme_bw()+
  scale_y_continuous(limits=c(0,560), expand=c(0,0), guide = guide_prism_minor())+
  xlab("Sample Count")+
  ylab("Persisting MAG Count")+
  annotate("text",x=9.6,y=530,label="n=1,093 total persisting MAGs")+
  ggtitle("Number of samples within an individual a persisting MAG is present in")+
  theme(axis.text.x = element_text(size=9.5, face="bold"), axis.text.y = element_text(size=11, face="bold"), plot.title = element_text(size=13, face="bold", hjust=0.5), axis.title = element_text(size=11.5, face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

</pre>

# Supplementary Figure 7C

<pre>
#------------------------Suppl Figure 7C: Persisting MAGs: Number of samples across individuals that carry a persisting MAG of a taxa------------------------

library(ggplot2)
library(ggprism)

df_samples_w_persistingMAGs_of_a_taxa <- data.frame(samples  = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21-30","31-40","41-50","51-60","61-80","81-100","101-119"),
                                         count = c("46","12","22","11","11","9","7","7","2","7","7","3","4","3","3","5","4","1","0","20","13","9","6","6","4","2")
)

df_samples_w_persistingMAGs_of_a_taxa$samples <- factor(df_samples_w_persistingMAGs_of_a_taxa$samples, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21-30","31-40","41-50","51-60","61-80","81-100","101-119"))
df_samples_w_persistingMAGs_of_a_taxa$count <- as.numeric(df_samples_w_persistingMAGs_of_a_taxa$count)

ggplot(data=df_samples_w_persistingMAGs_of_a_taxa, aes(x=samples, y=count)) +
  geom_bar(stat="identity",width=0.75, fill="#CCCCFF", color="black", size=0.65)+
  theme_bw()+
  scale_y_continuous(limits=c(0,50), expand=c(0,0), guide = guide_prism_minor())+
  xlab("Sample Count")+
  ylab("Taxa Count")+
  annotate("text",x=21.6,y=47.5,label="n=204 samples across individuals")+
  ggtitle("Number of samples across individuals that carry a persisting MAG of a taxon")+
  theme(axis.text.x = element_text(size=9.5, face="bold", angle = 270, vjust=0.5, hjust=0), axis.text.y = element_text(size=11, face="bold"), plot.title = element_text(size=13, face="bold", hjust=0.5), axis.title = element_text(size=11.5, face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

</pre>
