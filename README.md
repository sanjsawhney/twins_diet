# twins_diet
Code used for data visualization and analysis in Sawhney et al (Weaning accelerates and transforms within-host adaptation in the infant gut microbiome).

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
 # scale_y_continuous(limits=c(30,100),expand = c(0,3))+
  scale_y_continuous(limits=c(30,107), expand = c(0,0), breaks=c(40,60,80,100))+
  scale_fill_manual(values=c("#3399FF", "#33FF99", "#FF9900"))+
  theme(axis.title.y=element_blank())+
  #Remove fill= if not showing grouped boxplot
  #labs(y="Completeness (%)", fill="Assembly Type")+
  labs(y="Completeness (%)")+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, hide.ns = FALSE, label = "p.signif", test = "kruskal.test", tip.length = 0.005, label.y = c(98,99.5,101.5), vjust=.75, size =5)+
  theme(legend.position='none',panel.border = element_rect(color = "black", size=1),axis.ticks = element_line(colour = "black", size = 0.7), axis.title = element_blank(),axis.text = element_text(size = 24, face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
 # theme(legend.position='none', axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=24),axis.text.x = element_text(size = 16, face="bold"),axis.text.y = element_text(face="bold",size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

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
#  scale_y_continuous(limits=c(0,10),expand = c(0,0.2))+
  scale_y_continuous(limits=c(0,12.55), expand = c(0,0.2), breaks=c(0.0,2.5,5.0,7.5,10.0,12.5))+
  scale_fill_manual(values=c("#3399FF", "#33FF99", "#FF9900"))+
  theme(axis.title.y=element_blank())+
  #Remove fill= if not showing grouped boxplot
  #labs(y="Contamination (%)", fill="Assembly Type")+
  labs(y="Contamination (%)")+
  theme_bw()+
  stat_compare_means(comparisons = my_comparisons, hide.ns = TRUE, label = "p.signif", test = "kruskal.test", tip.length = 0.005, label.y = c(11,11.4,11.8), vjust=.75, size =5)+
  theme(legend.position='none',panel.border = element_rect(color = "black", size=1),axis.ticks = element_line(colour = "black", size = 0.7), axis.title = element_blank(),axis.text = element_text(size = 24, face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  #theme(legend.position='none', axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=24),axis.text.x = element_text(size = 16, face="bold"),axis.text.y = element_text(face="bold",size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

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
#  theme(legend.position='none', axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=24),axis.text.x = element_text(size = 16, face="bold"),axis.text.y = element_text(face="bold",size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

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
#  theme(legend.position='none', axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=24),axis.text.x = element_text(size = 16, face="bold"),axis.text.y = element_text(face="bold",size = 16),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),plot.margin = margin(10,3,3,5))

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
 # theme(legend.position='none', axis.title.x = element_blank(),axis.title.y = element_text(face="bold", size=24),axis.text.x = element_text(size = 16, face="bold"),axis.text.y = element_text(size = 16,face="bold"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

boxplot_dRep_MAG_contig_count

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

# Supplementary Figure 5E

<pre>
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
