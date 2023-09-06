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

# Figure 4C

<pre>
#------------------------Figure 4C1: Most abundant infant MAGs: Persisting vs. Transient (COUNT ONLY BARPLOT)------------------------
library(ggplot2)
library(ggprism)

#Read in dataframe that ONLY has COUNT info
df_most_abundant_infant_MAGs_count<-read.csv('230118_MostAbundantInfantMAGs_count_appearance.csv',
                                       sep=",",
                                       header = T)

#Order by phylogeny
df_most_abundant_infant_MAGs_count$Taxa <- factor(df_most_abundant_infant_MAGs_count$Taxa, levels=c("Bifidobacterium pseudocatenulatum","Bifidobacterium bifidum","Bifidobacterium longum","Bifidobacterium breve","Escherichia coli","Akkermansia muciniphila","Parabacteroides distasonis","Phocaeicola vulgatus","Bacteroides uniformis","Bacteroides fragilis","Enterococcus faecalis","Streptococcus salivarius","Faecalibacillus intestinalis","Erysipelatoclostridium ramosum","Clostridium innocuum","Veillonella parvula","Intestinibacter bartlettii","Flavonifractor plautii","Ruminococcus bromii","Ruminococcus bicirculans","Faecalibacterium prausnitzii","Gemmiger formicilis","Coprococcus eutactus","Anaerostipes hadrus","Anaerostipes caccae","Blautia wexlerae","Blautia massiliensis","Sellimonas intestinalis","Ruminococcus gnavus","Ruminococcus torques","Agathobacter rectalis"))

#Plot ONLY the COUNT data
ggplot(data=df_most_abundant_infant_MAGs_count, aes(x=Taxa, y=Count)) +
  geom_bar(stat="identity",width=0.5, color="black", size=0.65)+
  ylab("RG Count")+
  ggtitle("Taxa with greatest Reconstructed Genome diversity within infants")+
  scale_y_continuous(limits=c(0,80), expand=c(0,0), guide = guide_prism_minor())+
  theme_bw()+
  theme(axis.text.x = element_text(size=9, face="bold", angle = 90, vjust=0.5, hjust=1), axis.text.y = element_text(size=13, face="bold"), plot.title = element_text(size=15, face="bold", hjust=0.5), axis.title.y = element_text(size=14, face="bold"), axis.title.x=element_blank(), legend.title=element_blank(), legend.text=element_text(face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))


#------------------------Figure 4C2: Most abundant infant MAGs: Persisting vs. Transient (PERCENTAGE ONLY BARPLOT)------------------------
#Read in dataframe that ONLY has PERCENTAGE info
df_most_abundant_infant_MAGs_percentage<-read.csv('230118_MostAbundantInfantMAGs_percentage.csv',
                                                  sep=",",
                                                  header = T)

#Order by phylogeny
df_most_abundant_infant_MAGs_percentage$Taxa <- factor(df_most_abundant_infant_MAGs_percentage$Taxa, levels=c("Bifidobacterium pseudocatenulatum","Bifidobacterium bifidum","Bifidobacterium longum","Bifidobacterium breve","Escherichia coli","Akkermansia muciniphila","Parabacteroides distasonis","Phocaeicola vulgatus","Bacteroides uniformis","Bacteroides fragilis","Enterococcus faecalis","Streptococcus salivarius","Faecalibacillus intestinalis","Erysipelatoclostridium ramosum","Clostridium innocuum","Veillonella parvula","Intestinibacter bartlettii","Flavonifractor plautii","Ruminococcus bromii","Ruminococcus bicirculans","Faecalibacterium prausnitzii","Gemmiger formicilis","Coprococcus eutactus","Anaerostipes hadrus","Anaerostipes caccae","Blautia wexlerae","Blautia massiliensis","Sellimonas intestinalis","Ruminococcus gnavus","Ruminococcus torques","Agathobacter rectalis"))

#Plot ONLY the PERCENTAGE data
ggplot(data=df_most_abundant_infant_MAGs_percentage, aes(x=Taxa, y=Percent, fill=Appearance)) +
  geom_bar(stat="identity",width=0.75, color="black", size=0.65)+
  scale_fill_manual(values=c('#634E52','#BD8F3F'))+
  ylab("Percentage")+
  ggtitle("Taxa with greatest MAG diversity within infants")+
  scale_y_continuous(limits=c(0,100.1), expand=c(0,0), guide = guide_prism_minor())+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, face="bold", angle = 270, vjust=0.5, hjust=0), axis.text.y = element_text(size=11, face="bold"), plot.title = element_text(size=13, face="bold", hjust=0.5), axis.title.y = element_text(size=11.5, face="bold"), axis.title.x=element_blank(), legend.title=element_blank(), legend.text=element_text(face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

</pre>

# Figure 4D

<pre>
#------------------------Figure 4D: Top Taxa Avg first/last apps: All MAGs 90% BOXPLOT------------------------

#Create dataframe of average first and last appearance, and 2 x (first +/- CI, last +/- CI)
#Read in dataframe for **All MAGs** 90% BOXPLOT
df_taxa_first_last_appearance_AllMAGs_90CI_bounds<-read.csv('230110_Top10Taxa_AvgFirstLastAppearance_AllMAGs_90CIbounds.csv',
                                                           sep=",",
                                                           header = T)

#Order df by phylogeny and appearance
df_taxa_first_last_appearance_AllMAGs_90CI_bounds$Taxa <- factor(df_taxa_first_last_appearance_AllMAGs_90CI_bounds$Taxa, levels = c("Bacteroides uniformis", "Bacteroides fragilis", "Phocaeicola vulgatus", "Parabacteroides distasonis", "Bifidobacterium pseudocatenulatum", "Bifidobacterium bifidum", "Bifidobacterium longum", "Faecalibacterium prausnitzii", "Anaerostipes hadrus", "Blautia wexlerae"))
df_taxa_first_last_appearance_AllMAGs_90CI_bounds$Appearance <- factor(df_taxa_first_last_appearance_AllMAGs_90CI_bounds$Appearance, levels = c("Last","First"))


#Create dataframe of taxa, ymin, ymax
#Read in dataframe for **All MAGs** BOXPLOT
df_taxa_first_last_appearance_AllMAGs_90CI_segment<-read.csv('230110_Top10Taxa_AvgFirstLastAppearance_AllMAGs_90CIsegment.csv',
                                                            sep=",",
                                                            header = T)
#Order df by phylogeny
df_taxa_first_last_appearance_AllMAGs_90CI_segment$Taxa <- factor(df_taxa_first_last_appearance_AllMAGs_90CI_segment$Taxa, levels = c("Bacteroides uniformis", "Bacteroides fragilis", "Phocaeicola vulgatus", "Parabacteroides distasonis", "Bifidobacterium pseudocatenulatum", "Bifidobacterium bifidum", "Bifidobacterium longum", "Faecalibacterium prausnitzii", "Anaerostipes hadrus", "Blautia wexlerae"))


# **All MAGs** w/ 90% CI
ggplot(data=df_taxa_first_last_appearance_AllMAGs_90CI_bounds, aes(x=Taxa, y=MOL))+
  geom_boxplot(data=df_taxa_first_last_appearance_AllMAGs_90CI_bounds, aes(x=Taxa, y=MOL, fill=Appearance), position=position_dodge(0)) +
  scale_fill_manual(values=c(First="#33B31A80",Last="#B3331A80"))+
  geom_segment(data=df_taxa_first_last_appearance_AllMAGs_90CI_segment, aes(x=Taxa, xend=Taxa, y=First_Appearance_Upper90CI, yend=Last_Appearance_Lower90CI))+
  coord_flip()+
  ggtitle("All MAGs per Individual (90% CI)")+
  ylab("MOL")+
  scale_y_continuous(limits=c(0,105),expand = c(0,0))+
  theme_bw()+
  theme(axis.text.x = element_text(size=10, face="bold"), axis.text.y = element_text(size=8, face="bold"), plot.title=element_text(size=16, face="bold"), axis.title.y=element_blank(), axis.title.x=element_text(size=14, face="bold"), legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=0.75))

</pre>

# Figure 4E

<pre>

#Read dataframe
df_strain_cooccurrence_rate<-read.csv('230402_Strain_Cooccurrence_rate_multi_only.csv', sep=",", header = T)

#Order by phylogeny
df_strain_cooccurrence_rate$Taxa <- factor(df_strain_cooccurrence_rate$Taxa, levels=c("Bifidobacterium pseudocatenulatum","Bifidobacterium bifidum","Bifidobacterium longum","Bifidobacterium breve","Escherichia coli","Akkermansia muciniphila","Parabacteroides distasonis","Phocaeicola vulgatus","Bacteroides uniformis","Bacteroides fragilis","Enterococcus faecalis","Streptococcus salivarius","Faecalibacillus intestinalis","Erysipelatoclostridium ramosum","Clostridium innocuum","Veillonella parvula","Intestinibacter bartlettii","Flavonifractor plautii","Ruminococcus bromii","Ruminococcus bicirculans","Faecalibacterium prausnitzii","Gemmiger formicilis","Coprococcus eutactus","Anaerostipes hadrus","Anaerostipes caccae","Blautia wexlerae","Blautia massiliensis","Sellimonas intestinalis","Ruminococcus gnavus","Ruminococcus torques","Agathobacter rectalis"))

#Plot
ggplot(data=df_strain_cooccurrence_rate, aes(x=Taxa, y=Percentage)) +
  geom_bar(stat="identity",width=0.75, color="black", size=0.65)+
  ylab("Percentage of TPs\nwith multiple strains")+
  ggtitle("Strain Co-occurrence Rate")+
  scale_y_continuous(limits=c(0,42.5), expand=c(0,0), guide = guide_prism_minor())+
  theme_bw()+
  theme(axis.text.x = element_text(size=10, face="bold", angle = 90, vjust=0.5, hjust=1), axis.text.y = element_text(size=11, face="bold"), plot.title = element_text(size=16, face="bold", hjust=0.5), axis.title.y = element_text(size=12, face="bold"), axis.title.x=element_blank(), legend.title=element_blank(), legend.text=element_text(face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

</pre>

# Figure 4F

<pre>

#Read in dataframe that has Taxa count and MOL_range for EACH stool sample
df_straindiversity_per_MOLrange_taxa_20individuals_nozeros<-read.csv('230402_StrainDiversity_TaxaCt_perSample_20individuals_ISprof_IScomp_dRep.csv',sep=",",header = T)
  
df_straindiversity_per_MOLrange_taxa_20individuals_nozeros$MOL <- as.numeric(df_straindiversity_per_MOLrange_taxa_20individuals_nozeros$MOL)

#Plot **MOL** x count AS DOT PLOT - 20 individuals
ggplot(data=df_straindiversity_per_MOLrange_taxa_20individuals_nozeros, aes(x=MOL, y=Count)) +
  geom_point(size=2, position=position_jitter(h=0,w=0.75, seed=6), alpha=0.4)+
  #geom_point(size=1, position=position_jitter(h=0,w=0.5, seed=3),alpha=0.5)+
  geom_smooth(method = 'loess', color="black", fill="#996666")+
  scale_x_continuous(limits=c(-0.1,100), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(-0.1,9.25), expand=c(0,0), breaks=c(2,4,6,8))+
  theme_bw()+
  xlab("MOL")+
  ylab("Taxa Count")+
  ggtitle("Number of taxa with strain diversity per sample")+
  theme(axis.text = element_text(size=14, face="bold"), plot.title = element_text(size=17, face="bold", hjust=0.5), axis.title = element_text(size=15, face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))+
  theme(plot.margin = margin(0.2,0.5,0.1,0.25, "cm"))

</pre>

# Figure 5A

<pre>
#------------------------Figure 5A1: Shared-vs Not shared------------------------
library(ggplot2)
library(ggprism)

df_familyRGs<-read.csv('230413_FamilyRGs.csv',
                              sep=",",
                              header = T)


df_familyRGs$Group <- factor(df_familyRGs$Group, levels = c("Single individual","Multiple individuals"))


ggplot(data=df_familyRGs, aes(x=Comparison, y=Count, fill=Group)) +
  geom_bar(stat="identity",width=0.45, color="black")+
  scale_fill_manual(values=c("#DDDDDD","#86A3BF"))+
  scale_y_continuous(limits=c(0,2600), expand=c(0,0), guide = guide_prism_minor())+
  xlab("RGs (Family)")+
  ylab("Count")+
  theme_classic()+
  theme(axis.text = element_text(size=11, face="bold"), axis.title.y = element_text(size=13, face="bold"), axis.title.x = element_blank(), legend.title = element_blank(), legend.text = element_blank(), legend.position = "bottom")

#------------------------Figure 5A2: Sharing type------------------------

df_shared_familyRGs<-read.csv('230413_Shared_FamilyRGs.csv',
                       sep=",",
                       header = T)


df_shared_familyRGs$Group <- factor(df_shared_familyRGs$Group, levels = c("Family Triad","Dyad: Mother-Infant","Dyad: Infant-Infant"))


ggplot(data=df_shared_familyRGs, aes(x=Comparison, y=Count, fill=Group)) +
  geom_bar(stat="identity",width=0.45, color="black")+
  scale_fill_manual(name="Sharing type", values=c('#CB7122','#4E9885','#2B659D'))+
  scale_y_continuous(limits=c(0,800), expand=c(0,0), guide = guide_prism_minor())+
  xlab("Shared RGs")+
  theme_classic()+
  theme(axis.text = element_text(size=11, face="bold"), axis.title.y = element_text(size=13, face="bold"), axis.title.x = element_blank(), legend.title = element_blank(), legend.text = element_blank(), legend.position = "bottom")

</pre>

# Figure 5B

<pre>
#------------------------Figure 5B: Strain Sharing: Infant-Infant by Bin - INCIDENCE RATE------------------------


##PLOT INCIDENCE RATE


#Read in dataframe that ONLY has PERCENTAGE info
df_infant_infant_bins<-read.csv('230413_StrainSharing_Infant_Infant_Bins_Incidence.csv',
                                sep=",",
                                header = T)

df_infant_infant_bins$Bin <- as.factor(df_infant_infant_bins$Bin)
df_infant_infant_bins$MOL <- as.factor(df_infant_infant_bins$MOL)

#Order
df_infant_infant_bins$Taxa <- factor(df_infant_infant_bins$Taxa, levels=c("Bacteroidales","Bifidobacteriaceae","Enterococcaceae","Erysipelotrichales","Lachnospiraceae","Oscillospiraceae","Other Eubacteriales","Other"))
df_infant_infant_bins$MOL_elapsed <- factor(df_infant_infant_bins$MOL_elapsed, levels=c("0-6","6-9","9-12","12-15","15-36","36-94"))

#Plot ONLY the PERCENTAGE data
ggplot(data=df_infant_infant_bins, aes(x=MOL_elapsed, y=Incidence_per_individual, fill=Taxa)) +
  geom_bar(stat="identity",width=0.5, color="black", size=0.5)+
  scale_fill_manual(values=c('#ffffbf','#d53e4f','#7DC0A7','#ffffff','#e0e0e0','#878787','#4B86B8','#333333'))+
  xlab("MOL")+
  ylab("Incidence rate")+
  geom_vline(xintercept=3.5, size=0.75, color=c('#0066FF'), linetype="dashed")+
  ggtitle("New Strain Sharing Events per Month Between Twins")+
  scale_y_continuous(limits=c(0,1.07), expand=c(0,0), guide = guide_prism_minor())+
  theme_bw()+
  theme(axis.text = element_text(size=21, face="bold"), plot.title = element_blank(), axis.title = element_text(size=26, face="bold"), legend.title=element_blank(), legend.text=element_text(face="bold"), legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1.25))

</pre>

# Figure 5C

<pre>
#------------------------Figure 5C: Strain Sharing: Dyad/Triad Taxa representation------------------------

library(ggplot2)
library(ggprism)

###STACKED BARPLOT - COMBINED
#Read in dataframe with count for transient and persisting MAGs
df_dyad_triad_taxa<-read.csv('230131_Dyad_Triad_taxa.csv',
                                                  sep=",",
                                                  header = T)

#Order by phylogeny
df_dyad_triad_taxa$Taxa <- factor(df_dyad_triad_taxa$Taxa, levels=c("Bifidobacterium pseudocatenulatum","Bifidobacterium adolescentis","Bifidobacterium bifidum","Bifidobacterium longum","Bifidobacterium breve","Escherichia coli","Alistipes putredinis","Alistipes finegoldii","Parabacteroides distasonis","Phocaeicola vulgatus","Bacteroides uniformis","Bacteroides thetaiotaomicron","Bacteroides fragilis","Streptococcus thermophilus","Enterococcus faecalis","Faecalibacillus intestinalis","Erysipelatoclostridium ramosum","Clostridium innocuum","Veillonella parvula","Ruminococcus bromii","Ruminococcus bicirculans","Faecalibacterium prausnitzii","Anaerostipes hadrus","Anaerostipes caccae","Lacrimispora celerecrescens","Blautia wexlerae","Blautia massiliensis","Sellimonas intestinalis","Ruminococcus gnavus","Ruminococcus torques","Ruminococcus faecis"))

#Create STACKED BARPLOT - COMBINED
ggplot(data=df_dyad_triad_taxa, aes(x=Taxa, y=Count, fill=Type)) +
  geom_bar(stat="identity",width=0.75, color="black", size=0.65)+
  scale_fill_manual(name="Sharing type", values=c('#CB7122','#4E9885','#2B659D'))+
  ylab("Shared MAG Count")+
  ggtitle("Most common taxa representing intra-family strain sharing")+
  scale_y_continuous(limits=c(0,22.575), expand=c(0,0), guide = guide_prism_minor())+
  theme_bw()+
  theme(axis.text.x = element_text(size=8, face="bold", angle = 270, vjust=0.5, hjust=0), axis.text.y = element_text(size=11, face="bold"), plot.title = element_text(size=13, face="bold", hjust=0.5), axis.title.y = element_text(size=11.5, face="bold"), axis.title.x=element_blank(), legend.title=element_text(face="bold"), legend.text=element_text(face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

</pre>

# Figure 5D

<pre>
#------------------------Figure 5D1: Strain Sharing: Mother-Infant First Appearance: Bacteroidales------------------------

library(ggplot2)
library(ggprism)

#Read in dataframe with first appearance of a shared strain in mother and first infant
df_mother_infant_first_app_bacteroidales<-read.csv('230207_SharedMAG_MotherInfant_FirstAppearance_Bacteroidales.csv',
                                                   sep=",",
                                                   header = T)

#Calculate group means
library(plyr)
mean_mother_infant_first_app_bacteroidales <- ddply(df_mother_infant_first_app_bacteroidales, "Group", summarise, grp.mean=mean(MOL))
median_mother_infant_first_app_bacteroidales <- ddply(df_mother_infant_first_app_bacteroidales, "Group", summarise, grp.median=median(MOL))

#Create density plot
ggplot(data=df_mother_infant_first_app_bacteroidales, aes(x=MOL, fill=Group)) +
  geom_density(alpha=0.4)+
  scale_color_manual(values=c("#FFCC66", "black"))+
  scale_fill_manual(values=c("#F9DA78", "black"))+
  scale_x_continuous(limits=c(0,97), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(0,0.1500001), expand=c(0,0), guide = guide_prism_minor())+
  geom_vline(data=mean_mother_infant_first_app_bacteroidales, aes(xintercept=grp.mean, color=Group), size=0.75, linetype="dashed")+
  annotate("text",x=75,y=0.142,size=5,label="n=60 shared strains")+
  ggtitle("Bacteroidales")+
  ylab("Density")+
  theme_bw()+
  theme(axis.text = element_text(size=18, face="bold"), plot.title = element_text(size=21, face="bold", hjust=0.5), axis.title = element_text(size=20, face="bold"), legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

#------------------------Figure 5D2: Strain Sharing: Mother-Infant First Appearance: Lachnospiraceae------------------------

library(ggplot2)
library(ggprism)

#Read in dataframe with first appearance of a shared strain in mother and first infant
df_mother_infant_first_app_lachnospiraceae<-read.csv('230207_SharedMAG_MotherInfant_FirstAppearance_Lachnospiraceae.csv',
                                                   sep=",",
                                                   header = T)

#Calculate group means
library(plyr)
mean_mother_infant_first_app_lachnospiraceae <- ddply(df_mother_infant_first_app_lachnospiraceae, "Group", summarise, grp.mean=mean(MOL))
median_mother_infant_first_app_lachnospiraceae <- ddply(df_mother_infant_first_app_lachnospiraceae, "Group", summarise, grp.median=median(MOL))

#Create density plot
ggplot(data=df_mother_infant_first_app_lachnospiraceae, aes(x=MOL, fill=Group)) +
  geom_density(alpha=0.4)+
  scale_color_manual(values=c("#FFCC66", "black"))+
  scale_fill_manual(values=c("#F9DA78", "black"))+
  scale_x_continuous(limits=c(0,97), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(0,0.1500001), expand=c(0,0), guide = guide_prism_minor())+
  geom_vline(data=mean_mother_infant_first_app_lachnospiraceae, aes(xintercept=grp.mean, color=Group), size=0.75, linetype="dashed")+
  annotate("text",x=75,y=0.141,size=5,label="n=51 shared strains")+
  ggtitle("Lachnospiraceae")+
  ylab("Density")+
  theme_bw()+
  theme(axis.text = element_text(size=18, face="bold"), plot.title = element_text(size=21, face="bold", hjust=0.5), axis.title = element_text(size=20, face="bold"), legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

#------------------------Figure 5D3: Strain Sharing: Mother-Infant First Appearance: Oscillospiraceae------------------------

library(ggplot2)
library(ggprism)

#Read in dataframe with first appearance of a shared strain in mother and first infant
df_mother_infant_first_app_oscillospiraceae<-read.csv('230207_SharedMAG_MotherInfant_FirstAppearance_Oscillospiraceae.csv',
                                                      sep=",",
                                                      header = T)

#Calculate group means
library(plyr)
mean_mother_infant_first_app_oscillospiraceae <- ddply(df_mother_infant_first_app_oscillospiraceae, "Group", summarise, grp.mean=mean(MOL))
median_mother_infant_first_app_oscillospiraceae <- ddply(df_mother_infant_first_app_oscillospiraceae, "Group", summarise, grp.median=median(MOL))

#Create density plot
ggplot(data=df_mother_infant_first_app_oscillospiraceae, aes(x=MOL, fill=Group)) +
  geom_density(alpha=0.4)+
  scale_color_manual(values=c("#FFCC66", "black"))+
  scale_fill_manual(values=c("#F9DA78", "black"))+
  scale_x_continuous(limits=c(0,97), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(0,0.1500001), expand=c(0,0), guide = guide_prism_minor())+
  geom_vline(data=mean_mother_infant_first_app_oscillospiraceae, aes(xintercept=grp.mean, color=Group), size=0.75, linetype="dashed")+
  annotate("text",x=75,y=0.141,size=5,label="n=21 shared strains")+
  ggtitle("Oscillospiraceae")+
  ylab("Density")+
  theme_bw()+
  theme(axis.text = element_text(size=18, face="bold"), plot.title = element_text(size=21, face="bold", hjust=0.5), axis.title = element_text(size=20, face="bold"), legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

#------------------------Figure 5D4: Strain Sharing: Mother-Infant First Appearance: Bifidobacteriaceae------------------------

library(ggplot2)
library(ggprism)

#Read in dataframe with first appearance of a shared strain in mother and first infant
df_mother_infant_first_app_bifidobacteriaceae<-read.csv('230207_SharedMAG_MotherInfant_FirstAppearance_Bifidobacteriaceae.csv',
                                                        sep=",",
                                                        header = T)

#Calculate group means
library(plyr)
mean_mother_infant_first_app_bifidobacteriaceae <- ddply(df_mother_infant_first_app_bifidobacteriaceae, "Group", summarise, grp.mean=mean(MOL))
median_mother_infant_first_app_bifidobacteriaceae <- ddply(df_mother_infant_first_app_bifidobacteriaceae, "Group", summarise, grp.median=median(MOL))

#Create density plot
ggplot(data=df_mother_infant_first_app_bifidobacteriaceae, aes(x=MOL, fill=Group)) +
  geom_density(alpha=0.4)+
  scale_color_manual(values=c("#FFCC66", "black"))+
  scale_fill_manual(values=c("#F9DA78", "black"))+
  scale_x_continuous(limits=c(0,97), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(0,0.1500001), expand=c(0,0), guide = guide_prism_minor())+
  geom_vline(data=mean_mother_infant_first_app_bifidobacteriaceae, aes(xintercept=grp.mean, color=Group), size=0.75, linetype="dashed")+
  annotate("text",x=75,y=0.141,size=5,label="n=20 shared strains")+
  ggtitle("Bifidobacteriaceae")+
  ylab("Density")+
  theme_bw()+
  theme(axis.text = element_text(size=18, face="bold"), plot.title = element_text(size=21, face="bold", hjust=0.5), axis.title = element_text(size=20, face="bold"), legend.position="none", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

#------------------------Figure 5D: Strain Sharing: Paired samples Wilcoxon tests------------------------
bacteroidales_mother <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,6,6,7,23,23,24,36,36,36,6,7,36)
bacteroidales_infant <-c(9,36,6,6,6,6,6,6,6,6,6,6,6,7,7,7,9,10,10,11,11,11,11,12,14,14,14,17,18,20,23,23,23,33,33,35,35,35,37,37,52,66,93,95,97,97,7,46,10,10,37,97,97,52,9,16,22,6,7,35)
wilcox.test(bacteroidales_mother, bacteroidales_infant, paired = TRUE)
wilcox.test(bacteroidales_mother, bacteroidales_infant, paired = TRUE, alternative="less")

lachnospiraceae_mother <-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,6,6,6,6,6,7,7,7,7,24,27,36,36,7,7,7,7,7,7,24,24,24,36,36)
lachnospiraceae_infant <-c(6,6,6,9,10,10,10,10,10,12,13,14,16,22,23,35,36,47,66,76,95,12,12,12,16,21,24,13,13,17,24,36,10,12,13,17,7,6,23,12,7,7,7,7,6,6,24,22,23,35,34)
wilcox.test(lachnospiraceae_mother, lachnospiraceae_infant, paired = TRUE)
wilcox.test(lachnospiraceae_mother, lachnospiraceae_infant, paired = TRUE, alternative="less")

oscillospiraceae_mother <-c(0,0,0,0,0,0,0,1,1,6,7,7,36,24,24,36,36,0,7,24,36)
oscillospiraceae_infant <-c(6,10,20,33,36,36,47,36,73,14,12,13,24,22,23,34,35,10,37,24,35)
wilcox.test(oscillospiraceae_mother, oscillospiraceae_infant, paired = TRUE)
wilcox.test(oscillospiraceae_mother, oscillospiraceae_infant, paired = TRUE, alternative="less")

bifidobacteriaceae_mother <-c(0,0,0,0,0,0,0,0,1,24,36,24,24,36,1,7,7,7,7,36)
bifidobacteriaceae_infant <-c(2,2,6,6,6,7,10,10,7,33,41,2,7,24,1,7,7,7,7,35)
wilcox.test(bifidobacteriaceae_mother, bifidobacteriaceae_infant, paired = TRUE)
wilcox.test(bifidobacteriaceae_mother, bifidobacteriaceae_infant, paired = TRUE, alternative="less")

p.adjust(c(2.914e-09, 1.255e-05, 0.001894, 0.1738), method="fdr")
</pre>

# Figure 6A

<pre>

#INFANTS: ALL TAXA: Breadth-adjusted Aggregate popSNP ct by LoPthruTP2 (230311_01.1_Infant_breadth-adjAGGREGATEpopSNP_alltaxa.pdf)
df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa<-read.csv('230314_popSNP_ct_perMAG_byLOPthruTP2_infants_gs.csv', sep=",", header = T)

df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa$Group <- factor(df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa$Group, levels = c("Colonized before WN, timepoint before WN","Colonized before WN, timepoint after WN", "Colonized after WN, timepoint after WN"))
df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa$Initial_Colonization <- factor(df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa$Initial_Colonization, levels = c("Before weaning","After weaning"))

ggplot(data=df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa, aes(x=LoPthruTP2, y=Aggregate_Adjusted_popSNP_Count)) +
  #geom_point(size=1.5, alpha=0.5, fill="black")+
  geom_point(aes(color=Group), size=1, position=position_jitter(h=0.3,w=0.025), alpha=0.5)+
  scale_color_manual(values=c("#33CCFF","#0066FF","#CC9966","#0099FF","#CC9966"))+
  geom_smooth(aes(color=Initial_Colonization), method = 'loess')+
  scale_x_continuous(limits=c(-.085,8), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(-2,135), expand=c(0,0), guide = guide_prism_minor())+
  annotate("text",x=0.77,y=130,label="n = 763 MAGs (194 taxa)")+
  xlab("Length of persistence")+
  ylab("popSNPs since seeding")+
  ggtitle("INFANTS: Breadth- and Length-adjusted popSNP Count by Years since seeding")+
  theme_bw()+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=13,face="bold"), panel.grid.minor=element_blank(), legend.title=element_blank(), legend.position = "none", plot.title=element_text(size=15,face="bold",hjust=0.5), panel.grid.major=element_blank(), panel.border=element_rect(colour="black",size=1))

</pre>

# Figure 6B-C

Made in Prism 9.

# Figure 6D

<pre>
#------------------------Figure 6D: PCoA by ORF COG: BRAY-CURTIS: One point per Category, One entry per ORF: PERCENTAGE - No S------------------------

#load vegan for jaccard distance, ape for pcoa, ggplot2 for visualization
library(ape)
library(vegan)
library(ggplot2)
library(reshape2)
library(ggforce)

#Read in MAG_COG_ct.csv (.Rtab file)
##Change file name
df_category_COG_pct_noS_v1<-read.csv('230321_04_COGpercentage_byIndividualWeaningDiet_noS.csv', sep=",", header = T)
#Convert df to wide
df_category_COG_pct_noS_v1 <- dcast(df_category_COG_pct_noS_v1, Category ~ COG, fill=0)
#Make MAG column the Rownames
df_category_COG_pct_noS <- df_category_COG_pct_noS_v1[,-1]
rownames(df_category_COG_pct_noS) <- df_category_COG_pct_noS_v1[,1]

#Calculate bray distance. Set to matrix
bray_category_COG_pct_noS<-as.matrix(vegdist(df_category_COG_pct_noS, method='bray'))

#Make PCoA (correction is to account for negative eigenvalues) - can use "lingoes" or "cailliez"
#Look through pcoa_accessorygenome_corr on first use to gain an understanding of what $values, $vectors, and $vectors.cor mean
pca_bray_category_COG_pct_noS<-pcoa(bray_category_COG_pct_noS)

#Get PCoA vectors to plot ordination as axis x axis. Set to data frame
pcavectors_bray_category_COG_pct_noS<-as.data.frame(pca_bray_category_COG_pct_noS$vectors)

#Get % variance captured by each axis. Rel_corr_eig = Relative eigenvalues following correction method. Sums to 1.0000
rel_eigen_bray_category_COG_pct_noS<-pca_bray_category_COG_pct_noS$values$Relative_eig
rel_eigen_bray_category_COG_pct_noS

#Add cohort metadata
##Change csv name in both lines
write.csv(pcavectors_bray_category_COG_pct_noS[,c("Axis.1","Axis.2","Axis.3")],"230322_04_bray_Category_COG_pca_pct_noS.csv")
pcavectors_bray_category_COG_pct_noS=read.csv("230322_04_bray_Category_COG_pca_pct_noS.csv", sep=",")

pcavectors_bray_category_COG_pct_noS$Category <- factor(pcavectors_bray_category_COG_pct_noS$Category, levels=c("Breastfed_beforeWeaning","Intermediate_beforeWeaning","Formula_beforeWeaning","Breastfed_afterWeaning","Intermediate_afterWeaning","Formula_afterWeaning","Mother"))
pcavectors_bray_category_COG_pct_noS$Diet <- factor(pcavectors_bray_category_COG_pct_noS$Diet, levels=c("Breastfed","Intermediate","Formula","Mother"))
pcavectors_bray_category_COG_pct_noS$Weaning <- factor(pcavectors_bray_category_COG_pct_noS$Weaning, levels=c("Before","After","Mother"))

#plot
ggplot(
  pcavectors_bray_category_COG_pct_noS,aes(x=Axis.1, y=Axis.2))+
  geom_point(aes(fill=Diet, shape=Weaning), size=5, stroke=1)+
  ggtitle("Distribution of SNP-accruing ORFs across persisting strains by COG")+
  xlab("PCA1 (45.7%)")+ylab("PCA2 (32.8%)")+
  scale_fill_manual(values=c("#FCFDC6", "#C34E6E", "#21134E", "white"))+
  #scale_fill_manual(values=c("#FCFDC6","#F25F88","#5D35DF","#CFD0A3", "#C34E6E", "#21134E", "#ED774D"))+
  scale_shape_manual(values=c(21, 23, 22))+
  scale_color_manual(values=c("#0099FF","#CC9966"))+
  geom_mark_ellipse(aes(color = Weaning, filter = Weaning != 'Mother'), size = 0.75)+
  theme_bw()+theme(aspect.ratio = rel_eigen_bray_category_COG_pct_noS[2]/rel_eigen_bray_category_COG_pct_noS[1])+
  theme(axis.title=element_text(size=14,face="bold"), axis.text=element_blank(), axis.ticks=element_blank(), panel.border=element_rect(colour="black",fill=NA,size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(face="bold"), legend.text=element_text(face="bold"), plot.title = element_text(hjust = 0.25,size=16,face="bold"))+
  guides(fill=guide_legend(override.aes=list(shape=21)))
stat_ellipse()

</pre>

# Figure 6E

<pre>
#------------------------Figure 6E: Boxplot comparing Bray-Curtis distance between Pre-Weaning, Post-Weaning, and Mother mutated gene COG profiles------------------------

#read in boxplot csv
df_pairwise_boxplot<-read.csv('230322_bray_pairwise_boxplot_2groups.csv',
                              sep=",",
                              header = T)

library(ggplot2)
library(ggpubr)
library(ggprism)

df_pairwise_boxplot$Comparison_type___Family <- factor(df_pairwise_boxplot$Comparison_type___Family,levels = c("before___mother___same","before___mother___different","after___mother___same","after___mother___different"))
df_pairwise_boxplot$Family <- factor(df_pairwise_boxplot$Family,levels = c("Same","Different"))
df_pairwise_boxplot$Weaning <- factor(df_pairwise_boxplot$Weaning,levels = c("Pre","Post"))

ggplot(df_pairwise_boxplot, aes(Weaning, Bray_Curtis_distance)) +
  geom_boxplot(aes(colour = Weaning), outlier.shape = NA, lwd=0.75, width=0.55)+
  geom_point(aes(fill = Weaning), size = 5, alpha=0.75, shape=21, position = position_jitterdodge(jitter.width = 0.75, jitter.height = 0, seed=35))+
  scale_x_discrete(labels=c("Pre-Weaning Infant vs. Mother","Post-Weaning Infant vs. Mother"))+
  scale_color_manual(values=c("#0099FF","#CC9966"))+
  scale_fill_manual(values=c("#0099FF","#CC9966"))+
  scale_y_continuous(limits=c(0.18,0.4995), expand=c(0,0), guide = guide_prism_minor())+
  ylab("Bray-Curtis dissimilarity")+
  ggtitle("Distance between COG functional profiles of mutated genes")+
  stat_compare_means(method = "wilcox.test", aes(group = Weaning), comparisons=list(c("Pre", "Post")), label = "p.signif", label.x = 1.5, label.y = 0.45, bracket.size = 0.75, size=10, face="bold")+
  theme_bw()+
  theme(plot.title=element_text(size=17,face="bold", hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_text(size=13, face="bold"), axis.text = element_text(face = "bold",size = 12), legend.title=element_text(face="bold"), legend.text=element_text(face="bold"), panel.border=element_rect(colour="black",fill=NA,size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

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
#------------------------Suppl Figure 7C: Strain Sharing: Samples per shared MAG------------------------

library(ggplot2)
library(ggprism)

df_sharedMAGs_samples_per_MAG <- data.frame(samples  = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15-17","18-20","21-23","24-26","27-29"),
                                                       count = c("204","88","117","81","54","32","29","26","23","13","11","9","9","10","10","7","2","1")
)

df_sharedMAGs_samples_per_MAG$samples <- factor(df_sharedMAGs_samples_per_MAG$samples, levels = c("2","3","4","5","6","7","8","9","10","11","12","13","14","15-17","18-20","21-23","24-26","27-29"))
df_sharedMAGs_samples_per_MAG$count <- as.numeric(df_sharedMAGs_samples_per_MAG$count)

ggplot(data=df_sharedMAGs_samples_per_MAG, aes(x=samples, y=count)) +
  geom_bar(stat="identity",width=0.75, fill="#86A3BF", color="black", size=0.65)+
  theme_bw()+
  scale_y_continuous(limits=c(0,210), expand=c(0,0), guide = guide_prism_minor())+
  xlab("Intra-Family Sample Count")+
  ylab("Shared MAG Count")+
  annotate("text",x=15.8,y=201,label="n=726 total shared MAGs")+
  ggtitle("Appearances per Shared MAG within a Family")+
  theme(axis.text.x = element_text(size=8, face="bold"), axis.text.y = element_text(size=11, face="bold"), plot.title = element_text(size=13, face="bold", hjust=0.5), axis.title = element_text(size=11.5, face="bold"), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black", size=1))

</pre>

# Supplementary Figure 8

<pre>

#------------------------Suppl Figure 8: Aggregate popSNP by Length of persistence (Infant): FACETED BY TAXA, COLORED BY INDIVIDUAL------------------------

library(ggplot2)
library(ggprism)
library(ggpubr)
#library(ggExtra)

#TAXA5CT: Read in INFANT dataframe that has breadth-adjusted popSNP ct by LoPthruTP2
df_popSNP_ct_perMAG_byLOPthruTP2_infants_5ct<-read.csv('230311_popSNP_ct_perMAG_byLOPthruTP2_Taxa5ct_infants_diet.csv',
                                                   sep=",",
                                                   header = T)
  
df_popSNP_ct_perMAG_byLOPthruTP2_infants_5ct$Taxa <- factor(df_popSNP_ct_perMAG_byLOPthruTP2_infants_5ct$Taxa, levels = c("Faecalibacterium prausnitzii","Clostridium sp.","Bifidobacterium pseudocatenulatum","Oscillospiraceae bacterium","Bifidobacterium longum","Collinsella sp.","Bifidobacterium bifidum","Sellimonas intestinalis","[Ruminococcus] gnavus","Acutalibacteraceae bacterium","Enterococcus faecalis","Erysipelatoclostridium ramosum","Parabacteroides distasonis","Bacteroides uniformis","Blautia wexlerae","Bifidobacterium breve","Blautia sp.","[Clostridium] innocuum","Phocaeicola vulgatus","Bacteroides fragilis","Anaerostipes hadrus","Anaerostipes caccae","Veillonella parvula","Bacteroides thetaiotaomicron","Flavonifractor plautii","Blautia massiliensis","Alistipes finegoldii","Mediterraneibacter faecis","Ruminococcus bromii","Enterococcus avium","Coprococcus phoceensis","Akkermansia muciniphila","Alistipes onderdonkii","Lachnospiraceae bacterium","Eubacterium sp.","Streptococcus salivarius","Clostridia bacterium","Prevotella copri","Bacteroides ovatus","Bacteroides caccae","Faecalibacillus intestinalis","Lacticaseibacillus rhamnosus","Ruminococcus bicirculans","Alistipes putredinis","Dialister invisus","Dorea longicatena","Bifidobacterium adolescentis","Anaerovoracaceae bacterium","Blautia hansenii","Coprococcus eutactus","Escherichia coli","Blautia caecimuris","Subdoligranulum sp.","Hungatella hathewayi","Streptococcus sp."))

#PLOT: INFANTS, FACET: TAXA >= 5 MAGs: color by INDIVIDUAL (230312_03.1_Infant_breadth-adjAGGREGATEpopSNP_FACETEDbyTaxa5ct.pdf)
plot_popSNPs_byLOPthruTP2_ALLtaxa_infants_barebones_individual <- ggplot(data=df_popSNP_ct_perMAG_byLOPthruTP2_infants_5ct, aes(x=LoPthruTP2, y=Aggregate_Adjusted_popSNP_Count)) +
  geom_point(aes(color=Individual), size=1.25, position=position_jitter(h=0.1,w=0.025), alpha=0.5)+
  geom_smooth(method = 'lm', se=FALSE, color="black", size=0.5)+
  stat_regline_equation(aes(label =  paste(..rr.label..)), size = 3)+
  xlab("Duration of persistence")+
  ylab("Breadth-adjusted popSNPs since seeding")+
  ggtitle("INFANTS: Breadth-adjusted popSNP count by Years since seeding")+
  theme_bw()+
  theme(axis.text=element_text(size=8,face="bold"), axis.title=element_text(size=13,face="bold"), panel.grid.minor=element_blank(), legend.title=element_blank(), plot.title=element_text(size=15,face="bold",hjust=0.5), panel.grid.major=element_blank(), panel.border=element_rect(colour="black",size=1))

plot_popSNPs_byLOPthruTP2_ALLtaxa_infants_barebones_individual
plot_popSNPs_byLOPthruTP2_ALLtaxa_infants_barebones_individual + facet_wrap(~ Taxa, ncol=8, scales='free') + theme(strip.text.x = element_text(size=5.5, face = "bold"))

</pre>

# Supplementary Figure 9

<pre>

#------------------------Suppl Figure 9: Aggregate popSNP by Length of persistence (Infant): FACETED BY TAXA, COLORED BY WEANING STATUS------------------------

library(ggplot2)
library(ggprism)
library(ggpubr)
#library(ggExtra)

#TAXA5CT: Read in INFANT dataframe that has breadth-adjusted popSNP ct by LoPthruTP2
df_popSNP_ct_perMAG_byLOPthruTP2_infants_5ct<-read.csv('230311_popSNP_ct_perMAG_byLOPthruTP2_Taxa5ct_infants_diet.csv',
                                                   sep=",",
                                                   header = T)
  
df_popSNP_ct_perMAG_byLOPthruTP2_infants_5ct$Taxa <- factor(df_popSNP_ct_perMAG_byLOPthruTP2_infants_5ct$Taxa, levels = c("Faecalibacterium prausnitzii","Clostridium sp.","Bifidobacterium pseudocatenulatum","Oscillospiraceae bacterium","Bifidobacterium longum","Collinsella sp.","Bifidobacterium bifidum","Sellimonas intestinalis","[Ruminococcus] gnavus","Acutalibacteraceae bacterium","Enterococcus faecalis","Erysipelatoclostridium ramosum","Parabacteroides distasonis","Bacteroides uniformis","Blautia wexlerae","Bifidobacterium breve","Blautia sp.","[Clostridium] innocuum","Phocaeicola vulgatus","Bacteroides fragilis","Anaerostipes hadrus","Anaerostipes caccae","Veillonella parvula","Bacteroides thetaiotaomicron","Flavonifractor plautii","Blautia massiliensis","Alistipes finegoldii","Mediterraneibacter faecis","Ruminococcus bromii","Enterococcus avium","Coprococcus phoceensis","Akkermansia muciniphila","Alistipes onderdonkii","Lachnospiraceae bacterium","Eubacterium sp.","Streptococcus salivarius","Clostridia bacterium","Prevotella copri","Bacteroides ovatus","Bacteroides caccae","Faecalibacillus intestinalis","Lacticaseibacillus rhamnosus","Ruminococcus bicirculans","Alistipes putredinis","Dialister invisus","Dorea longicatena","Bifidobacterium adolescentis","Anaerovoracaceae bacterium","Blautia hansenii","Coprococcus eutactus","Escherichia coli","Blautia caecimuris","Subdoligranulum sp.","Hungatella hathewayi","Streptococcus sp."))

#PLOT: INFANTS, FACET: TAXA >= 5 MAGs: color by DIET (230313_03.3_Infant_breadth-adjAGGREGATEpopSNP_FACETEDbyTaxa5ct_DIET.pdf)
plot_popSNPs_byLOPthruTP2_ALLtaxa_infants_barebones_diet <- ggplot(data=df_popSNP_ct_perMAG_byLOPthruTP2_infants_5ct, aes(x=LoPthruTP2, y=Aggregate_Adjusted_popSNP_Count)) +
  geom_point(aes(color=Group), size=1.25, position=position_jitter(h=0.1,w=0.025), alpha=0.5)+
  scale_color_manual(values=c("#33CCFF","#0066FF","#CC9966"))+
  geom_smooth(method = 'lm', se=FALSE, color="black", size=0.5)+
  stat_regline_equation(aes(label =  paste(..rr.label..)), size = 3)+
  xlab("Duration of persistence")+
  ylab("Breadth-adjusted popSNPs since seeding")+
  ggtitle("INFANTS: Breadth-adjusted popSNP count by Years since seeding")+
  theme_bw()+
  theme(axis.text=element_text(size=8,face="bold"), axis.title=element_text(size=13,face="bold"), panel.grid.minor=element_blank(), legend.title=element_blank(), plot.title=element_text(size=15,face="bold",hjust=0.5), panel.grid.major=element_blank(), panel.border=element_rect(colour="black",size=1))

plot_popSNPs_byLOPthruTP2_ALLtaxa_infants_barebones_diet
plot_popSNPs_byLOPthruTP2_ALLtaxa_infants_barebones_diet + facet_wrap(~ Taxa, ncol=8, scales='free') + theme(strip.text.x = element_text(size=5.5, face = "bold"))

</pre>

# Supplementary Figure 10

<pre>

#------------------------Suppl Figure 10: Aggregate popSNP by Length of persistence (Mother): FACETED BY TAXA, COLORED BY INDIVIDUAL------------------------

#TAXA5CT: Read in MOTHER dataframe that has breadth-adjusted popSNP ct by LoPthruTP2
df_popSNP_ct_perMAG_byLOPthruTP2_mothers_5ct<-read.csv('230311_popSNP_ct_perMAG_byLOPthruTP2_Taxa5ct_mothers.csv',
                                                       sep=",",
                                                       header = T)
df_popSNP_ct_perMAG_byLOPthruTP2_mothers_5ct$Taxa <- factor(df_popSNP_ct_perMAG_byLOPthruTP2_mothers_5ct$Taxa, levels = c("Acutalibacteraceae bacterium","Phocaeicola vulgatus","[Ruminococcus] torques","Lachnospiraceae bacterium","Faecalibacterium prausnitzii","Oscillospiraceae bacterium","Bacteroides uniformis","Blautia massiliensis","Anaerostipes hadrus","Clostridium sp.","Bifidobacterium longum","Parabacteroides distasonis","Collinsella sp.","Dorea longicatena","Coprococcus comes","Alistipes putredinis","Blautia sp.","Agathobacter rectalis","Anaerobutyricum hallii","Bifidobacterium adolescentis","Ruminococcus sp."))

#PLOT: MOTHERS
plot_popSNPs_byLOPthruTP2_ALLtaxa_mothers_barebones <- ggplot(data=df_popSNP_ct_perMAG_byLOPthruTP2_mothers_5ct, aes(x=LoPthruTP2, y=Aggregate_Adjusted_popSNP_Count)) +
  geom_point(aes(color=Individual), size=1.5, position=position_jitter(h=0.3,w=0.02, seed=8), alpha=0.5)+
  geom_smooth(method = 'lm', se=FALSE, color="black", size=0.5)+
  stat_regline_equation(aes(label =  paste(..rr.label.., sep = "~~~~")), size = 4)+
  xlab("Duration of persistence")+
  ylab("popSNPs accrued")+
  ggtitle("MOTHERS: popSNP Accumulation by Time")+
  theme_bw()+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=13,face="bold"), panel.grid.minor=element_blank(), legend.title=element_blank(), plot.title=element_text(size=15,face="bold",hjust=0.5), panel.grid.major=element_blank(), panel.border=element_rect(colour="black",size=1))

plot_popSNPs_byLOPthruTP2_ALLtaxa_mothers_barebones + facet_wrap(~ Taxa, ncol=5, scales='free') + theme(strip.text.x = element_text(face = "bold"))

</pre>

# Supplementary Figure 11A

<pre>
#------------------------Suppl Figure 11A1: Aggregate popSNP by Length of persistence - Mothers, linear regression------------------------

library(ggplot2)
library(ggprism)
library(ggpubr)

#MOTHERS: ALL TAXA: Breadth-adjusted Aggregate popSNP ct by LoPthruTP2 (230311_02.1_Mother_breadth-adjAGGREGATEpopSNP_alltaxa.pdf)

df_popSNP_ct_perMAG_byLOPthruTP2_mothers_alltaxa<-read.csv('230314_popSNP_ct_perMAG_byLOPthruTP2_mothers_gs.csv',
                                                   sep=",",
                                                   header = T)
  
ggplot(data=df_popSNP_ct_perMAG_byLOPthruTP2_mothers_alltaxa, aes(x=LoPthruTP2, y=Aggregate_Adjusted_popSNP_Count)) +
  geom_point(size=1, position=position_jitter(h=0.2,w=0.01, seed=7), alpha=0.25, fill="black")+
  geom_smooth(method = 'lm', color="black", fill="#996666")+
  stat_regline_equation(aes(label =  paste(..eq.label..,..rr.label.., sep = "~~~~")),label.x = c(0.02), label.y = c(35.25), size = 4.75)+
  scale_x_continuous(limits=c(-.02,3.1), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(-0.5,40), expand=c(0,0), guide = guide_prism_minor())+
  annotate("text",x=0.41,y=38,label="n = 330 MAGs (130 taxa)", size=4.75)+
  #annotate("text",x=0.4,y=45.75,label="n = 140 MAGs (21 taxa with > 4 persisting MAGs)")+
  xlab("Length of persistence")+
  ylab("Aggregate popSNPs since seeding")+
  ggtitle("Generalized mutation rate in mothers (Local regression)")+
  theme_bw()+
  theme(axis.text=element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"), panel.grid.minor=element_blank(), legend.title=element_blank(), plot.title=element_text(size=19,face="bold",hjust=0.5), panel.grid.major=element_blank(), panel.border=element_rect(colour="black",size=1))


#------------------------Suppl Figure 11A2: Aggregate popSNP by Length of persistence - Infants, linear regression------------------------

#INFANTS: ALL TAXA: Breadth-adjusted Aggregate popSNP ct by LoPthruTP2 (230311_01.1_Infant_breadth-adjAGGREGATEpopSNP_alltaxa.pdf)
df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa<-read.csv('230314_popSNP_ct_perMAG_byLOPthruTP2_infants_gs.csv', sep=",", header = T)

ggplot(data=df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa, aes(x=LoPthruTP2, y=Aggregate_Adjusted_popSNP_Count)) +
  #geom_point(size=1.5, alpha=0.5, fill="black")+
  geom_point(shape=21, size=1, position=position_jitter(h=0.3,w=0.025), alpha=0.25, fill="black")+
  geom_smooth(method = 'lm', color="black", fill="#996666")+
  stat_regline_equation(aes(label =  paste(..eq.label..,..rr.label.., sep = "~~~~")),label.x = c(0.05), label.y = c(119), size = 4.75)+
  scale_x_continuous(limits=c(-.085,8), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(-2,135), expand=c(0,0), guide = guide_prism_minor())+
  annotate("text",x=1.07,y=128,label="n = 763 MAGs (194 taxa)", size = 4.75)+
  xlab("Length of persistence")+
  ylab("Aggregate popSNPs since seeding")+
  ggtitle("Generalized mutation rate in infants (Linear regression)")+
  theme_bw()+
  theme(axis.text=element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"), panel.grid.minor=element_blank(), legend.title=element_blank(), plot.title=element_text(size=19,face="bold",hjust=0.5), panel.grid.major=element_blank(), panel.border=element_rect(colour="black",size=1))

</pre>

# Supplementary Figure 11B

<pre>
#------------------------Suppl Figure 11B1: Aggregate popSNP by Length of persistence - Mothers, local regression------------------------

library(ggplot2)
library(ggprism)
library(ggpubr)

#MOTHERS: ALL TAXA: Breadth-adjusted Aggregate popSNP ct by LoPthruTP2 (230311_02.1_Mother_breadth-adjAGGREGATEpopSNP_alltaxa.pdf)

df_popSNP_ct_perMAG_byLOPthruTP2_mothers_alltaxa<-read.csv('230314_popSNP_ct_perMAG_byLOPthruTP2_mothers_gs.csv',
                                                   sep=",",
                                                   header = T)
ggplot(data=df_popSNP_ct_perMAG_byLOPthruTP2_mothers_alltaxa, aes(x=LoPthruTP2, y=Aggregate_Adjusted_popSNP_Count)) +
  geom_point(size=1, position=position_jitter(h=0.2,w=0.01, seed=7), alpha=0.25, fill="black")+
  geom_smooth(method = 'loess', color="black", fill="#3366CC")+
  scale_x_continuous(limits=c(-.02,3.1), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(-0.5,40), expand=c(0,0), guide = guide_prism_minor())+
  annotate("text",x=0.41,y=38,label="n = 330 MAGs (130 taxa)", size=4.75)+
  xlab("Length of persistence")+
  ylab("Aggregate popSNPs since seeding")+
  ggtitle("Generalized mutation rate in mothers (Local regression)")+
  theme_bw()+
  theme(axis.text=element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"), panel.grid.minor=element_blank(), legend.title=element_blank(), plot.title=element_text(size=19,face="bold",hjust=0.5), panel.grid.major=element_blank(), panel.border=element_rect(colour="black",size=1))

#------------------------Suppl Figure 11B2: Aggregate popSNP by Length of persistence - Infants, local regression------------------------


#INFANTS: ALL TAXA: Breadth-adjusted Aggregate popSNP ct by LoPthruTP2 (230311_01.1_Infant_breadth-adjAGGREGATEpopSNP_alltaxa.pdf)
df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa<-read.csv('230314_popSNP_ct_perMAG_byLOPthruTP2_infants_gs.csv', sep=",", header = T)

ggplot(data=df_popSNP_ct_perMAG_byLOPthruTP2_infants_alltaxa, aes(x=LoPthruTP2, y=Aggregate_Adjusted_popSNP_Count)) +
  geom_point(shape=21, size=1, position=position_jitter(h=0.3,w=0.025), alpha=0.25, fill="black")+
  geom_smooth(method = 'loess', color="black", fill="#3366CC")+
  scale_x_continuous(limits=c(-.085,8), expand=c(0,0), guide = guide_prism_minor())+
  scale_y_continuous(limits=c(-2,135), expand=c(0,0), guide = guide_prism_minor())+
  annotate("text",x=1.07,y=128,label="n = 763 MAGs (194 taxa)", size = 4.75)+
  xlab("Length of persistence")+
  ylab("Aggregate popSNPs since seeding")+
  ggtitle("Generalized mutation rate in infants (Local regression)")+
  theme_bw()+
  theme(axis.text=element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"), panel.grid.minor=element_blank(), legend.title=element_blank(), plot.title=element_text(size=19,face="bold",hjust=0.5), panel.grid.major=element_blank(), panel.border=element_rect(colour="black",size=1))

</pre>

# Supplementary Figure 12

Made in Prism 9.

# Supplementary Figure 13

<pre>
#------------------------Suppl Figure 13: Boxplot comparing Bray-Curtis distance between Pre-Weaning, Post-Weaning, and Mother mutated gene COG profiles, split by family------------------------

#read in boxplot csv
df_pairwise_boxplot<-read.csv('230322_bray_pairwise_boxplot_2groups.csv',
                              sep=",",
                              header = T)

library(ggplot2)
library(ggpubr)
library(ggprism)

df_pairwise_boxplot$Comparison_type___Family <- factor(df_pairwise_boxplot$Comparison_type___Family,levels = c("before___mother___same","before___mother___different","after___mother___same","after___mother___different"))
df_pairwise_boxplot$Family <- factor(df_pairwise_boxplot$Family,levels = c("Same","Different"))
df_pairwise_boxplot$Weaning <- factor(df_pairwise_boxplot$Weaning,levels = c("Pre","Post"))

ggplot(df_pairwise_boxplot, aes(Weaning, Bray_Curtis_distance)) +
  geom_boxplot(aes(colour = Family), outlier.shape = NA, lwd=0.75)+
  geom_point(aes(fill = Family, shape = Weaning), size = 5, alpha=0.75, position = position_jitterdodge(jitter.width = 0.75, jitter.height = 0, seed=10))+
  scale_x_discrete(labels=c("Pre-Weaning Infant vs. Mother","Post-Weaning Infant vs. Mother"))+
  scale_color_manual(values=c("#ED774D","#263C57"))+
  scale_fill_manual(values=c("#ED774D","#263C57"))+
  scale_shape_manual(values=c(21, 23))+
  ylab("Bray-Curtis dissimilarity")+
  ggtitle("Distance between COG functional profiles of mutated genes")+
  stat_compare_means(method = "kruskal.test", aes(group = Family), label = "p.signif", label.y = 0.45)+
  theme_bw()+
  theme(plot.title=element_text(size=17,face="bold", hjust=0.5), axis.title.x=element_blank(), axis.title.y=element_text(size=13, face="bold"), axis.text = element_text(face = "bold",size = 12), legend.title=element_text(face="bold"), legend.text=element_text(face="bold"), panel.border=element_rect(colour="black",fill=NA,size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  theme(axis.title=element_text(size=14,face="bold"), axis.text=element_blank(), axis.ticks=element_blank(), panel.border=element_rect(colour="black",fill=NA,size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title=element_text(face="bold"), legend.text=element_text(face="bold"), plot.title = element_text(hjust = 0.25,size=16,face="bold"))

</pre>

