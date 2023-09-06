# twins_diet
Code used for data visualization and analysis in Sawhney et al (Weaning accelerates and transforms within-host adaptation in the infant gut microbiome).

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
