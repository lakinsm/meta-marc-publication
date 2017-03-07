library(ggplot2)
library(RColorBrewer)
library(data.table)
library(scales)
library(arm)

setwd('/home/lakinsm/Documents/morleyBioinformatics/GigaScience_HMM_Publication/mmarc_production_repositories/meta-marc-publication/')

rdata <- data.table(read.csv('analytic_data/mismatch_analytic_data.csv',
                             col.names=c('TestSet', 'SampleName', 'DantasSampleSet', 'TruthLabel', 'Pipeline', 'HMMGroup',
                                         'AnnotationLevel', 'NodeName', 'MinorAlleleCount', 'ReadCount')))

rdata_class <- rdata[!is.na(NodeName) & AnnotationLevel == 'Class' & (HMMGroup != 'groupII' | is.na(HMMGroup))]

rdata_class[, AvgMinorAlleleCount := ( MinorAlleleCount / ReadCount )]

df1 <- data.frame(a = c(1, 1, 2, 2), b = c(48, 50, 50, 48))
df2 <- data.frame(a = c(2, 2, 3, 3), b = c(52, 54, 54, 52))
df3 <- data.frame(a = c(2, 2, 4, 4), b = c(56, 58, 58, 56))


png('graphs/ncba_mismatch_by_pipeline.png', width=1200, height=900)
g <- ggplot(rdata_class, aes(x=Pipeline, y=AvgMinorAlleleCount)) + geom_boxplot(outlier.shape = NA) +
    #geom_jitter(width=0.3, height=0, size=1, alpha=0.1) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x=element_text(size=20),
          strip.text.y=element_text(size=20, angle=180),
          axis.text.y=element_text(size=22),
          axis.text.x=element_text(size=22, angle=45, hjust=1),
          axis.title.x=element_text(size=26),
          axis.title.y=element_text(size=26),
          panel.margin=unit(0.1, "lines"),
          plot.title=element_text(size=30, hjust=0.5),
          legend.text=element_text(size=18),
          legend.title=element_blank()) +
    xlab('\nPipeline') + ylab('Average Non-major Allele Count per Read\n') +
    ylim(c(-2, 60)) +
    geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 51, label = "***", size = 8) +
    geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 55, label = "***", size = 8) +
    geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 3, y = 59, label = "***", size = 8)
print(g + ggtitle(paste('Average Read-wise Genetic Variation by Pipeline\n')))
dev.off()

png('graphs/ncba_mismatch_by_pipeline_notext.png', width=1200, height=900)
print(g)
dev.off()

rdata_class[, FamilyMethod := ( ifelse(Pipeline == 'MMARC', 'HMMFASTQ', 'Alignment'))]



glmer(AvgMinorAlleleCount ~ I(FamilyMethod) + ( 1 | NodeName ), family = 'gaussian', data = rdata_class)

rcast <- dcast(rdata_class, SampleName + NodeName ~ FamilyMethod, value.var='AvgMinorAlleleCount', fill=NA, fun=mean)
rcast2 <- dcast(rdata_class, SampleName + NodeName ~ Pipeline, value.var='AvgMinorAlleleCount', fill=NA, fun=mean)

# Overal HMM vs alignment
wilcox.test(rcast$Alignment, rcast$HMMFASTQ)

# Pairwise
wilcox.test(rcast2$MMARC, rcast2[['MMARC-assembled']])$p.value * 3
wilcox.test(rcast2$MMARC, rcast2[['Resfams']])$p.value * 3
wilcox.test(rcast2$MMARC, rcast2[['Alignment']])$p.value * 3





