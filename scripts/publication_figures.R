## Exploratory analysis and publication figure generation for
## Meta-MARC HMM data.  Test set data for Soil, Pediatric, and NCBA1
## are in long format in the analytic_data/mmarc_publication_analytic_data.csv
## file, and the simulation output data is in the
## analytic_data/mmarc_publication_simulation_results.csv file.

library(ggplot2)
library(RColorBrewer)
library(data.table)

setwd("~/Documents/morleyBioinformatics/GigaScience_HMM_Publication/publication")

abund_data <- data.table(read.csv('analytic_data/mmarc_publication_analytic_data.csv'))
ncba_metadata <- data.table(read.csv('sample_lists/PRJNA292471_Metadata.csv'))

pattern <- unique(abund_data$TruthLabel[!is.na(abund_data$TruthLabel)])
replace <- as.character(c(
    "betalactams",
    "Phenicol",
    "betalactams",
    "Glycopeptides",       
    "betalactams",
    "Aminoglycosides",
    "betalactams",
    "betalactams",
    "betalactams",
    "Tetracyclines",
    "Trimethoprim",
    "Trimethoprim",
    "Tetracyclines",
    "betalactams",
    "betalactams",
    "betalactams",
    "Tetracyclines",
    "Tetracyclines",
    "betalactams",
    "betalactams",
    "Phenicol",
    "betalactams",
    "betalactams",
    "Glycopeptides",
    "betalactams",
    "Aminoglycosides",
    "Trimethoprim",
    "Trimethoprim"
))

transl_table <- data.table(TruthLabel=pattern, AnalyticTruthLabel=replace)
setkey(transl_table, TruthLabel)

group2_data <- abund_data[HMMGroup == 'groupII', ]
g2filter_data <- abund_data[HMMGroup != 'groupII' | is.na(HMMGroup), ]

class_data <- g2filter_data[AnnotationLevel == 'Class', ]

dantas_class_data <- class_data[TestSet != 'NCBA', ]
setkey(dantas_class_data, TruthLabel)

dantas_class_data <- transl_table[dantas_class_data]

norm_dantas_class_data <- dantas_class_data[, NormAbundance := (100 * Abundance / sum(Abundance) ), by=c('SampleName', 'Pipeline')]


# ## Dantas, at the sample level
# dantas_sample_names <- as.character(unique(norm_dantas_class_data$SampleName))
# for( s in 1:length(dantas_sample_names)) {
#     temptruth <- as.character(unique(norm_dantas_class_data[SampleName == dantas_sample_names[s]]$AnalyticTruthLabel))
#     temp <- norm_dantas_class_data[SampleName == dantas_sample_names[s] & NormAbundance >= 1, .SD, .SDcols=c('NodeName', 'NormAbundance', 'Pipeline')]
#     temp2 <- temp[, sum(NormAbundance), by=c('Pipeline', 'NodeName')]
#     png(paste('graphs/dantas_samples/', dantas_sample_names[s], '.png', collapse='', sep=''), width=900, height=900)
#     g <- ggplot(temp2, aes(x=NodeName, y=V1, fill=Pipeline)) +
#         geom_bar(stat='identity', position='dodge') +
#         theme(panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               strip.text.x=element_text(size=17),
#               strip.text.y=element_text(size=17, angle=0),
#               axis.text.y=element_text(size=16),
#               axis.text.x=element_text(size=16, angle=45, hjust=1),
#               axis.title.x=element_text(size=20),
#               axis.title.y=element_text(size=20),
#               legend.position="bottom",
#               panel.margin=unit(0.1, "lines"),
#               plot.title=element_text(size=20, hjust=0.5),
#               legend.text=element_text(size=14),
#               legend.title=element_blank()) +
#         ggtitle(paste('AMR Class Abundance by Pipeline\nSample ', dantas_sample_names[s], ', Truth Label: ', temptruth)) +
#         xlab('\nAMR Drug Class') + ylab('Relative Abundance (%)\n')
#     print(g)
#     dev.off()
# }
# 
## Dantas all classes
# temp <- norm_dantas_class_data[, .SD, .SDcols=c('NodeName', 'NormAbundance', 'Pipeline')]
# temp2 <- temp[, sum(NormAbundance) / length(unique(dantas_class_data$SampleName)), by=c('Pipeline', 'NodeName')]
# #temp2[, TotalRelativeAbundance := (SumNormAbundance), by='Pipeline']
# temp3 <- temp3[V1 > 1, ]
# temp3 <- rbind(temp3, data.table(Pipeline=c('Alignment', 'Resfams', 'Alignment'), NodeName=c(rep('Aminocoumarins', 2), 'Glycopeptides'), V1=rep(0, 3)))
# #temp2$TotalRelativeAbundance <- temp2$TotalRelativeAbundance / length(unique(dantas_class_data$SampleName))
# png('graphs/dantas_class_all_samples.png', width=1200, height=900)
# g <- ggplot(temp3, aes(x=NodeName, y=V1, fill=Pipeline)) +
#     geom_bar(stat='identity', position='dodge') +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           strip.text.x=element_text(size=20),
#           strip.text.y=element_text(size=20, angle=0),
#           axis.text.y=element_text(size=20),
#           axis.text.x=element_text(size=20, angle=45, hjust=1),
#           axis.title.x=element_text(size=24),
#           axis.title.y=element_text(size=24),
#           legend.position="bottom",
#           panel.margin=unit(0.1, "lines"),
#           plot.title=element_text(size=30, hjust=0.5),
#           legend.text=element_text(size=18),
#           legend.title=element_blank()) +
#     xlab('\nAMR Drug Class') + ylab('Mean Relative Abundance (%)\n') +
#     scale_fill_discrete(drop=F) + scale_fill_brewer(palette = 'Set2')
# print(g + ggtitle(paste('Mean Sample-wise AMR Class Abundance by Pipeline\nSoil and Pediatric Test Sets')))
# dev.off()
# 
# png('graphs/dantas_class_all_samples_notitle.png', width=1200, height=900)
# print(g)
# dev.off()


## MMARC performance on Dantas all classes
mmarc_dantas_class_data <- dantas_class_data[Pipeline == 'MMARC', ]
performance_mmarc_dantas_class_data <-
    mmarc_dantas_class_data[, PerformanceCount := ( 
        ifelse(NodeName == AnalyticTruthLabel | NodeName == 'Multi-drug resistance', Abundance, 0) )]
performance_mmarc_dantas_class_data[, PerformancePercent := ( 100 * sum(PerformanceCount) / sum(Abundance) ), by='SampleName']
performance_mmarc_dantas_class_data <- performance_mmarc_dantas_class_data[, .SD, .SDcols=c('SampleName', 'TestSet', 'PerformancePercent', 'AnalyticTruthLabel')]
setkey(performance_mmarc_dantas_class_data, SampleName)
performance_mmarc_dantas_class_data <- unique(performance_mmarc_dantas_class_data)


png('graphs/mmarc_dantas_performance_class_all_samples.png', width=1200, height=900)
g <- ggplot(performance_mmarc_dantas_class_data, aes(x=AnalyticTruthLabel, y=PerformancePercent)) +
    geom_boxplot() + geom_jitter(width=0.1, height=0.1, size=1) + facet_wrap(~TestSet) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x=element_text(size=20),
          strip.text.y=element_text(size=20, angle=0),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=20, angle=45, hjust=1),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.position="bottom",
          panel.margin=unit(0.1, "lines"),
          plot.title=element_text(size=30, hjust=0.5),
          legend.text=element_text(size=18),
          legend.title=element_blank()) +
    xlab('\nFunctional Metagenomic Truth Label') + ylab('Percent of Reads On-Target\n')
print(g + ggtitle(paste('Meta-MARC Percent of Hits On-Target by Truth Label\nSoil and Pediatric Test Sets')))
dev.off()

png('graphs/mmarc_dantas_performance_class_all_samples_notitle.png', width=1200, height=900)
print(g)
dev.off()





