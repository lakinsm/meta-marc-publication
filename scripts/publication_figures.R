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
    "Phenicol",
    "betalactams",
    "betalactams",
    "betalactams",
    "Glycopeptides",
    "betalactams",
    "Aminoglycosides",
    "betalactams",
    "betalactams",
    "Tetracyclines",
    "Tetracyclines",
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


## Dantas, at the sample level
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

# Dantas all classes
temp <- norm_dantas_class_data[, .SD, .SDcols=c('NodeName', 'NormAbundance', 'Pipeline')]
temp2 <- temp[, sum(NormAbundance) / length(unique(dantas_class_data$SampleName)), by=c('Pipeline', 'NodeName')]
#temp2[, TotalRelativeAbundance := (SumNormAbundance), by='Pipeline']
temp3 <- temp2[V1 > 1, ]
temp3 <- rbind(temp3, data.table(Pipeline=c('Alignment', 'Resfams', 'Alignment'), NodeName=c(rep('Aminocoumarins', 2), 'Glycopeptides'), V1=rep(0, 3)))
#temp2$TotalRelativeAbundance <- temp2$TotalRelativeAbundance / length(unique(dantas_class_data$SampleName))
png('graphs/dantas_class_all_samples.png', width=1200, height=900)
g <- ggplot(temp3, aes(x=NodeName, y=V1, fill=Pipeline)) +
    geom_bar(stat='identity', position='dodge') +
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
    xlab('\nAMR Drug Class') + ylab('Mean Relative Abundance (%)\n') +
    scale_fill_discrete(drop=F) + scale_fill_brewer(palette = 'Set2')
print(g + ggtitle(paste('Mean Sample-wise AMR Class Abundance by Pipeline\nSoil and Pediatric Test Sets')))
dev.off()

png('graphs/dantas_class_all_samples_notitle.png', width=1200, height=900)
print(g)
dev.off()


## MMARC performance on Dantas all classes
performance_mmarc_dantas_class_data <-
    dantas_class_data[, PerformanceCount := ( 
        ifelse(NodeName == AnalyticTruthLabel | NodeName == 'Multi-drug resistance', Abundance, 0) )]
performance_mmarc_dantas_class_data[, PerformanceTotalCount := ( sum(PerformanceCount) ), by=c('SampleName', 'Pipeline')]
performance_mmarc_dantas_class_data[, TotalCountsOverall := ( sum(Abundance) ), by=c('SampleName', 'Pipeline')]
performance_mmarc_dantas_class_data[, PerformancePercent := ( 100 * sum(PerformanceCount) / sum(Abundance) ), by=c('SampleName', 'Pipeline')]
performance_mmarc_dantas_class_data <- performance_mmarc_dantas_class_data[, .SD, .SDcols=c('Pipeline','TestSet', 'SampleName', 'AnalyticTruthLabel', 'PerformanceTotalCount', 'TotalCountsOverall', 'PerformancePercent')]
setkeyv(performance_mmarc_dantas_class_data, c('SampleName', 'Pipeline'))
performance_mmarc_dantas_class_data <- unique(performance_mmarc_dantas_class_data)

dantas_sample_names <- as.character(unique(norm_dantas_class_data$SampleName))
pipelines <- as.character(unique(norm_dantas_class_data$Pipeline))

truth_labels <- performance_mmarc_dantas_class_data[, .SD, .SDcols=c('SampleName', 'AnalyticTruthLabel', 'TestSet')]
setkey(truth_labels, SampleName)
truth_labels <- unique(truth_labels)

zero_entries = data.table(Pipeline=character(),
                          TestSet=character(),
                          SampleName=character(),
                          AnalyticTruthLabel=character(),
                          PerformanceTotalCount=numeric(),
                          TotalCountsOverall=numeric(),
                          PerformancePercent=numeric())
for( s in 1:length(dantas_sample_names) ) {
    for( p in 1:length(pipelines) ) {
        label_name <- as.character(truth_labels$AnalyticTruthLabel[truth_labels$SampleName == dantas_sample_names[s]])
        set_name <- as.character(truth_labels$TestSet[truth_labels$SampleName == dantas_sample_names[s]])
        temp <- performance_mmarc_dantas_class_data[Pipeline == pipelines[p] &
                                                        AnalyticTruthLabel == label_name &
                                                        SampleName == dantas_sample_names[s],]
        if( nrow(temp) == 0 ) {
            zero_entries <- rbind(zero_entries, data.table(Pipeline=pipelines[p],
                                                           TestSet=set_name,
                                                           SampleName=dantas_sample_names[s],
                                                           AnalyticTruthLabel=label_name,
                                                           PerformanceTotalCount=0,
                                                           TotalCountsOverall=0,
                                                           PerformancePercent=0))
        }
    }
}

performance_mmarc_dantas_class_data <- rbind(performance_mmarc_dantas_class_data, zero_entries)


png('graphs/mmarc_dantas_performance_class_all_samples.png', width=1200, height=900)
g <- ggplot(performance_mmarc_dantas_class_data, aes(x=AnalyticTruthLabel, y=PerformancePercent)) +
    geom_boxplot() + geom_jitter(width=0.1, height=0, size=1) +
    facet_grid(Pipeline ~ TestSet, switch='y') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x=element_text(size=20),
          strip.text.y=element_text(size=20, angle=180),
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
print(g + ggtitle(paste('Percent of Hits On-Target for each Truth Label, by Test Set and Pipeline\n')))
dev.off()

png('graphs/mmarc_dantas_performance_class_all_samples_notitle.png', width=1200, height=900)
print(g)
dev.off()

write.table(performance_mmarc_dantas_class_data,
            'graphs/mmarc_performance_class_all_samples_table.csv',
            sep=',', col.names = T, row.names = F)


## Using the data from above, graph the number of reads found with alignment vs. MMARC
## for each truth label category
temp <- performance_mmarc_dantas_class_data
colnames(temp) <- c('Pipeline', 'TestSet', 'SampleName', 'AnalyticTruthLabel',
                    'NumberReadsOnTarget', 'TotalReadsIdentified', 'PerformancePercent')

temp <- temp[Pipeline != 'Resfams', .SD, .SDcols=c('Pipeline', 'TestSet', 'AnalyticTruthLabel',
                                                   'NumberReadsOnTarget', 'TotalReadsIdentified')]
temp <- temp[, lapply(.SD, sum), by=c('Pipeline', 'TestSet', 'AnalyticTruthLabel')]
# performance_raw <- melt(temp,
#                         id.vars=c('Pipeline', 'TestSet', 'AnalyticTruthLabel'),
#                         measure.vars=c('NumberReadsOnTarget', 'TotalReadsIdentified'),
#                         variable.name='AbundanceCategory',
#                         value.name='Abundance')
# performance_raw <- performance_raw[Pipeline != 'Resfams', ]

png('graphs/mmarc_dantas_abundance_class_all_samples.png', width=1200, height=900)
g <- ggplot(temp, aes(x=AnalyticTruthLabel, fill=Pipeline)) +
    geom_bar(data=temp, aes(y=TotalReadsIdentified), position='dodge', stat='identity') +
    geom_point(data=temp, aes(y=NumberReadsOnTarget), position=position_dodge(width=0.9), shape=3, size=6) +
    facet_wrap(~TestSet, nrow=1, ncol=2, scale='free_y') +
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
    xlab('\nFunctional Metagenomic Truth Label') + ylab('Number of Reads Identified\n') +
    scale_fill_brewer(palette = 'Set2')
print(g + ggtitle(paste('Total and On-Target Abundance for each Truth Label, by Test Set and Pipeline\n')))
dev.off()

png('graphs/mmarc_dantas_abundance_class_all_samples_notitle.png', width=1200, height=900)
print(g)
dev.off()

############################
############################

## Simulation data
sim_data <- data.table(read.csv('analytic_data/mmarc_publication_simulation_results.csv'))
class_sim_data <- sim_data[Level == 'Classes']
class_sim_data$Level <- gsub('Classes', 'Class', class_sim_data$Level)
colnames(class_sim_data) <- c('NodeName', 'AnnotationLevel', 'TruePositives', 'FalsePositives', 'FalseNegatives',
                              'TrueNegatives', 'Sensitivity', 'Specificity')
class_sim_data$NodeName <- gsub('MLS', 'Macrolides Lincosamides and Streptogramins', class_sim_data$NodeName)
class_sim_data <- class_sim_data[, .SD, .SDcols=!c('AnnotationLevel', 'Specificity', 'Sensitivity')]
setkey(class_sim_data, NodeName)

class_sim_data[, Sensitivity := ( 100 * TruePositives / (TruePositives + FalseNegatives) )]
class_sim_data[, Specificity := ( 100 * TrueNegatives / (TrueNegatives + FalsePositives) )]

write.table(class_sim_data, 'graphs/mmarc_class_simulation_metrics.csv', sep=',', row.names=F, col.names=T)

#########################
#########################
## NCBA data with HMMs and alignment

ncba_data <- abund_data[TestSet == 'NCBA', ]
colnames(ncba_metadata)[4] <- 'SampleName'
setkey(ncba_data, SampleName)
setkey(ncba_metadata, SampleName)

ncba_data <- ncba_metadata[ncba_data]





