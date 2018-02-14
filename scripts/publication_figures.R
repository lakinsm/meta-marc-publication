## Exploratory analysis and publication figure generation for
## Meta-MARC HMM data.  Test set data for Soil, Pediatric, and NCBA1
## are in long format in the analytic_data/mmarc_publication_analytic_data.csv
## file, and the simulation output data is in the
## analytic_data/mmarc_publication_simulation_results.csv file.

library(ggplot2)
library(RColorBrewer)
library(data.table)
library(scales)

setwd("/mnt/manuscripts/CommunicationsBiology_HMM_Publication/meta-marc-publication/")

abund_data <- data.table(read.csv('analytic_data/mmarc_publication_analytic_data.csv'))
ncba_metadata <- data.table(read.csv('metadata/PRJNA292471_metadata.csv'))



colnames(abund_data) <- c('TestSet', 'SampleName', 'DantasSampleSet', 'TruthLabel', 'Pipeline', 'HMMGroup',
                          'AnnotationLevel', 'NodeName', 'Abundance')


abund_data$NodeName <- gsub('MLS', 'Macrolides, Lincosamides,\nStreptogrammins', as.character(abund_data$NodeName))
abund_data$NodeName <- gsub('betalactams', 'beta-lactams', as.character(abund_data$NodeName))
abund_data$NodeName <- gsub('Cationic antimicrobial peptides',
                            'Cationic Antimicrobial Peptides', as.character(abund_data$NodeName))
abund_data$NodeName <- gsub('Multi-drug resistance',
                            'Multi-drug Resistance', as.character(abund_data$NodeName))

pattern <- unique(abund_data$TruthLabel[!is.na(abund_data$TruthLabel)])
replace <- as.character(c(
    "beta-lactams",
    "Phenicol",
    "beta-lactams",
    "Glycopeptides",
    "beta-lactams",
    "Aminoglycosides",
    "beta-lactams",
    "beta-lactams",
    "beta-lactams",
    "Tetracyclines",
    "Trimethoprim",
    "Trimethoprim",
    "Tetracyclines",
    "Phenicol",
    "beta-lactams",
    "Glycopeptides",
    "beta-lactams",
    "Tetracyclines",
    "beta-lactams",
    "Tetracyclines",
    "Trimethoprim",
    "beta-lactams",
    "beta-lactams",
    "beta-lactams",
    "beta-lactams",
    "Aminoglycosides"
))



transl_table <- data.table(TruthLabel=pattern, AnalyticTruthLabel=replace)
setkey(transl_table, TruthLabel)

group2_data <- abund_data[HMMGroup == 'groupII', ]
g2filter_data <- abund_data[HMMGroup != 'groupII' | is.na(HMMGroup), ]

class_data <- g2filter_data[AnnotationLevel == 'Class', ]

dantas_class_data <- class_data[TestSet != 'NCBA', ]
setkey(dantas_class_data, TruthLabel)

dantas_class_data <- transl_table[dantas_class_data]

#norm_dantas_class_data <- dantas_class_data[, NormAbundance := (100 * Abundance / sum(Abundance) ), by=c('SampleName', 'Pipeline')]


## Dantas, at the sample level
# dantas_sample_names <- as.character(unique(dantas_class_data$SampleName))
# for( s in 1:length(dantas_sample_names)) {
#     temptruth <- as.character(unique(dantas_class_data[SampleName == dantas_sample_names[s]]$AnalyticTruthLabel))
#     temp <- dantas_class_data[SampleName == dantas_sample_names[s] & Abundance >= 100, .SD, .SDcols=c('NodeName', 'Abundance', 'Pipeline')]
#     temp2 <- temp[, sum(Abundance), by=c('Pipeline', 'NodeName')]
#     
#     
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
temp <- dantas_class_data[, .SD, .SDcols=c('NodeName', 'Abundance', 'Pipeline')]
temp <- temp[, SumAbundance := (lapply(.SD, sum)), by=c('NodeName', 'Pipeline')]
setkeyv(temp, c('NodeName', 'Pipeline'))
temp <- unique(temp[, .SD, .SDcols=c('NodeName', 'Pipeline', 'SumAbundance')])
#temp2[, TotalRelativeAbundance := (SumNormAbundance), by='Pipeline']

threshold <- c()
for( x in 1:length(unique(temp$NodeName))) {
    ltemp <- temp[NodeName == as.character(unique(temp$NodeName)[x])]
    if( any(ltemp$SumAbundance > 2e+06) ) {
        threshold <- c(threshold, as.character(unique(temp$NodeName)[x]))
    }
}
temp3 <- temp[NodeName %in% threshold]

pipenames <- c('MMARC', 'MMARC-assembled', 'Alignment', 'Resfams')
for( x in 1:length(unique(temp3$NodeName)) ) {
    ltemp <- temp3[NodeName == as.character(unique(temp3$NodeName)[x])]
    for( p in pipenames ) {
        if( ! p %in% unique(as.character(ltemp$Pipeline)) ) {
            temp3 <- rbind(temp3, data.table(NodeName=as.character(unique(temp3$NodeName)[x]),
                                             Pipeline=p,
                                             SumAbundance=0))
        }
    }
}

png('graphs/dantas_class_all_samples.png', width=1200, height=900)
g <- ggplot(temp3, aes(x=NodeName, y=SumAbundance, fill=Pipeline)) +
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
    xlab('\nAMR Drug Class') + ylab('Abundance\n') +
    scale_fill_brewer(palette = 'Set2', drop=F) + scale_y_continuous(labels = comma)
print(g + ggtitle(paste('AMR Class Abundance by Pipeline\nSoil and Pediatric Test Sets')))
dev.off()

png('graphs/dantas_class_all_samples_notitle.png', width=1200, height=900)
print(g)
dev.off()

performance_mmarc_dantas_class_data <-
    dantas_class_data[, PerformanceCount := (
        ifelse(NodeName == AnalyticTruthLabel | NodeName == 'Multi-drug Resistance', Abundance, 0) )]

performance_mmarc_dantas_class_data[, PerformanceTotalCount := ( sum(PerformanceCount) ), by=c('SampleName', 'Pipeline')]
performance_mmarc_dantas_class_data[, TotalCountsOverall := ( sum(Abundance) ), by=c('SampleName', 'Pipeline')]
performance_mmarc_dantas_class_data[, PerformancePercent := ( 100 * sum(PerformanceCount) / sum(Abundance) ), by=c('SampleName', 'Pipeline')]
performance_mmarc_dantas_class_data <- performance_mmarc_dantas_class_data[, .SD, .SDcols=c('Pipeline','TestSet', 'SampleName', 'AnalyticTruthLabel', 'PerformanceTotalCount', 'TotalCountsOverall', 'PerformancePercent')]
setkeyv(performance_mmarc_dantas_class_data, c('SampleName', 'Pipeline'))
performance_mmarc_dantas_class_data <- unique(performance_mmarc_dantas_class_data)

dantas_sample_names <- as.character(unique(dantas_class_data$SampleName))

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
    for( p in 1:length(pipenames) ) {
        label_name <- as.character(truth_labels$AnalyticTruthLabel[truth_labels$SampleName == dantas_sample_names[s]])
        set_name <- as.character(truth_labels$TestSet[truth_labels$SampleName == dantas_sample_names[s]])
        temp <- performance_mmarc_dantas_class_data[Pipeline == pipenames[p] &
                                                        AnalyticTruthLabel == label_name &
                                                        SampleName == dantas_sample_names[s],]
        if( nrow(temp) == 0 ) {
            zero_entries <- rbind(zero_entries, data.table(Pipeline=pipenames[p],
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

performance_mmarc_dantas_class_data_plot <- performance_mmarc_dantas_class_data
levels(performance_mmarc_dantas_class_data_plot$Pipeline) <- c('Alignment', 'Meta-MARC HTS Reads', 'Meta-MARC Assembled', 'Resfams')

png('graphs/mmarc_dantas_performance_class_all_samples.png', width=1200, height=900)
g <- ggplot(performance_mmarc_dantas_class_data_plot, aes(x=AnalyticTruthLabel, y=PerformancePercent)) +
    geom_boxplot() + geom_jitter(width=0.1, height=0, size=1) +
    facet_grid(Pipeline ~ TestSet, switch='y') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x=element_text(size=20),
          strip.text.y=element_text(size=20, angle=180),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=20, angle=45, hjust=1),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.position="bottom",
          panel.margin=unit(0.1, "lines"),
          plot.title=element_text(size=30, hjust=0.5),
          legend.text=element_text(size=18),
          legend.title=element_blank()) + scale_y_continuous(labels = comma) +
    xlab('\nFunctional Metagenomic Truth Label') + ylab('Percent of Reads On-Target\n')
print(g + ggtitle(paste('Percent of Hits On-Target for each Truth Label\nby Test Set and Method\n')))
dev.off()

png('graphs/mmarc_dantas_performance_class_all_samples_notitle.png', width=1200, height=900)
print(g)
dev.off()

write.table(performance_mmarc_dantas_class_data,
            'graphs/mmarc_performance_class_all_samples_table.csv',
            sep=',', col.names = T, row.names = F)

## MMARC performance on Dantas all classes
 performance_mmarc_dantas_class_data2 <-
     dantas_class_data[, PerformanceCount := ( 
         ifelse(NodeName == AnalyticTruthLabel, Abundance, 0) )]

performance_mmarc_dantas_class_data2[, PerformanceTotalCount := ( sum(PerformanceCount) ), by=c('SampleName', 'Pipeline')]
performance_mmarc_dantas_class_data2[, TotalCountsOverall := ( sum(Abundance) ), by=c('SampleName', 'Pipeline')]
performance_mmarc_dantas_class_data2[, PerformancePercent := ( 100 * sum(PerformanceCount) / sum(Abundance) ), by=c('SampleName', 'Pipeline')]
performance_mmarc_dantas_class_data2 <- performance_mmarc_dantas_class_data2[, .SD, .SDcols=c('Pipeline','TestSet', 'SampleName', 'AnalyticTruthLabel', 'PerformanceTotalCount', 'TotalCountsOverall', 'PerformancePercent')]
setkeyv(performance_mmarc_dantas_class_data2, c('SampleName', 'Pipeline'))
performance_mmarc_dantas_class_data2 <- unique(performance_mmarc_dantas_class_data2)

dantas_sample_names <- as.character(unique(dantas_class_data$SampleName))

truth_labels2 <- performance_mmarc_dantas_class_data2[, .SD, .SDcols=c('SampleName', 'AnalyticTruthLabel', 'TestSet')]
setkey(truth_labels2, SampleName)
truth_labels2 <- unique(truth_labels2)

zero_entries = data.table(Pipeline=character(),
                          TestSet=character(),
                          SampleName=character(),
                          AnalyticTruthLabel=character(),
                          PerformanceTotalCount=numeric(),
                          TotalCountsOverall=numeric(),
                          PerformancePercent=numeric())
for( s in 1:length(dantas_sample_names) ) {
    for( p in 1:length(pipenames) ) {
        label_name <- as.character(truth_labels2$AnalyticTruthLabel[truth_labels2$SampleName == dantas_sample_names[s]])
        set_name <- as.character(truth_labels2$TestSet[truth_labels2$SampleName == dantas_sample_names[s]])
        temp <- performance_mmarc_dantas_class_data2[Pipeline == pipenames[p] &
                                                        AnalyticTruthLabel == label_name &
                                                        SampleName == dantas_sample_names[s],]
        if( nrow(temp) == 0 ) {
            zero_entries <- rbind(zero_entries, data.table(Pipeline=pipenames[p],
                                                           TestSet=set_name,
                                                           SampleName=dantas_sample_names[s],
                                                           AnalyticTruthLabel=label_name,
                                                           PerformanceTotalCount=0,
                                                           TotalCountsOverall=0,
                                                           PerformancePercent=0))
        }
    }
}

performance_mmarc_dantas_class_data2 <- rbind(performance_mmarc_dantas_class_data2, zero_entries)

performance_mmarc_dantas_class_data2_plot <- performance_mmarc_dantas_class_data2
levels(performance_mmarc_dantas_class_data2_plot$Pipeline) <- c('Alignment', 'Meta-MARC HTS Reads', 'Meta-MARC Assembled', 'Resfams')

png('graphs/mmarc_dantas_performance_class_all_samples_exact_match.png', width=1200, height=900)
g <- ggplot(performance_mmarc_dantas_class_data2_plot, aes(x=AnalyticTruthLabel, y=PerformancePercent)) +
    geom_boxplot(outlier.shape=NA) + geom_jitter(width=0.1, height=0, size=1) +
    facet_grid(Pipeline ~ TestSet, switch='y') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x=element_text(size=20),
          strip.text.y=element_text(size=20, angle=180),
          axis.text.y=element_text(size=16),
          axis.text.x=element_text(size=20, angle=45, hjust=1),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.position="bottom",
          panel.margin=unit(0.1, "lines"),
          plot.title=element_text(size=30, hjust=0.5),
          legend.text=element_text(size=18),
          legend.title=element_blank()) + scale_y_continuous(labels = comma) +
    xlab('\nFunctional Metagenomic Truth Label') + ylab('Percent of Reads On-Target\n')
print(g + ggtitle(paste('Percent of Hits On-Target for each Truth Label\nby Test Set and Pipeline\n')))
dev.off()

png('graphs/mmarc_dantas_performance_class_all_samples_exact_match_notitle.png', width=1200, height=900)
print(g)
dev.off()

write.table(performance_mmarc_dantas_class_data2,
            'graphs/mmarc_performance_class_all_samples_table_exact_match.csv',
            sep=',', col.names = T, row.names = F)


## Using the data from above, graph the number of reads found with alignment vs. MMARC
## for each truth label category
temp <- performance_mmarc_dantas_class_data
colnames(temp) <- c('Pipeline', 'TestSet', 'SampleName', 'AnalyticTruthLabel',
                    'NumberReadsOnTarget', 'TotalReadsIdentified', 'PerformancePercent')

temp <- temp[, .SD, .SDcols=c('Pipeline', 'TestSet', 'AnalyticTruthLabel',
                                                   'NumberReadsOnTarget', 'TotalReadsIdentified')]
temp <- temp[, lapply(.SD, sum), by=c('Pipeline', 'TestSet', 'AnalyticTruthLabel')]
# performance_raw <- melt(temp,
#                         id.vars=c('Pipeline', 'TestSet', 'AnalyticTruthLabel'),
#                         measure.vars=c('NumberReadsOnTarget', 'TotalReadsIdentified'),
#                         variable.name='AbundanceCategory',
#                         value.name='Abundance')
# performance_raw <- performance_raw[Pipeline != 'Resfams', ]

temp_plot <- temp

levels(temp_plot$Pipeline) <- c('Alignment', 'Meta-MARC HTS Reads', 'Meta-MARC Assembled', 'Resfams')

png('graphs/mmarc_dantas_abundance_class_all_samples.png', width=1200, height=900)
g <- ggplot(temp_plot, aes(x=AnalyticTruthLabel, fill=Pipeline)) +
    geom_bar(data=temp_plot, aes(y=TotalReadsIdentified), position='dodge', stat='identity') +
    geom_point(data=temp_plot, aes(y=NumberReadsOnTarget), position=position_dodge(width=0.9), shape=3, size=4) +
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
    xlab('\nFunctional Metagenomic Truth Label') + ylab('Number of Reads Identified\nTotal (bar) and True Positives (crosshair)\n') +
    scale_fill_brewer(palette = 'Set2') + scale_y_continuous(labels = comma)
print(g + ggtitle(paste('Total and On-Target Abundance for each Truth Label\nby Test Set and Pipeline\n')))
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

ncba_data <- abund_data[TestSet == 'NCBA' & (HMMGroup != 'groupII' | is.na(HMMGroup)), ]
colnames(ncba_metadata)[4] <- 'SampleName'
ncba_metadata$SampleName <- as.character(ncba_metadata$SampleName)
setkey(ncba_data, SampleName)
setkey(ncba_metadata, SampleName)

ncba_data <- ncba_metadata[ncba_data]
ncba_data <- ncba_data[, .SD, .SDcols=!c('BiosampleAccession', 'CattlePenID', 'TestSet', 'DantasSampleSet',
                                         'TruthLabel', 'HMMGroup')]

ncba_class_data <- ncba_data[AnnotationLevel == 'Class', ]
ncba_class_data <- ncba_class_data[, TotalAbundance :=( lapply(.SD, sum) ), by=c('Pipeline', 'SampleName', 'NodeName'), .SDcols='Abundance']
ncba_class_data <- ncba_class_data[, .SD, .SDcols=c('SampleGroup', 'SampleName', 'Pipeline', 'NodeName', 'TotalAbundance')]
setkeyv(ncba_class_data, c('NodeName', 'SampleName', 'Pipeline'))
ncba_class_data <- unique(ncba_class_data)
#norm_ncba_class_data <- ncba_class_data[, NormAbundance := ( 100 * TotalAbundance / sum(TotalAbundance) ), by=c('Pipeline', 'SampleName')]

temp <- ncba_data[AnnotationLevel == 'Class', ]
samplegroup_ncba_class_data <- temp[, GroupAbundance := ( lapply(.SD, sum) ), by=c('SampleGroup', 'Pipeline', 'NodeName'), .SDcols='Abundance']
samplegroup_ncba_class_data <- samplegroup_ncba_class_data[, .SD, .SDcols=!'Abundance']
setkeyv(samplegroup_ncba_class_data, c('SampleGroup', 'Pipeline', 'NodeName'))
samplegroup_ncba_class_data <- unique(samplegroup_ncba_class_data)

temp2 <- ncba_class_data[, .SD, .SDcols=c('Pipeline', 'NodeName', 'TotalAbundance')]
temp2 <- temp2[, NodeAbundance := ( lapply(.SD, sum)), by=c('Pipeline', 'NodeName'), .SDcols='TotalAbundance']
temp2 <- temp2[, .SD, .SDcols=!'TotalAbundance']
setkeyv(temp2, c('Pipeline', 'NodeName'))
temp2 <- unique(temp2)

# coeff = 0.01
# mmarc_threshold <- coeff * sum(temp2$NodeAbundance[temp2$Pipeline == 'MMARC'])
# align_threshold <- coeff * sum(temp2$NodeAbundance[temp2$Pipeline == 'Alignment'])
# temp2 <- temp2[!(NodeAbundance < mmarc_threshold & Pipeline == 'MMARC') & !(NodeAbundance < align_threshold & Pipeline == 'Alignment'), ]

threshold <- c()
for( x in 1:length(unique(temp2$NodeName))) {
    ltemp <- temp2[NodeName == as.character(unique(temp2$NodeName)[x])]
    if( any(ltemp$NodeAbundance > 1e+06) ) {
        threshold <- c(threshold, as.character(unique(temp2$NodeName)[x]))
    }
}
temp2 <- temp2[NodeName %in% threshold]


pipelines <- unique(as.character(temp2$Pipeline))
nodenames <- unique(as.character(temp2$NodeName))


zero_entries = data.table(Pipeline=character(),
                          NodeName=character(),
                          NodeAbundance=numeric())
for( p in 1:length(pipelines) ) {
    for( n in 1:length(nodenames) ) {
        local_temp <- temp2[Pipeline == pipelines[p] & NodeName == nodenames[n],]
        if( nrow(local_temp) == 0 ) {
            zero_entries <- rbind(zero_entries, data.table(Pipeline=pipelines[p],
                                                           NodeName=nodenames[n],
                                                           NodeAbundance=0))
        }
    }
}

temp2 <- rbind(temp2, zero_entries)

temp2_plot <- temp2

levels(temp2_plot$Pipeline) <- c('Alignment', 'Meta-MARC HTS Reads', 'Meta-MARC Assembled', 'Resfams')

png('graphs/ncba_class_all_samples.png', width=1200, height=900)
g <- ggplot(temp2_plot, aes(x=NodeName, y=NodeAbundance, fill=Pipeline)) +
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
    xlab('\nAMR Drug Class') + ylab('Number of Reads Classified\n') +
    scale_fill_discrete(drop=F) + scale_fill_brewer(palette = 'Set2') + scale_y_continuous(labels = comma)
print(g + ggtitle(paste('AMR Class Abundance by Pipeline\nPRJNA292471 Metagenomic Data Set')))
dev.off()

png('graphs/ncba_class_all_samples_notitle.png', width=1200, height=900)
print(g)
dev.off()

#########################
#########################
## Alignment coverage histograms for each test set

cov_data <- data.table(read.csv('analytic_data/alignment_cov.csv'))
cov_data$TestSet <- gsub('ncba', 'PRJNA292471', cov_data$TestSet)
cov_data$TestSet <- gsub('soil', 'Soil', cov_data$TestSet)
cov_data$TestSet <- gsub('pediatric', 'Pediatric', cov_data$TestSet)

png('graphs/ncba_alignment_all_samples.png', width=1200, height=900)
g <- ggplot(cov_data, aes(x=Coverage)) +
    geom_histogram(binwidth=0.03) + 
    geom_vline(xintercept=0.8, color='red', size=1) +
    facet_wrap(~TestSet, ncol=1) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.x=element_text(size=20),
          strip.text.y=element_text(size=20, angle=0),
          axis.text.y=element_text(size=20),
          axis.text.x=element_text(size=20),
          axis.title.x=element_text(size=24),
          axis.title.y=element_text(size=24),
          legend.position="bottom",
          panel.margin=unit(0.1, "lines"),
          plot.title=element_text(size=30, hjust=0.5),
          legend.text=element_text(size=18),
          legend.title=element_blank()) +
    xlab('\nPercent Gene Coverage') + ylab('Count\n') + scale_y_continuous(labels = comma)
print(g + ggtitle(paste('Alignment Coverage Distribution for All Genes by Data Set')))
dev.off()

png('graphs/ncba_alignment_all_samples_notitle.png', width=1200, height=900)
print(g)
dev.off()

###################################################
###################################################
## Data analysis steps ##
#########################
pdata <- performance_mmarc_dantas_class_data

## Number of samples each method classified more reads correctly pairwise comparison
table_performance_data <- data.table(Comparison=c('AlignmentVsMMARC', 'AlignmentVsMMARC-assembled',
                                                  'AlignmentVsResfams', 'MMARCVsMMARC-assembled',
                                                  'MMARCVsResfams', 'MMARCVsAlignment',
                                                  'MMARC-assembledVsMMARC', 'MMARC-assembledVsAlignment',
                                                  'MMARC-assembledVsResfams', 'ResfamsVsMMARC',
                                                  'ResfamsVsAlignment', 'ResfamsVsMMARC-assembled'),
                                     MoreTP=rep(0, 12),
                                     MoreAccurate=rep(0,12))
for( s in dantas_sample_names ) {
    temp <- pdata[SampleName == s]
    for( c in 1:nrow(table_performance_data) ) {
        members <- unlist(strsplit(table_performance_data$Comparison[c], 'Vs'))
        if( as.numeric(temp$NumberReadsOnTarget[temp$Pipeline == members[1]]) > as.numeric(temp$NumberReadsOnTarget[temp$Pipeline == members[2]]) ) {
            table_performance_data[c, MoreTP := table_performance_data[c, MoreTP] + 1]
        }
        if( as.numeric(temp$PerformancePercent[temp$Pipeline == members[1]]) > as.numeric(temp$PerformancePercent[temp$Pipeline == members[2]]) ) {
            table_performance_data[c, MoreAccurate := table_performance_data[c, MoreAccurate] + 1]
        }
    }
}

table_performance_data[, MoreTP := 100 * (MoreTP / length(dantas_sample_names)) ]
table_performance_data[, MoreAccurate := 100 * (MoreAccurate / length(dantas_sample_names)) ]

write.table(table_performance_data, 'graphs/pairwise_performance_comparison_dantas.csv', sep=',',
            row.names = F, col.names = T)


## How many times did each method perform the best with a given metric?
methods <- c('Alignment', 'MMARC', 'MMARC-assembled', 'Resfams')
truthlabels <- unique(as.character(pdata$AnalyticTruthLabel))
top_performance_data <- data.table(Pipeline=character(),
                                   TruthLabel=character(),
                                   MostTP=numeric(),
                                   MostAccurate=numeric())
for( m in methods ) {
    for( t in truthlabels ) {
        top_performance_data <- rbind(top_performance_data, data.table(Pipeline=m,
                                                                       TruthLabel=t,
                                                                       MostTP=0,
                                                                       MostAccurate=0))
    }
}
for( s in dantas_sample_names ) {
    temp <- pdata[SampleName == s]
    for( t in unique(as.character(temp$AnalyticTruthLabel)) ) {
        temp2 <- temp[AnalyticTruthLabel == t]
        topTP <- as.character(temp2$Pipeline[temp2$NumberReadsOnTarget == max(temp2$NumberReadsOnTarget)])
        topAcc <- as.character(temp2$Pipeline[temp2$PerformancePercent == max(temp2$PerformancePercent)])
        top_performance_data[Pipeline == topTP & TruthLabel == t,
                             MostTP := top_performance_data[Pipeline == topTP & TruthLabel == t, MostTP] + 1]
        top_performance_data[Pipeline == topAcc & TruthLabel == t,
                             MostAccurate := top_performance_data[Pipeline == topAcc & TruthLabel == t,
                                                                  MostAccurate] + 1]
    }
}

top_performance_data[, MostTP := 100 * (MostTP / length(dantas_sample_names)) ]
top_performance_data[, MostAccurate := 100 * (MostAccurate / length(dantas_sample_names)) ]

mtop <- melt(top_performance_data, id=1:2, measure=c('MostTP', 'MostAccurate'),
             variable.name='MeasureVar', value.name='SampleCount')

mtop$MeasureVar <- gsub('MostAccurate', 'Highest Proportion Classified Correctly',
                        as.character(mtop$MeasureVar))
mtop$MeasureVar <- gsub('MostTP', 'Highest Number of True Positives',
                        as.character(mtop$MeasureVar))

png('graphs/dantas_top_pipeline_tpaccuracy_comparison.png', width=1200, height=900)
g <- ggplot(mtop, aes(x=TruthLabel, y=SampleCount, fill=Pipeline)) +
    geom_bar(stat='identity', position='dodge') + 
    facet_wrap(~MeasureVar, ncol=1) +
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
    xlab('\nFunctional Metagenomic Truth Label') + ylab('Number of Soil and Pediatric Samples\n') +
    scale_fill_brewer(palette = 'Set2', drop=F)
print(g + ggtitle(paste('Samplewise Pipeline Performance by Category\nTrue Positives and Accuracy')))
dev.off()

png('graphs/dantas_top_pipeline_tpaccuracy_comparison_notitle.png', width=1200, height=900)
print(g)
dev.off()













