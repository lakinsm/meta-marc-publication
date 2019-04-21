library(ggplot2)
library(data.table)
library(tidyr)

df = data.table(read.csv('/mnt/datasets/hmm_simulations/2019_publication/analytic_data/experiment5_results.csv',
                         header=F))

colnames(df) = c("Method", "GeneID", "FractionRandomized", 
                 "GeneCopyNumber", "ContiguousRandomization",
                 "Level", "NodeName", "Classified", "Total")
df[, GeneID := NULL]

model = df[Level == 3]
model[, ClassifiedCollapsed := sum(Classified), 
                  by=c('Method', 'FractionRandomized', 
                       'GeneCopyNumber', 'ContiguousRandomization', 
                       'NodeName')]
model[, TotalCollapsed := sum(Total), 
      by=c('Method', 'FractionRandomized', 
           'GeneCopyNumber', 'ContiguousRandomization', 
           'NodeName')]
model[, PercentRecovered := (ClassifiedCollapsed / TotalCollapsed)]
model[, Level := rep('Model', nrow(model))]

group = df[Level == 2]
group[, ClassifiedCollapsed := sum(Classified), 
                  by=c('Method', 'FractionRandomized', 
                       'GeneCopyNumber', 'ContiguousRandomization', 
                       'NodeName')]
group[, TotalCollapsed := sum(Total), 
      by=c('Method', 'FractionRandomized', 
           'GeneCopyNumber', 'ContiguousRandomization', 
           'NodeName')]
group[, PercentRecovered := (ClassifiedCollapsed / TotalCollapsed)]
group[, Level := rep('Group', nrow(group))]

mech = df[Level == 1]
mech[, ClassifiedCollapsed := sum(Classified), 
                by=c('Method', 'FractionRandomized', 
                     'GeneCopyNumber', 'ContiguousRandomization', 
                     'NodeName')]
mech[, TotalCollapsed := sum(Total), 
     by=c('Method', 'FractionRandomized', 
          'GeneCopyNumber', 'ContiguousRandomization', 
          'NodeName')]
mech[, PercentRecovered := (ClassifiedCollapsed / TotalCollapsed)]
mech[, Level := rep('Mechanism', nrow(mech))]

class = df[Level == 0]
class[, ClassifiedCollapsed := sum(Classified), 
                  by=c('Method', 'FractionRandomized', 
                       'GeneCopyNumber', 'ContiguousRandomization', 
                       'NodeName')]
class[, TotalCollapsed := sum(Total), 
      by=c('Method', 'FractionRandomized', 
           'GeneCopyNumber', 'ContiguousRandomization', 
           'NodeName')]
class[, PercentRecovered := (ClassifiedCollapsed / TotalCollapsed)]
class[, Level := rep('Class', nrow(class))]

all_collapsed = rbind(model, group, mech, class)

method_transl = c("alignment"="Alignment",
                  "mmarc_hts"="Meta-MARC HTS Reads",
                  "mmarc_assembly"="Meta-MARC Assembly",
                  "resfams"="Resfams")
all_collapsed$Method = as.character(method_transl[as.character(all_collapsed$Method)])
all_collapsed$Method = factor(all_collapsed$Method, 
                              levels=c("Alignment", 
                                       "Meta-MARC HTS Reads",
                                       "Meta-MARC Assembly",
                                       "Resfams"),
                              ordered=T)

random_transl = c("true"="Noncontiguous",
                  "false"="Contiguous")
all_collapsed$ContiguousRandomization = random_transl[all_collapsed$ContiguousRandomization]
all_collapsed$ContiguousRandomization = factor(all_collapsed$ContiguousRandomization,
                                               levels=c("Contiguous",
                                                        "Noncontiguous"),
                                               ordered=T)

all_collapsed$Level = factor(all_collapsed$Level,
                             levels=c("Class", "Mechanism", "Group", "Model"),
                             ordered=T)

all_collapsed[, `:=` (Classified = NULL,
                      ClassifiedCollapsed = NULL,
                      Total = NULL,
                      TotalCollapsed = NULL,
                      NodeName = NULL)]

setkey(all_collapsed, Method, FractionRandomized, GeneCopyNumber,
       ContiguousRandomization, Level)
df2 = all_collapsed[CJ(Method, FractionRandomized,
                 GeneCopyNumber, ContiguousRandomization,
                 Level, unique=T)]
df2[is.na(PercentRecovered)]$PercentRecovered = 0

df2 = df2[!(Method == "Resfams" & Level %in% c("Mechanism",
                                               "Group",
                                               "Model"))]
df2 = df2[!(Method == "Alignment" & Level == "Model")]

gene_transl = c("5"="5 Gene Copies",
  "25" = "25 Gene Copies",
  "50" = "50 Gene Copies")
df2$GeneCopyNumber = as.character(gene_transl[as.character(df2$GeneCopyNumber)])
df2$GeneCopyNumber = factor(as.character(df2$GeneCopyNumber),
                                         levels=c("5 Gene Copies",
                                                  "25 Gene Copies",
                                                  "50 Gene Copies"),
                                         ordered=T)


png('/mnt/datasets/hmm_simulations/2019_publication/graphs/experiment5_levels_randomization.png', width=1200, height=900)

plotting = df2[GeneCopyNumber == "5 Gene Copies"]
# plotting = df2[GeneCopyNumber == 5]
all_levels = ggplot(plotting, aes(x=FractionRandomized, y=PercentRecovered, color=Method)) +
  facet_wrap(Level~ContiguousRandomization, nrow=4, ncol=2) +
  theme(strip.text.x=element_text(size=20),
        strip.text.y=element_text(size=20, angle=0),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=20, angle=45, hjust=1),
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        legend.position="bottom",
        legend.text=element_text(size=18),
        legend.title=element_blank()) +
  scale_y_continuous(limits = c(-0.3, 1), breaks=seq(0, 1, .2),
                     labels=c(seq(0, 1, .2))) +
 geom_smooth(method='loess', size=1.5) +
  ylab("Simulated Genes Classified Correctly (%)\n") +
  xlab("Gene Randomization (%)")
print(all_levels)
dev.off()

plotting2 = df2[Level=="Class"]
png('/mnt/datasets/hmm_simulations/2019_publication/graphs/experiment5_copynumber_randomization.png', width=1200, height=900)
copynumber = ggplot(plotting2, aes(x=FractionRandomized, y=PercentRecovered, color=Method)) +
  facet_wrap(GeneCopyNumber~ContiguousRandomization, nrow=3, ncol=2) +
  theme(strip.text.x=element_text(size=20),
        strip.text.y=element_text(size=20, angle=0),
        axis.text.y=element_text(size=20),
        axis.text.x=element_text(size=20, angle=45, hjust=1),
        axis.title.x=element_text(size=24),
        axis.title.y=element_text(size=24),
        legend.position="bottom",
        legend.text=element_text(size=18),
        legend.title=element_blank()) +
  scale_y_continuous(limits = c(-0.3, 1), breaks=seq(0, 1, .2),
                     labels=c(seq(0, 1, .2))) +
  geom_smooth(method='loess', size=1.5) +
  ylab("Simulated Genes Classified Correctly (%)\n") +
  xlab("Gene Randomization (%)")
print(copynumber)
  
dev.off()

