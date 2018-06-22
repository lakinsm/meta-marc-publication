library(ggplot2)
library(data.table)


bench_data <- data.table(
  Method=c('Alignment', 'Meta-MARC HTS', 'Meta-MARC Assembly', 'Resfams'),
  WallClockTime=c(2.51, 267.04+82.32+506.09, 469.34+133.6+882.2, 24.3),
  MemUse=c(23076, 36160+41736+44860, 305396+258088+314428, 59652)
)

bench_data[, MemUseMB := ( MemUse / 1000 )]
bench_data[,MemUse:=NULL]

mbdata <- melt(bench_data, id=1)
colnames(mbdata) <- c('Method', 'Resource', 'Value')

transl_table = c('WallClockTime'='CPU Use',
                 'MemUseMB'='Memory Use')

mbdata[['Resource']] <- transl_table[mbdata[['Resource']]]

g <- ggplot(mbdata, aes(x=Method, y=Value)) +
  facet_wrap(~Resource, ncol=1, scales = 'free_y') + geom_bar(stat='identity') +
  ylab('Peak Memory Use (MB) and Wall Clock Time (seconds)\n') +
  xlab('\nMethod') +
  theme(
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=20),
    axis.title.x=element_text(size=24),
    axis.title.y=element_text(size=24),
    strip.text=element_text(size=22),
    panel.margin=unit(0.1, "lines")
  )
  
png('/mnt/manuscripts/CommunicationsBiology_HMM_Publication/meta-marc-publication/graphs/benchmarking.png', width=1200, height=900)
print(g)
dev.off()
