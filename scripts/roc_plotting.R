library(data.table)
library(ggplot2)
library(gridExtra)

calculate_trap_auc_pr <- function(precision, recall) {
  if(recall[1] > recall[length(recall)]) {
    precision <- rev(precision)
    recall <- rev(recall)
  }
  cumulative_auc <- 0
  for(i in 2:length(precision)) {
    width <- abs(recall[i] - recall[i-1])
    lower_point <- ifelse(precision[i] > precision[i-1], precision[i-1], precision[i])
    rect_area <- width * lower_point
    tri_area <- width * abs(precision[i-2] - precision[i]) / 2
    total_area <- rect_area + tri_area
    if(length(total_area) > 0) {
      cumulative_auc <- cumulative_auc + total_area
    } else {
      next
    }
  }
  return(cumulative_auc)
}

calculate_trap_auc_roc <- function(sens, invspec) {
  if(sens[1] > sens[length(sens)]) {
    sens <- rev(sens)
    invspec <- rev(invspec)
  }
  cumulative_auc <- 0
  for(i in 2:length(sens)) {
    width <- abs(invspec[i] - invspec[i-1])
    lower_point <- ifelse(sens[i] > sens[i-1], sens[i-1], 
                          sens[i])
    rect_area <- width * lower_point
    tri_area <- width * abs(sens[i-2] - sens[i]) / 2
    total_area <- rect_area + tri_area
    if(length(total_area) > 0) {
      cumulative_auc <- cumulative_auc + total_area
    } else {
      next
    }
  }
  return(cumulative_auc)
}


add_corner_cases <- function(dt) {
  new_row <- dt[1]
  new_row[['Evalue']] <- Inf
  new_row[['Sensitivity']] <- 1
  new_row[['OneMinusSpecificity']] <- 1
  new_row[['Specificity']] <- 0
  new_row[['Precision']] <- 0
  new_row[['Recall']] <- 1
  
  return_dt <- rbind(dt, new_row)
  
  new_row[['Evalue']] <- 0
  new_row[['Sensitivity']] <- 0
  new_row[['OneMinusSpecificity']] <- 0
  new_row[['Specificity']] <- 1
  new_row[['Precision']] <- 1
  new_row[['Recall']] <- 0
  
  return_dt <- rbind(new_row, return_dt)
  
  return(return_dt)
}

#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

roc <- data.table(read.csv('/mnt/manuscripts/CommunicationsBiology_HMM_Publication/meta-marc-publication/analytic_data/roc_data.csv', header=T))

roc[, Sensitivity :=( TP / (TP + FN) )]
roc[, OneMinusSpecificity :=( 1 - (TN / (TN + FP)) )]
roc[, Specificity :=( TN / (TN + FP) )]
roc[, Precision :=( TP / (TP + FP) )]
roc[, Recall :=( TP / (TP + FN) )]

roc[is.nan(Precision)][['Precision']] <- 1

full_roc <- roc[, add_corner_cases(.SD), by=c('NodeName', 'Level')]

meanroc <- roc[, lapply(.SD, mean), .SDcols=c("Sensitivity", "OneMinusSpecificity", "Precision", "Recall"), by=c("Level", "Evalue")]
medianroc <- roc[, lapply(.SD, median), .SDcols=c("Sensitivity", "OneMinusSpecificity", "Precision", "Recall"), by=c("Level", "Evalue")]

meanroc <- rbind(data.table(
  Level=c('Model', 'Group', 'Mechanism', 'Class'),
  Evalue=c(0, 0, 0, 0),
  Sensitivity=c(1, 1, 1, 1),
  OneMinusSpecificity=c(1, 1, 1, 1),
  Precision=c(0, 0, 0, 0),
  Recall=c(1, 1, 1, 1)
), meanroc)

medianroc <- rbind(data.table(
  Level=c('Model', 'Group', 'Mechanism', 'Class'),
  Evalue=c(0, 0, 0, 0),
  Sensitivity=c(1, 1, 1, 1),
  OneMinusSpecificity=c(1, 1, 1, 1),
  Precision=c(0, 0, 0, 0),
  Recall=c(1, 1, 1, 1)
), medianroc)

meanroc_plot <- ggplot(meanroc, aes(x=OneMinusSpecificity, y=Sensitivity, group=Level, color=Level)) +
  geom_line(size=1.8) + ylim(0, 1) + xlim(0, 1) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=20),
    axis.title.x=element_text(size=24),
    axis.title.y=element_text(size=24),
    legend.position="bottom",
    panel.margin=unit(0.1, "lines"),
    plot.title=element_text(size=30, hjust=0.5),
    legend.text=element_text(size=18),
    legend.title=element_blank()
  ) + xlab('1 - Specificity')
png('/mnt/manuscripts/CommunicationsBiology_HMM_Publication/meta-marc-publication/graphs/mean_roc.png', width=1200, height=900)
print(meanroc_plot)
dev.off()

mean_precisionrecall_plot <- ggplot(meanroc, aes(x=Recall, y=Precision, group=Level, color=Level)) +
  geom_line(size=1.8) + xlim(c(0, 1)) + ylim(c(0, 1)) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y=element_text(size=20),
    axis.text.x=element_text(size=20),
    axis.title.x=element_text(size=24),
    axis.title.y=element_text(size=24),
    legend.position="none",
    panel.margin=unit(0.1, "lines"),
    plot.title=element_text(size=30, hjust=0.5),
    legend.text=element_text(size=18),
    legend.title=element_blank()
  )
png('/mnt/manuscripts/CommunicationsBiology_HMM_Publication/meta-marc-publication/graphs/mean_pr.png', width=1200, height=900)
print(mean_precisionrecall_plot)
dev.off()

roc_pr_legend <- g_legend(meanroc_plot)

png('/mnt/manuscripts/CommunicationsBiology_HMM_Publication/meta-marc-publication/graphs/mean_roc_pr.png', width=1200, height=900)
roc_pr_plot <- grid.arrange(arrangeGrob(
  meanroc_plot + theme(legend.position="none"),
  mean_precisionrecall_plot,
  nrow=1
), roc_pr_legend, nrow=2, heights=c(8, 1))
print(roc_pr_plot)
dev.off()

pr_auc <- full_roc[, calculate_trap_auc_pr(Precision, Recall), by=c('Level', 'NodeName')]
roc_auc <- full_roc[, calculate_trap_auc_roc(Sensitivity, OneMinusSpecificity), by=c('Level', 'NodeName')]

print(tapply(pr_auc$V1, pr_auc$Level, summary))
print(tapply(roc_auc$V1, roc_auc$Level, summary))
