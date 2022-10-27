library(ggplot2)
library(ggforce) # On Ubuntu 22.04 LTS: sudo apt -y install libfontconfig1-dev
library(ggpubr)
library(tidyverse)

options(scipen = 0)

stats_df <- read.table('~/git/BiWFA-paper/evaluation/data/statistics.simulations.tsv.gz', sep = '\t', header = F)
names(stats_df) <- c('error', 'seq.length', 'seq', 'mode', 'replicate', 'value', 'statistic')
# Compute averages for time and maximums for memory
stats_time_df <- stats_df %>%
  filter(statistic == 'time_ns') %>%
  group_by(error, seq.length, seq, mode, statistic) %>%
  dplyr::summarize(value = mean(value, na.rm=TRUE), num.replicates = n())
stats_memory_df <- stats_df %>%
  filter(statistic == 'memory_kb') %>%
  group_by(error, seq.length, seq, mode, statistic) %>%
  dplyr::summarize(value = max(value), num.replicates = n())
stats_df <- bind_rows(stats_time_df, stats_memory_df)
rm(stats_time_df, stats_memory_df)

scores_df <- read.table('~/git/BiWFA-paper/evaluation/data/scores.simulations.tsv.gz', sep = '\t', header = F)
names(scores_df) <- c('error', 'seq.length', 'seq', 'mode', 'score')

statsWithMetadata_df <- merge(
  stats_df,
  scores_df,
  all.x=T,
  by.x = c('error', 'seq.length', 'seq', 'mode'), by.y = c('error', 'seq.length', 'seq', 'mode')
)
statsWithMetadata_df$score[is.na(statsWithMetadata_df$score)] <- 1

# Force column's type
statsWithMetadata_df$error <- as.factor(statsWithMetadata_df$error)
statsWithMetadata_df$mode <- as.factor(statsWithMetadata_df$mode)

levels(statsWithMetadata_df$error)[match("0.001",levels(statsWithMetadata_df$error))] <- "Error rate = 0.1%"
levels(statsWithMetadata_df$error)[match("0.01",levels(statsWithMetadata_df$error))] <- "Error rate = 1%"
levels(statsWithMetadata_df$error)[match("0.05",levels(statsWithMetadata_df$error))] <- "Error rate = 5%"
levels(statsWithMetadata_df$error)[match("0.1",levels(statsWithMetadata_df$error))] <- "Error rate = 10%"
levels(statsWithMetadata_df$error)[match("0.2",levels(statsWithMetadata_df$error))] <- "Error rate = 20%"
levels(statsWithMetadata_df$error)[match("0.3",levels(statsWithMetadata_df$error))] <- "Error rate = 30%"

# Rename algorithms and change their order
levels(statsWithMetadata_df$mode)[match("wfa-ultralow",levels(statsWithMetadata_df$mode))] <- "BiWFA"
levels(statsWithMetadata_df$mode)[match("wfa-high",levels(statsWithMetadata_df$mode))] <- "WFA-high"

statsWithMetadata_df$mode <- factor(
  statsWithMetadata_df$mode,
  levels = c("WFA-high", "BiWFA")
)

# Plots
# Memory
px <- ggplot(
  statsWithMetadata_df[statsWithMetadata_df$statistic != 'time_ns', ],
  aes(x = seq.length, y = value / 1024, color = mode, alpha=I(0.5))
) +
  geom_point() +
  facet_wrap (
    ~error,
    ncol = 2
  ) +
  ggtitle('Maximum memory consumption') +
  guides(color=guide_legend(title="Algorithm")) +
  xlab("Length") + ylab("") +  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "none",
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    
    axis.title=element_text(size=18),
    
    #axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1, size=15),
    axis.title.x = element_blank(),
    axis.text.x=element_blank(), #remove x axis labels
    #axis.ticks.x=element_blank(), #remove x axis ticks
    
    axis.text.y=element_text(size=15),
    #axis.text.y=element_blank(),  #remove y axis labels
    #axis.ticks.y=element_blank()  #remove y axis ticks
    
    strip.text = element_text(size=15)
  ) +
  scale_y_continuous(
    trans='log10',
    #labels = scales::comma,
    limits=c(1, 10000),
    n.breaks = 6,
    labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "NA" = "")
  ) +
  scale_x_continuous(
    #trans='log10'
    #labels = scales::comma,
    limits=c(250, 20000),
    n.breaks = 12
  ) + 
  scale_color_manual(values=c("#39B600", "#FF62BC"))
px

# Time
py <- ggplot(
  statsWithMetadata_df[statsWithMetadata_df$statistic == 'time_ns', ],
  aes(x = seq.length, y = value / 1000, color = mode, alpha=I(0.5))
) +
  geom_point() +
  facet_wrap (
    ~error,
    ncol = 2
  ) +
  ggtitle('Average execution time') +
  guides(color=guide_legend(title="Algorithm")) +
  xlab("Length") + ylab("Microseconds") +  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "none",
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    
    axis.title=element_text(size=18),
    
    axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1, size=15),
    #axis.title.x = element_blank(),
    #axis.text.x=element_blank(), #remove x axis labels
    #axis.ticks.x=element_blank(), #remove x axis ticks
    
    axis.text.y=element_text(size=15),
    #axis.text.y=element_blank(),  #remove y axis labels
    #axis.ticks.y=element_blank()  #remove y axis ticks
    
    strip.text = element_text(size=15)
  ) +
  scale_y_continuous(
    trans='log10',
    #labels = scales::comma,
    #limits=c(1, 1000000),
    n.breaks = 6,
    #labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
  ) +
  scale_x_continuous(
    #trans='log10'
    #labels = scales::comma,
    limits=c(250, 20000),
    n.breaks = 12
  ) + 
  scale_color_manual(values=c("#39B600", "#FF62BC"))
py

# Plot both and save the image
pxy <- ggpubr::ggarrange(
  px, py,
  align='v',
  labels=c('A', 'B'),
  heights=c(1, 1.2),
  font.label = list(size = 18),
  legend = "top", # legend position,
  common.legend = T,
  nrow = 2
) + 
  theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) # to avoid cutting labels
pxy
ggsave(plot = pxy, paste0('FigureS3', '.pdf'), width = 30, height = 40, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
ggsave(plot = pxy, paste0('FigureS3', '.png'), width = 30, height = 40, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)


if (FALSE) {
  a <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'time_s', ] #%>%
  #dplyr::filter(seq.length <= 10000) #%>%
  #dplyr::filter(grepl("_n4",seq))
  px1 <- ggplot(
    a,
    aes(x = -score, y = value + 0.001, color = mode)
  ) +
    geom_point() +
    ggtitle('Execution time') +
    xlab("Score") + ylab("Seconds") +  theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1)
    ) +
    scale_y_continuous(
      #trans='log10'
      #labels = scales::comma,
      #limits=c(0, 0.01),
      #n.breaks = 6,
      #labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
    ) +
    scale_x_continuous(
      #trans='log10'
      #labels = scales::comma,
      #limits=c(1, 1000),
      n.breaks = 12
    )# +
  #coord_equal()+
  #theme(aspect.ratio=1)
  px1
  ggsave(plot = px1, paste0('ScoreVsTime', '.pdf'), width = 20, height = 20, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
  
  
  
  
  
  a <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'time_s', ]
  a$seq.length <- as.factor(a$seq.length)
  ggplot(
    a,
    aes(x = seq.length, y = value + 0.001)
  ) +
    geom_boxplot(aes(fill=mode, color=error))
  
  
  b <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'memory_kb', ]
  py <- ggplot(
    b,
    aes(x = -score, y = value / 1024, color = mode, alpha=I(1/5))
  ) +
    geom_point() +
    ggtitle('Memory consumption') +
    xlab("Score") + ylab("MB") +  theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1)
    ) +
    scale_y_continuous(
      #trans='log10',
      #labels = scales::comma,
      #limits=c(1000, 500000),
      #n.breaks = 6#,
      #labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
    ) +
    scale_x_continuous(
      #trans='log10',
      #labels = scales::comma,
      #limits=c(1000, 500000),
      #n.breaks = 6
    ) + 
    coord_equal()+
    theme(aspect.ratio=1)
  py
  ggsave(plot = py, paste0('ScoreVsMemory', '.pdf'), width = 20, height = 20, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
}
