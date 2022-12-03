library(ggplot2)
library(ggforce) # On Ubuntu 22.04 LTS: sudo apt -y install libfontconfig1-dev
library(ggpubr)
library(tidyverse)

options(scipen = 9)

stats_df <- read.table('~/git/BiWFA-paper/evaluation/data/statistics_all.tsv.gz', sep = '\t', header = F)
names(stats_df) <- c('set', 'seq', 'mode', 'replicate', 'value', 'statistic', 'error.rate')

# Compute averages for time and maximums for memory
#stats_time_df <- stats_df %>%
#  filter(statistic != 'memory_kb') %>% #time_s or time_ns
#  group_by(set, seq, mode, statistic) %>%
#  dplyr::summarize(value = mean(value, na.rm=TRUE), error.rate = mean(error.rate, na.rm=TRUE), num.replicates = n())
#stats_memory_df <- stats_df %>%
#  filter(statistic == 'memory_kb') %>%
#  group_by(set, seq, mode, statistic) %>%
#  dplyr::summarize(value = max(value), error.rate = mean(error.rate, na.rm=TRUE), num.replicates = n())
#stats_df <- bind_rows(stats_time_df, stats_memory_df)
#rm(stats_time_df, stats_memory_df)

lengths_df <- read.table('~/git/BiWFA-paper/evaluation/data/lengths_all.tsv.gz', sep = '\t', header = F)
names(lengths_df) <- c('set', 'seq', 'query.length', 'target.length')
scores_df <- read.table('~/git/BiWFA-paper/evaluation/data/scores_all.tsv.gz', sep = '\t', header = T)
names(scores_df) <- c('set', 'seq', 'mode', 'score')


statsWithMetadata_df <- merge(
  stats_df,
  lengths_df,
  all.x=T,
  by.x = c('set', 'seq'), by.y = c('set', 'seq')
)
statsWithMetadata_df <- merge(
  statsWithMetadata_df,
  scores_df,
  all.x=T,
  by.x = c('set', 'seq', 'mode'), by.y = c('set', 'seq', 'mode')
)
statsWithMetadata_df$score[is.na(statsWithMetadata_df$score)] <- 1
statsWithMetadata_df$statistic[statsWithMetadata_df$statistic == ''] <- NaN

# Force column's type
statsWithMetadata_df$set <- as.factor(statsWithMetadata_df$set)
statsWithMetadata_df$mode <- as.factor(statsWithMetadata_df$mode)

# Rename datasets
levels(statsWithMetadata_df$set)[match("ont_regions",levels(statsWithMetadata_df$set))] <- "ONT PromethION reads vs CHM13 v1.1"
levels(statsWithMetadata_df$set)[match("ONT",levels(statsWithMetadata_df$set))] <- "ONT Ultra Long"
statsWithMetadata_df$set <- factor(
  statsWithMetadata_df$set,
  levels = c("ONT PromethION reads vs CHM13 v1.1", "ONT Ultra Long")
)

# Rename algorithms and change their order
levels(statsWithMetadata_df$mode)[match("edlib",levels(statsWithMetadata_df$mode))] <- "edlib (edit distance)"
levels(statsWithMetadata_df$mode)[match("bitpal-scored",levels(statsWithMetadata_df$mode))] <- "bitpal (score only)"
levels(statsWithMetadata_df$mode)[match("biwfa",levels(statsWithMetadata_df$mode))] <- "BiWFA"
levels(statsWithMetadata_df$mode)[match("wfa-high",levels(statsWithMetadata_df$mode))] <- "WFA-high"
levels(statsWithMetadata_df$mode)[match("wfa-low",levels(statsWithMetadata_df$mode))] <- "WFA-low"
levels(statsWithMetadata_df$mode)[match("wfa-med",levels(statsWithMetadata_df$mode))] <- "WFA-med"
levels(statsWithMetadata_df$mode)[match("wfalm-lowmem",levels(statsWithMetadata_df$mode))] <- "wfalm-low"
levels(statsWithMetadata_df$mode)[match("wfalm-rec",levels(statsWithMetadata_df$mode))] <- "wfalm-recursive"

statsWithMetadata_df$mode <- factor(
  statsWithMetadata_df$mode,
  levels = c("edlib (edit distance)", "bitpal (score only)", "ksw2-extz2-sse", "WFA-high", "WFA-med", "WFA-low",  "wfalm", "wfalm-low", "wfalm-recursive", "BiWFA")
)

################################################################################
# Long sequences
z <- statsWithMetadata_df[statsWithMetadata_df$set %in% c(
  "ONT PromethION reads vs CHM13 v1.1", "ONT Ultra Long"
),]
z <- z[
  (z$set == 'ONT PromethION reads vs CHM13 v1.1' & z$query.length > 10000 & z$target.length > 10000) |
   z$set == 'ONT Ultra Long' & z$query.length > 500000 & z$target.length > 500000,]
levels(z$set)[match("ONT PromethION reads vs CHM13 v1.1",levels(z$set))] <- "ONT PromethION reads vs CHM13 v1.1 > 10 kbps"
levels(z$set)[match("ONT Ultra Long",levels(z$set))] <- "ONT Ultra Long > 500 kbps"

# Plot memory use
x <- z[z$statistic == 'memory_kb' & !is.nan(z$value), ]
px <- ggplot(x, aes(x = mode, y = value / 1024, fill=mode)) +
  geom_boxplot() +
  scale_y_continuous(
    trans='log10',
    #labels = scales::comma,
    limits=c(1, 1000000),
    n.breaks = 6,
    labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
  ) +
  facet_wrap (
    ~set,
    ncol = 2
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "none",
    
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
  guides(fill=guide_legend(title="Algorithm")) +
  ggtitle('Maximum memory consumption') +
  xlab("Algorithm") + ylab("") + geom_vline(xintercept=2.5)
px

# Plot runtime
y <- z[z$statistic == 'time_ms' & !is.nan(z$value), ]
py <- ggplot(y, aes(x = mode, y = value / 1000, fill=mode)) +
  geom_boxplot() +
  scale_y_continuous(
    trans='log10',
    labels = scales::comma
  ) +
  facet_wrap (
    ~set,
    ncol = 2
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "none",
    
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
  guides(fill=guide_legend(title="Algorithm")) +
  ggtitle('Execution time') +
  xlab("Algorithm") + ylab("Seconds") + geom_vline(xintercept=2.5)
py

# Plot both and save the image
library(ggpubr)
pxy <- ggpubr::ggarrange(
  px, py,
  align='v',
  labels=c('A', 'B'),
  heights=c(1, 1.4),
  font.label = list(size = 18),
  legend = "none", # legend position,
  common.legend = T,
  nrow = 2
) + 
  theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) # to avoid cutting labels
pxy
ggsave(plot = pxy, paste0('Figure2', '.pdf'), width = 30, height = 20, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
ggsave(plot = pxy, paste0('Figure2', '.png'), width = 30, height = 20, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
################################################################################

################################################################################
# Sequences <= 10kbps
z <- statsWithMetadata_df[statsWithMetadata_df$set %in% c(
  "ONT PromethION reads vs CHM13 v1.1", "ONT Ultra Long"
),]
z <- z[z$query.length <= 10000 & z$target.length <= 10000,]
levels(z$set)[match("ONT PromethION reads vs CHM13 v1.1",levels(z$set))] <- "ONT PromethION reads vs CHM13 v1.1 <= 10 kbps"
levels(z$set)[match("ONT Ultra Long",levels(z$set))] <- "ONT Ultra Long <= 10 kbps"

# Plot memory use
x <- z[z$statistic == 'memory_kb' & !is.nan(z$value), ]
px <- ggplot(x, aes(x = mode, y = value / 1024, fill=mode)) +
  geom_boxplot() +
  scale_y_continuous(
    trans='log10',
    #labels = scales::comma,
    limits=c(1, 1000),
    n.breaks = 6,
    labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "NA" = "")
  ) +
  facet_wrap (
    ~set,
    ncol = 2
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "none",
    
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
  guides(fill=guide_legend(title="Algorithm")) +
  ggtitle('Maximum memory consumption') +
  xlab("Algorithm") + ylab("") + geom_vline(xintercept=2.5)
px

# Plot runtime
y <- z[z$statistic == 'time_ms' & !is.nan(z$value), ]
py <- ggplot(y, aes(x = mode, y = value, fill=mode)) +
  geom_boxplot() +
  scale_y_continuous(
    trans='log10',
    labels = scales::comma
  ) +
  facet_wrap (
    ~set,
    ncol = 2
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "none",
    
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
  guides(fill=guide_legend(title="Algorithm")) +
  ggtitle('Average execution time') +
  xlab("Algorithm") + ylab("Milliseconds") + geom_vline(xintercept=2.5)
py

# Plot both and save the image
pxy <- ggpubr::ggarrange(
  px, py,
  align='v',
  labels=c('A', 'B'),
  heights=c(1, 1.4),
  font.label = list(size = 18),
  legend = "none", # legend position,
  common.legend = T,
  nrow = 2
) + 
  theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) # to avoid cutting labels
pxy
ggsave(plot = pxy, paste0('FigureS1', '.pdf'), width = 30, height = 20, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
ggsave(plot = pxy, paste0('FigureS1', '.png'), width = 30, height = 20, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
################################################################################


################################################################################
# WFA-high vs BiWFA

# Memory
z <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'memory_kb' & !is.nan(statsWithMetadata_df$value), ]
z <- z[z$mode %in% c('WFA-high', 'BiWFA'),]# & z$num.replicates == 100,]
z <- z[z$set %in% c('ONT Ultra Long'),]
z $error.rate <- as.numeric(z$error.rate)
px <- ggplot(
  z,
  aes(x = (query.length + target.length) / 2, y = value / 1024, color = error.rate, shape=mode, alpha=I(1/2))
) +
  geom_point(size=2) +
  #facet_wrap (
  #  ~set,
  #  ncol = 2
  #) +
  ggtitle('Maximum memory consumption') +
  xlab("Length (base pairs)") + ylab("") +  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "top",
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
    
    strip.text = element_text(size=15),
  ) + guides(color=guide_legend(title="Error rate (%)"), shape=guide_legend(title="Algorithm")) +
  scale_y_continuous(
    trans='log10',
    #labels = scales::comma,
    limits=c(1, 10000),
    n.breaks = 6,
    labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "NA" = "")
  ) +
  scale_x_continuous(
    trans='log10',
    #labels = scales::comma,
    #limits=c(min((z$query.length + z$target.length)/2), 500000),
    n.breaks = 12
  ) + 
  scale_color_gradient(low = "green", high = "red") #+ scale_color_manual(values=c("#39B600", "#FF62BC"))
px
if(FALSE) {
  library(dplyr)
  zz <- z %>%
    group_by(set, seq, statistic, query.length, target.length, score) %>%
    summarise_at(vars(value), list(name = diff))
  ggplot(
    zz,
    aes(x = (query.length + target.length) / 2, y = name, alpha=I(1/3))
  ) +
    geom_point() +
    ggtitle('Memory consumption\nWFA-high - BiWFA') +
    xlab("Length") + ylab("") +  theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size=18),
      legend.position = "none",
      
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
      limits=c(1, 1000000),
      n.breaks = 6,
      labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
    ) +
    scale_x_continuous(
      #trans='log10',
      #labels = scales::comma,
      limits=c(min((z$query.length + z$target.length)/2), 500000),
      n.breaks = 12
    )
}

# Time
z <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'time_ms' & !is.nan(statsWithMetadata_df$value), ]
z <- z[z$mode %in% c('WFA-high', 'BiWFA'),]# & z$num.replicates == 100,]
z $error.rate <- as.numeric(z$error.rate)
py <- ggplot(
  z,
  aes(
    x = (query.length + target.length) / 2,
    y = value,
    color = mode, shape = mode, alpha=I(1/2))
) +
  geom_point(size=2) +
  facet_wrap (
    ~set,
    ncol = 2
  ) +
  ggtitle('Average execution time') +
  xlab("Length (base pairs)") + ylab("Milliseconds") +  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "top",
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
    
    strip.text = element_text(size=15),
  ) + guides(color=guide_legend(title="Algorithm")) +
  scale_y_continuous(
    trans='log10'
    #labels = scales::comma,
    #limits=c(1, 1000000),
    #n.breaks = 6,
    #labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
  ) +
  scale_x_continuous(
    trans='log10',
    #labels = scales::comma,
    #limits=c(min((z$query.length + z$target.length)/2), 500000),
    n.breaks = 12
  ) + scale_color_manual(values=c("#39B600", "#FF62BC"))
  #scale_color_gradient(low = "green", high = "red") #
py

zz <- z %>%
  spread(key = mode, value = value) %>%
  mutate(BiWFA_div_WFA = log2(BiWFA/`WFA-high`))
pzz <- ggplot(
  zz,
  aes(
    (query.length + target.length) / 2,
    y = BiWFA_div_WFA,
    alpha=I(0.5)#, color = error.rate
    )
) +
  geom_point(size=1.5) +
  facet_wrap (
    ~set,
    ncol = 2
  ) +
  ggtitle('Average execution time ratio') +
  guides(color=guide_legend(title="Algorithm")) +
  xlab("Length") + ylab(expression("log"[2]("BiWFA / WFA-high"))) +  theme_bw() +
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
    #trans='log10',
    #labels = scales::comma,
    limits=c(-1, 2.5),
    n.breaks = 8,
    #labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
  ) +
  scale_x_continuous(
    #trans='log10'
    #labels = scales::comma,
    limits=c(100, 350000),
    n.breaks = 18
  ) + geom_hline(yintercept=0) #+ scale_color_gradient(low = "green", high = "red") #+ 
  #scale_color_manual(values=c("#39B600", "#FF62BC")) +
  
pzz
ggsave(plot = pzz, paste0('FigureSX', '.pdf'), width = 40, height = 25, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
ggsave(plot = pzz, paste0('FigureSX', '.png'), width = 40, height = 25, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)

if (FALSE) {
  library(dplyr)
  zz <- z %>%
    group_by(set, seq, statistic, query.length, target.length, score) %>%
    summarise_at(vars(value), list(name = diff))
  ggplot(
    zz,
    aes(x = (query.length + target.length) / 2, y = name, alpha=I(1/3))
  ) +
    geom_point() +
    ggtitle('Execution time\nWFA-high - BiWFA') +
    xlab("Length") + ylab("Seconds") +  theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size=18),
      legend.position = "none",
      
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
      #trans='log10'
      #labels = scales::comma,
      #limits=c(1, 1000000),
      #n.breaks = 6,
      #labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
    ) +
    scale_x_continuous(
      #trans='log10',
      #labels = scales::comma,
      limits=c(min((z$query.length + z$target.length)/2), 500000),
      n.breaks = 12
    )
}

# Plot both and save the image
pxy <- ggpubr::ggarrange(
  px, py,
  align='v',
  labels=c('A', 'B'),
  heights=c(1, 1.25),
  font.label = list(size = 18),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
) + 
  theme(plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm")) # to avoid cutting labels
pxy
ggsave(plot = pxy, paste0('FigureS2', '.pdf'), width = 30, height = 20, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
ggsave(plot = pxy, paste0('FigureS2', '.png'), width = 30, height = 20, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
################################################################################


# Statistics
for (current_stat in unique(statsWithMetadata_df$statistic[statsWithMetadata_df$statistic != 'NaN'])) {
  print(current_stat)
  for (current_set in unique(statsWithMetadata_df$set)) {
    print(current_set)
    tmp_df <- statsWithMetadata_df[statsWithMetadata_df$statistic == current_stat & !is.nan(statsWithMetadata_df$value) & statsWithMetadata_df$set == current_set, ]
    tmp_df <- tmp_df[tmp_df$query.length > 10000 & tmp_df$target.length > 10000,]
    
    mean_BiWFA <- mean(tmp_df[tmp_df$mode == 'BiWFA',]$value)
    
    for (current_mode in unique(statsWithMetadata_df$mode)) {
      mean_XXX <- mean(tmp_df[tmp_df$mode == current_mode,]$value)
      
      print(paste0(current_stat, " - ", current_set, " --> ", current_mode, "/BiWFA", ": ", mean_XXX/mean_BiWFA))
      print(paste0(current_stat, " - ", current_set, " --> ", "BiWFA/", current_mode, ": ", mean_BiWFA/mean_XXX))
      print('-------------------')
    }
  }
}

current_stat = 'memory_kb'
# 'ONT PromethION reads vs CHM13 v1.1', 'ONT Ultra Long'
current_set = 'ONT Ultra Long'
current_mode = 'WFA-high'
tmp_df <- statsWithMetadata_df[statsWithMetadata_df$statistic == current_stat & !is.nan(statsWithMetadata_df$value) & statsWithMetadata_df$set == current_set, ]
tmp_df <- tmp_df[tmp_df$query.length > 10000 & tmp_df$target.length > 10000,]
mean_BiWFA <- mean(tmp_df[tmp_df$mode == 'BiWFA',]$value)
mean_XXX <- mean(tmp_df[tmp_df$mode == current_mode,]$value)
mean_XXX/mean_BiWFA


# OOMs
statsWithMetadata_df[statsWithMetadata_df$statistic == 'memory_kb' & is.nan(statsWithMetadata_df$value),]

x <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'memory_kb', ]
x <- x[x$query.length >= 10000 & x$target.length >= 10000,] # To clean outliers
length(unique(x[x$set == 'ONT PromethION reads vs CHM13 v1.1',]$seq))
mean(
  mean(x[x$set == 'ONT PromethION reads vs CHM13 v1.1' & x$mode == 'edlib (edit distance)', ]$query.length),
  mean(x[x$set == 'ONT PromethION reads vs CHM13 v1.1' & x$mode == 'edlib (edit distance)', ]$target.length)
)

length(unique(x[x$set == 'ONT Ultra Long > 500kbps',]$seq))
mean(
  mean(x[x$set == 'ONT Ultra Long > 500kbps' & x$mode == 'edlib (edit distance)', ]$query.length),
  mean(x[x$set == 'ONT Ultra Long > 500kbps' & x$mode == 'edlib (edit distance)', ]$target.length)
)


if (FALSE) {
  x <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'memory_kb' & !is.nan(statsWithMetadata_df$value), ]
  mean_BiWFA_ONTUL_mem <- mean(x[x$set == 'ONT Ultra Long > 500kbps' & x$mode == 'BiWFA',]$value)
  mean_WFAhigh_ONTUL_mem <- mean(x[x$set == 'ONT Ultra Long > 500kbps' & x$mode == 'WFA-high',]$value)
  mean_WFAlow_ONTUL_mem <- mean(x[x$set == 'ONT Ultra Long > 500kbps' & x$mode == 'WFA-low',]$value)
  mean_wfalm_low_ONTUL_mem <- mean(x[x$set == 'ONT Ultra Long > 500kbps' & x$mode == 'wfalm-low',]$value)
  mean_wfalm_rec_ONTUL_mem <- mean(x[x$set == 'ONT Ultra Long > 500kbps' & x$mode == 'wfalm-recursive',]$value)
  mean_edlib_ONTUL_mem <- mean(x[x$set == 'ONT Ultra Long > 500kbps' & x$mode == 'edlib (edit distance)',]$value)
  mean_bitpal_ONTUL_mem <- mean(x[x$set == 'ONT Ultra Long > 500kbps' & x$mode == 'bitpal (score only)',]$value)
  mean_WFAhigh_ONTUL_mem/mean_BiWFA_ONTUL_mem
  mean_WFAlow_ONTUL_mem/mean_BiWFA_ONTUL_mem
  mean_wfalm_low_ONTUL_mem/mean_BiWFA_ONTUL_mem
  mean_wfalm_rec_ONTUL_mem/mean_BiWFA_ONTUL_mem
  mean_BiWFA_ONTUL_mem/mean_edlib_ONTUL_mem
  mean_BiWFA_ONTUL_mem/mean_bitpal_ONTUL_mem
  
  
  mean_BiWFA_ION_mem <- mean(x[x$set == 'ONT PromethION reads vs CHM13 v1.1' & x$mode == 'BiWFA',]$value)
  mean_WFAhigh_ION_mem <- mean(x[x$set == 'ONT PromethION reads vs CHM13 v1.1' & x$mode == 'WFA-high',]$value)
  mean_edlib_ION_mem <- mean(x[x$set == 'ONT PromethION reads vs CHM13 v1.1' & x$mode == 'edlib (edit distance)',]$value)
  mean_bitpal_ION_mem <- mean(x[x$set == 'ONT PromethION reads vs CHM13 v1.1' & x$mode == 'bitpal (score only)',]$value)
  mean_WFAhigh_ION_mem/mean_BiWFA_ION_mem
  mean_BiWFA_ION_mem/mean_edlib_ION_mem
  mean_BiWFA_ION_mem/mean_bitpal_ION_mem
  
  
  y <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'time_s' & !is.nan(statsWithMetadata_df$value), ]
  mean_BiWFA_ONTUL_time <- mean(y[y$set == 'ONT Ultra Long > 500kbps' & y$mode == 'BiWFA',]$value)
  mean_WFAhigh_ONTUL_time <- mean(y[y$set == 'ONT Ultra Long > 500kbps' & y$mode == 'WFA-high',]$value)
  mean_WFAlow_ONTUL_time <- mean(y[y$set == 'ONT Ultra Long > 500kbps' & y$mode == 'WFA-low',]$value)
  mean_wfalm_low_ONTUL_time <- mean(y[y$set == 'ONT Ultra Long > 500kbps' & y$mode == 'wfalm-low',]$value)
  mean_wfalm_rec_ONTUL_time <- mean(y[y$set == 'ONT Ultra Long > 500kbps' & y$mode == 'wfalm-recursive',]$value)
  mean_edlib_ONTUL_time <- mean(y[y$set == 'ONT Ultra Long > 500kbps' & y$mode == 'edlib (edit distance)',]$value)
  mean_bitpal_ONTUL_time <- mean(y[y$set == 'ONT Ultra Long > 500kbps' & y$mode == 'bitpal (score only)',]$value)
  mean_WFAhigh_ONTUL_time/mean_BiWFA_ONTUL_time
  mean_WFAlow_ONTUL_time/mean_BiWFA_ONTUL_time
  mean_wfalm_low_ONTUL_time/mean_BiWFA_ONTUL_time
  mean_wfalm_rec_ONTUL_time/mean_BiWFA_ONTUL_time
  mean_BiWFA_ONTUL_time/mean_edlib_ONTUL_time
  mean_BiWFA_ONTUL_time/mean_bitpal_ONTUL_time
  
  mean_BiWFA_ION_time <- mean(y[y$set == 'ONT PromethION reads vs CHM13 v1.1' & y$mode == 'BiWFA',]$value)
  mean_WFAhigh_ION_time <- mean(y[y$set == 'ONT PromethION reads vs CHM13 v1.1' & y$mode == 'WFA-high',]$value)
  mean_edlib_ION_time <- mean(y[y$set == 'ONT PromethION reads vs CHM13 v1.1' & y$mode == 'edlib (edit distance)',]$value)
  mean_bitpal_ION_time <- mean(y[y$set == 'ONT PromethION reads vs CHM13 v1.1' & y$mode == 'bitpal (score only)',]$value)
  mean_WFAhigh_ION_time/mean_BiWFA_ION_time
  mean_BiWFA_ION_time/mean_edlib_ION_time
  mean_BiWFA_ION_time/mean_bitpal_ION_time
  
  x_wide_df <- x %>% select(-score) %>%pivot_wider(names_from = mode, values_from = value) %>% View()
}

