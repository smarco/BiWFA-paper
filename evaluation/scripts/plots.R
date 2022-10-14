library(ggplot2)
library(ggforce) # On Ubuntu 22.04 LTS: sudo apt -y install libfontconfig1-dev
library(ggpubr)

options(scipen = 9)

stats_df <- read.table('~/git/BiWFA-paper/evaluation/data/statistics_all.tsv', sep = '\t', header = F)
names(stats_df) <- c('set', 'seq', 'mode', 'value', 'statistic')
lengths_df <- read.table('~/git/BiWFA-paper/evaluation/data/lengths_all.tsv', sep = '\t', header = F)
names(lengths_df) <- c('set', 'seq', 'query.length', 'target.length')
scores_df <- read.table('~/git/BiWFA-paper/evaluation/data/scores_all.tsv', sep = '\t', header = T)
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

# Force column's type
statsWithMetadata_df$set <- as.factor(statsWithMetadata_df$set)
statsWithMetadata_df$mode <- as.factor(statsWithMetadata_df$mode)

# Rename datasets
levels(statsWithMetadata_df$set)[match("ont_regions",levels(statsWithMetadata_df$set))] <- "ONT PromethION reads vs CHM13 v1.1"
levels(statsWithMetadata_df$set)[match("ONT_UL",levels(statsWithMetadata_df$set))] <- "ONT Ultra Long > 500kbps"
levels(statsWithMetadata_df$set)[match("ONT_UL_SHORT",levels(statsWithMetadata_df$set))] <- "ONT Ultra Long <= 10kbps"

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
# Sequences >= 10kbps

# Plot memory use
x <- statsWithMetadata_df[statsWithMetadata_df$set %in% c(
  "ONT PromethION reads vs CHM13 v1.1", "ONT Ultra Long > 500kbps"
),]
x <- x[x$statistic == 'memory_kb' & !is.nan(x$value), ]
x <- x[x$query.length >= 10000 & x$target.length >= 10000,] # To clean outliers
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
  ggtitle('Memory consumption') +
  xlab("Algorithm") + ylab("") + geom_vline(xintercept=2.5)
px 

# Plot runtime
y <- statsWithMetadata_df[statsWithMetadata_df$set %in% c(
  "ONT PromethION reads vs CHM13 v1.1", "ONT Ultra Long > 500kbps"
),]
y <- y[y$statistic == 'time_s' & !is.nan(y$value), ]
y <- y[y$query.length >= 10000 & y$target.length >= 10000,] # To clean outliers
py <- ggplot(y, aes(x = mode, y = value + 0.001, fill=mode)) +
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
################################################################################

################################################################################
# Sequences <= 10kbps

# Plot memory use
x <- statsWithMetadata_df[statsWithMetadata_df$set %in% c(
  "ONT PromethION reads vs CHM13 v1.1", "ONT Ultra Long <= 10kbps"
),]
x <- x[x$statistic == 'memory_kb' & !is.nan(x$value), ]
x <- x[x$query.length <= 10000 & x$target.length <= 10000,]
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
  ggtitle('Memory consumption') +
  xlab("Algorithm") + ylab("") + geom_vline(xintercept=2.5)
px

# Plot runtime
y <- statsWithMetadata_df[statsWithMetadata_df$set %in% c(
  "ONT PromethION reads vs CHM13 v1.1", "ONT Ultra Long > 500kbps", "ONT Ultra Long <= 10kbps"
),]
y <- y[y$statistic == 'time_s' & !is.nan(y$value), ]
y <- y[y$query.length <= 10000 & y$target.length <= 10000,]
py <- ggplot(y, aes(x = mode, y = value + 0.001, fill=mode)) +
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
ggsave(plot = pxy, paste0('FigureS1', '.png'), width = 30, height = 20, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
################################################################################


################################################################################
# WFA-high vs BiWFA

# Time
z <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'time_s' & !is.nan(statsWithMetadata_df$value), ]
levels(z$set)[match("ONT_UL_OTHER",levels(z$set))] <- "ONT Ultra Long"
#levels(z$set)[match("ONT Ultra Long <= 10kbps",levels(z$set))] <- "ONT Ultra Long"
#levels(z$set)[match("ONT Ultra Long > 500kbps",levels(z$set))] <- "ONT Ultra Long"
z <- z[z$mode %in% c('WFA-high', 'BiWFA') & z$set %in% c('ONT Ultra Long'),]
px <- ggplot(
  z,
  aes(x = (query.length + target.length) / 2, y = value + 0.001, color = mode, alpha=I(1/3))
) +
  geom_point() +
  #facet_wrap (
  #  ~set,
  #  ncol = 2
  #) +
  ggtitle('Execution time') +
  xlab("Length (base pairs)") + ylab("Seconds") +  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "top",
    
    axis.title=element_text(size=18),
    
    #axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1, size=15),
    axis.title.x = element_blank(),
    axis.text.x=element_blank(), #remove x axis labels
    #axis.ticks.x=element_blank(), #remove x axis ticks
    
    axis.text.y=element_text(size=15),
    #axis.text.y=element_blank(),  #remove y axis labels
    #axis.ticks.y=element_blank()  #remove y axis ticks
    
    strip.text = element_text(size=15),
    
    legend.title = element_text(size=15),
    legend.text = element_text(size=14)
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
    limits=c(min((z$query.length + z$target.length)/2), 500000),
    n.breaks = 12
  ) + 
  scale_color_manual(values=c("#39B600", "#FF62BC"))
px 
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

# Memory
z <- statsWithMetadata_df[statsWithMetadata_df$statistic != 'time_s' & !is.nan(statsWithMetadata_df$value), ]
levels(z$set)[match("ONT_UL_OTHER",levels(z$set))] <- "ONT Ultra Long"
#levels(z$set)[match("ONT Ultra Long <= 10kbps",levels(z$set))] <- "ONT Ultra Long"
#levels(z$set)[match("ONT Ultra Long > 500kbps",levels(z$set))] <- "ONT Ultra Long"
z <- z[z$mode %in% c('WFA-high', 'BiWFA') & z$set %in% c('ONT Ultra Long'),]
py <- ggplot(
  z,
  aes(x = (query.length + target.length) / 2, y = value / 1024, color = mode, alpha=I(1/3))
) +
  geom_point() +
  #facet_wrap (
  #  ~set,
  #  ncol = 2
  #) +
  ggtitle('Memory consumption') +
  xlab("Length (base pairs)") + ylab("") +  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.position = "top",
    
    axis.title=element_text(size=18),
    
    axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1, size=15),
    #axis.title.x = element_blank(),
    #axis.text.x=element_blank(), #remove x axis labels
    #axis.ticks.x=element_blank(), #remove x axis ticks
    
    axis.text.y=element_text(size=15),
    #axis.text.y=element_blank(),  #remove y axis labels
    #axis.ticks.y=element_blank()  #remove y axis ticks
    
    strip.text = element_text(size=15),
    
    legend.title = element_text(size=15),
    legend.text = element_text(size=14),
  ) + guides(color=guide_legend(title="Algorithm")) +
  scale_y_continuous(
    trans='log10',
    #labels = scales::comma,
    limits=c(1, 1000000),
    n.breaks = 6,
    labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
  ) +
  scale_x_continuous(
    trans='log10',
    #labels = scales::comma,
    limits=c(min((z$query.length + z$target.length)/2), 500000),
    n.breaks = 12
  ) + 
  scale_color_manual(values=c("#39B600", "#FF62BC"))
py

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
################################################################################
