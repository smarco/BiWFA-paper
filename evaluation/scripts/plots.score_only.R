library(ggplot2)
library(ggforce) # On Ubuntu 22.04 LTS: sudo apt -y install libfontconfig1-dev
library(ggpubr)

options(scipen = 9)

stats_df <- read.table('~/git/BiWFA-paper/evaluation/data/statistics.alignment_vs_score_only.tsv', sep = '\t', header = F)
names(stats_df) <- c('set', 'seq', 'mode', 'value', 'statistic')
lengths_df <- read.table('~/git/BiWFA-paper/evaluation/data/statistics.alignment_vs_score_only.tsv', sep = '\t', header = F)
names(lengths_df) <- c('set', 'seq', 'query.length', 'target.length')

statsWithMetadata_df <- merge(
  stats_df,
  lengths_df,
  all.x=T,
  by.x = c('set', 'seq'), by.y = c('set', 'seq')
)

# Force column's type
statsWithMetadata_df$set <- as.factor(statsWithMetadata_df$set)
statsWithMetadata_df$mode <- as.factor(statsWithMetadata_df$mode)

# Rename datasets
levels(statsWithMetadata_df$set)[match("ONT_UL_OTHER",levels(statsWithMetadata_df$set))] <- "Full traceback"
levels(statsWithMetadata_df$set)[match("ONT_UL_OTHER_SCORE_ONLY",levels(statsWithMetadata_df$set))] <- "Score only"

# Rename algorithms and change their order
levels(statsWithMetadata_df$mode)[match("edlib",levels(statsWithMetadata_df$mode))] <- "edlib (edit distance)"
levels(statsWithMetadata_df$mode)[match("bitpal-scored",levels(statsWithMetadata_df$mode))] <- "bitpal (score only)"
levels(statsWithMetadata_df$mode)[match("wfa-ultralow",levels(statsWithMetadata_df$mode))] <- "BiWFA"
levels(statsWithMetadata_df$mode)[match("wfa-high",levels(statsWithMetadata_df$mode))] <- "WFA-high"
#levels(statsWithMetadata_df$mode)[match("wfa-low",levels(statsWithMetadata_df$mode))] <- "WFA-low"
#levels(statsWithMetadata_df$mode)[match("wfa-med",levels(statsWithMetadata_df$mode))] <- "WFA-med"
#levels(statsWithMetadata_df$mode)[match("wfalm-lowmem",levels(statsWithMetadata_df$mode))] <- "wfalm-low"
#levels(statsWithMetadata_df$mode)[match("wfalm-rec",levels(statsWithMetadata_df$mode))] <- "wfalm-recursive"

statsWithMetadata_df$mode <- factor(
  statsWithMetadata_df$mode,
  levels = c("edlib (edit distance)", "bitpal (score only)", "ksw2-extz2-sse", "WFA-high", "WFA-med", "WFA-low",  "wfalm", "wfalm-low", "wfalm-recursive", "BiWFA")
)

# Plot memory use
x <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'memory_kb' & !is.nan(statsWithMetadata_df$value), ]
px <- ggplot(x, aes(x = set, y = value / 1024, fill=mode)) +
  geom_boxplot() +
  scale_y_continuous(
    trans='log10',
    #labels = scales::comma,
    limits=c(1, 1000000),
    n.breaks = 6,
    labels=c("NA" = "", "1" = "1 MB", "10" = "10 MB", "100" = "100 MB", "1000" = "1 GB", "10000" = "10 GB", "100000" = "100 GB", "1000000" = "1 TB", "NA" = "")
  ) +
  facet_wrap (
    ~mode,
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
  xlab("Algorithm") + ylab("") + 
  scale_fill_manual(values=c("#39B600", "#FF62BC"))
px 

# Plot runtime
y <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'time_s' & !is.nan(statsWithMetadata_df$value), ]
py <- ggplot(y, aes(x = set, y = value + 0.001, fill=mode)) +
  geom_boxplot() +
  scale_y_continuous(
    trans='log10',
    labels = scales::comma
  ) +
  facet_wrap (
    ~mode,
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
  xlab("Algorithm") + ylab("Seconds") + 
  scale_fill_manual(values=c("#39B600", "#FF62BC"))
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
ggsave(plot = pxy, paste0('FigureS4', '.pdf'), width = 15, height = 20, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
ggsave(plot = pxy, paste0('FigureS4', '.png'), width = 15, height = 20, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)

library(dplyr)
statsWithMetadata_df[statsWithMetadata_df$statistic == 'time_s' & !is.nan(statsWithMetadata_df$value), ] %>%
  dplyr::select(set, mode, value) %>%
  dplyr::filter(mode == 'BiWFA') %>%
  dplyr::group_by(set) %>%
  dplyr::summarise_at(vars(value), list(name = mean))
