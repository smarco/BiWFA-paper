library(ggplot2)
library(ggforce)
options(scipen = 9)

stats_df <- read.table('statistics_all.tsv', sep = '\t', header = F)
names(stats_df) <- c('set', 'seq', 'mode', 'value', 'statistic')
lengths_df <- read.table('lengths_all.tsv', sep = '\t', header = F)
names(lengths_df) <- c('set', 'seq', 'query.length', 'target.length')
scores_df <- read.table('scores_all.tsv', sep = '\t', header = T)
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

# Rename datasets
levels(statsWithMetadata_df$set)[match("ont_regions",levels(statsWithMetadata_df$set))] <- "ONT PromethION reads vs CHM13 v1.1"
levels(statsWithMetadata_df$set)[match("ONT_UL",levels(statsWithMetadata_df$set))] <- "ONT Ultra Long > 500kbps"

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


# Plot memory use
x <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'memory_kb' & !is.nan(statsWithMetadata_df$value), ]
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
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1)
  ) +
  guides(fill=guide_legend(title="Algorithm")) +
  ggtitle('Memory consumption') +
  xlab("Algorithm") + ylab("")

# Plot runtime
y <- statsWithMetadata_df[statsWithMetadata_df$statistic == 'time_s' & !is.nan(statsWithMetadata_df$value), ]
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
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top",
    axis.text.x = element_text(angle = 30, vjust = 1.0, hjust=1)
  ) +
  guides(fill=guide_legend(title="Algorithm")) +
  ggtitle('Execution time') +
  xlab("Algorithm") + ylab("Seconds")
py

# Plot both and save the image
library(ggpubr)
pxy <- ggpubr::ggarrange(px, py, align='hv', labels=c('A', 'B'),legend = "right", # legend position,
                         common.legend = T, nrow = 2)
ggsave(plot = pxy, paste0('Figure1', '.pdf'), width = 30, height = 20, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
