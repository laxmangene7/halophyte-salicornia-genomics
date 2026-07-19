############################################################
# Chr6A population genomics comparison
# Salicornia sp. France (12 random accessions) vs S. procumbens
############################################################

library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)

############################################################
# COLORS
############################################################

plot_colors <- c(
  "Salicornia sp France" = "#1F4ED8",
  "Salicornia procumbens" = "#E66100"
)

############################################################
# OSCA gene start (optional vertical line)
############################################################

OSCA_START_MB <- 35807352 / 1e6

############################################################
# 1. TAJIMA'S D
############################################################

tajima_france <- read.table(
  "S_sp.France.subset1.random.sample.called.on.S.ramosissima..50kb.window.Tajima.D.txt",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia sp France"
  )

tajima_proc <- read.table(
  "S_procumbens.sites.with.S_ramosissima.ref.50kb.window_Chr6A_tajima_D.txt",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia procumbens"
  )

tajima_all <- bind_rows(
  tajima_france,
  tajima_proc
) %>%
  filter(CHROM == "chr6")

############################################################
# Tajima's D plot
############################################################

p1 <- ggplot(
  tajima_all,
  aes(
    x = POS_MB,
    y = TajimaD,
    color = Dataset
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.5
  ) +
  geom_smooth(
    method = "loess",
    span = 0.25,
    se = FALSE,
    linewidth = 1.2
  ) +
  
  ## Optional: OSCA gene position
  geom_vline(
    xintercept = OSCA_START_MB,
    linetype = "dashed",
    linewidth = 0.8,
    color = "steelblue"
  ) +
  
  scale_color_manual(values = plot_colors) +
  
  labs(
    title = "Tajima's D",
    x = "Chromosome 6A position (Mb)",
    y = "Tajima's D",
    color = NULL
  ) +
  
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold",
      size = 18
    ),
    axis.text = element_text(size = 14),
    legend.position = "top",
    legend.text = element_text(size = 14)
  )



############################################################
# Save Tajima's D plot
############################################################

ggsave(
  filename = "Chr6A_TajimasD_subset_vs_procumbens.pdf",
  plot = p1,
  width = 10,
  height = 5,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "Chr6A_TajimasD_subset_vs_procumbens.svg",
  plot = p1,
  width = 10,
  height = 5,
  units = "in"
)

