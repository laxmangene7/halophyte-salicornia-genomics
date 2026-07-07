################## osca region  ####################################


############################################################
# Chr6 population genomics comparison
# Panels:
# 1. Tajima's D
# 2. Nucleotide diversity (π)
############################################################

library(ggplot2)
library(dplyr)
library(data.table)
library(patchwork)

############################################################
# OSCA candidate region
############################################################

osca_start_bp <- 35807352
osca_end_bp   <- 35815194

osca_start_mb <- osca_start_bp / 1e6
osca_end_mb   <- osca_end_bp / 1e6

############################################################
# COLORS
############################################################

plot_colors <- c(
  "Salicornia sp France" = "#1F4ED8",
  "Salicornia procumbens" = "#E66100"
)

############################################################
# 1. TAJIMA'S D
############################################################

tajima_france <- read.table(
  "S_sp.France.sites.with.S_ramosissima.ref.50kb.window_Chr6A_tajimaD.txt",
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
  
  geom_rect(
    aes(
      xmin = osca_start_mb,
      xmax = osca_end_mb,
      ymin = -Inf,
      ymax = Inf
    ),
    inherit.aes = FALSE,
    fill = "grey40",
    alpha = 0.2
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
  
  scale_color_manual(
    values = plot_colors
  ) +
  
  labs(
    title = "Tajima's D",
    y = "Tajima's D"
  ) +
  
  theme_bw(base_size = 18) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = 0.5
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      face = "bold",
      size = 18
    ),
    
    axis.text = element_text(size = 14),
    
    legend.position = "top",
    
    legend.title = element_text(
      face = "bold",
      size = 16
    ),
    
    legend.text = element_text(size = 14)
  )




############################################################
# 2. NUCLEOTIDE DIVERSITY (PI)
############################################################

pi_france <- fread(
  "S_sp.France.sites.with.S_ramosissima.ref.50kb.window_Chr6A_pi",
  header = TRUE,
  check.names = TRUE,
  data.table = FALSE
) %>%
  mutate(
    POS_MB = BIN_END / 1e6,
    Dataset = "Salicornia sp France"
  )

pi_proc <- fread(
  "S_procumbens.sites.with.S_ramosissima.ref.50kb.window_Chr6A_pi",
  header = TRUE,
  check.names = TRUE,
  data.table = FALSE
) %>%
  mutate(
    POS_MB = BIN_END / 1e6,
    Dataset = "Salicornia procumbens"
  )

pi_all <- bind_rows(
  pi_france,
  pi_proc
) %>%
  filter(
    CHROM == "chr6",
    !is.na(PI)
  )

############################################################
# PI plot
############################################################

p2 <- ggplot(
  pi_all,
  aes(
    x = POS_MB,
    y = PI,
    color = Dataset
  )
) +
  
  geom_rect(
    aes(
      xmin = osca_start_mb,
      xmax = osca_end_mb,
      ymin = -Inf,
      ymax = Inf
    ),
    inherit.aes = FALSE,
    fill = "grey40",
    alpha = 0.2
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
  
  scale_color_manual(
    values = plot_colors
  ) +
  
  labs(
    title = expression(
      Nucleotide~diversity~(pi)
    ),
    y = expression(pi)
  ) +
  
  theme_bw(base_size = 18) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 20,
      hjust = 0.5
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      face = "bold",
      size = 18
    ),
    
    axis.text = element_text(size = 14),
    
    legend.position = "none"
  )





############################################################
# Combine panels
############################################################

combined_plot <-
  p1 /
  p2 +
  plot_layout(
    heights = c(1,1)
  )

############################################################
# Save PDF
############################################################

ggsave(
  "Chr6_population_genomics_test.tajimas.D.pi.pdf",
  combined_plot,
  width = 16,
  height = 16,
  dpi = 600
)

############################################################
# Save SVG
############################################################

ggsave(
  "Chr6_population_genomics_panels.ramosissima.svg",
  combined_plot,
  width = 16,
  height = 16
)

############################################################
# Print
############################################################

print(combined_plot)
