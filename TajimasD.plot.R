

#### merged plot ##### chr6A #################################
library(ggplot2)
library(dplyr)

# ---------- Read dataset 1 ----------
tajima_euro <- read.table(
  "Chr6A_F_MAF0.01_Miss20_Het20_Fr.ref.FrancA.flipped.D",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(tajima_euro) <- c("CHROM", "BIN_START", "BIN_END", "TajimaD")

tajima_euro <- tajima_euro %>%
  mutate(
    POS_KB = BIN_START / 1000,
    Dataset = "Tetraploid europaea"
  )

# ---------- Read dataset 2 ----------
tajima_perc <- read.table(
  "Chr6A_F_MAF0.01_Miss20_Het20-salicornia.percumbens_Flipped.clean.D",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(tajima_perc) <- c("CHROM", "BIN_START", "BIN_END", "TajimaD")

tajima_perc <- tajima_perc %>%
  mutate(
    POS_KB = BIN_START / 1000,
    Dataset = "Salicornia perumbens"
  )

# ---------- Merge datasets ----------
tajima_all <- bind_rows(tajima_euro, tajima_perc)

# ---------- Plot ----------
ggplot(tajima_all, aes(x = POS_KB, y = TajimaD, color = Dataset)) +
  
  geom_point(alpha = 0.8, size = 1.4) +
  
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = FALSE,
    linewidth = 0.9
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  scale_color_manual(
    values = c(
      "Tetraploid europaea" = "#1F4ED8",   # blue
      "Salicornia perumbens" = "#E66100"   # orange
    )
  ) +
  
  labs(
    title = "Comparison of Tajima's D across Chr6A",
    x = "Genomic Position (kb)",
    y = "Tajima's D",
    color = "Dataset",
    caption = "LOESS-smoothed Tajima’s D profiles along Chr6A"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# ---------- Save ----------
ggsave(
  "Chr6A_TajimaD_merged_france.tetraploid_vs_percumbens.flipped.orientation.svg",
  width = 6,
  height = 4.5,
  dpi = 300
)




##################### merged all 9 chromosomes ################################################



library(ggplot2)
library(dplyr)

# -----------------------------
# Read tetraploid europaea data
# -----------------------------
tajima_euro <- read.table(
  "F_MAF0.01_Miss20_Het20-tetraploid.euro.Fr.ref.FranA.vcf.recode.Tajima.D",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Tetraploid europaea"
  )

# -----------------------------
# Read Salicornia perumbens data
# -----------------------------
tajima_perc <- read.table(
  "F_MAF0.01_Miss20_Het20-salicornia.percumbens.ref.FranA.tajima.Tajima.D",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia perumbens"
  )

# -----------------------------
# Merge datasets
# -----------------------------
tajima_all <- bind_rows(tajima_euro, tajima_perc) %>%
  
  # keep only Chr1A–Chr9A and order properly
  filter(CHROM %in% paste0("Chr", 1:9, "A")) %>%
  mutate(
    CHROM = factor(CHROM, levels = paste0("Chr", 1:9, "A"))
  )

# -----------------------------
# Plot merged Tajima's D
# -----------------------------
ggplot(tajima_all, aes(x = POS_MB, y = TajimaD, color = Dataset)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = FALSE,
    linewidth = 0.8
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  facet_wrap(~ CHROM, ncol = 3, scales = "free_x") +
  
  scale_color_manual(
    values = c(
      "Tetraploid europaea" = "#1F4ED8",   # blue
      "Salicornia perumbens" = "#E66100"   # orange
    )
  ) +
  
  labs(
    title = "Tajima’s D across A-subgenome chromosomes (Chr1A–Chr9A)",
    x = "Genomic position (Mb)",
    y = "Tajima’s D",
    color = "Dataset",
    caption = "LOESS-smoothed Tajima’s D profiles; dashed line indicates neutral expectation (D = 0)"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  )

# -----------------------------
# Save figure
# -----------------------------
ggsave(
  "TajimaD_Chr1A_Chr9A_merged_euro_vs_percumbens.png",
  width = 13,
  height = 8,
  dpi = 600
)




## increase font size ##  chr1A-chr9A
# ============================================
# Tajima's D comparison:
# Salicornia sp (France) vs S. percumbens
# A-subgenome (Chr1A–Chr9A)
# ============================================

library(ggplot2)
library(dplyr)
library(grid)   # for unit()

# -----------------------------
# 1. Read Salicornia sp France
# -----------------------------
tajima_euro <- read.table(
  "F_MAF0.01_Miss20_Het20-tetraploid.euro.Fr.ref.FranA.vcf.recode.Tajima.flipped.clean1.D",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia sp France"
  )

# -----------------------------
# 2. Read Salicornia percumbens
# -----------------------------
tajima_perc <- read.table(
  "F_MAF0.01_Miss20_Het20-salicornia.percumbens_flipped.Asub.D",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia percumbens"
  )

# -----------------------------
# 3. Merge and filter chromosomes
# -----------------------------
tajima_all <- bind_rows(tajima_euro, tajima_perc) %>%
  filter(CHROM %in% paste0("Chr", 1:9, "A")) %>%
  mutate(
    CHROM = factor(CHROM, levels = paste0("Chr", 1:9, "A"))
  )

# -----------------------------
# 4. Plot
# -----------------------------
p <- ggplot(tajima_all, aes(x = POS_MB, y = TajimaD, color = Dataset)) +
  
  geom_point(alpha = 0.6, size = 1.3) +
  
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = FALSE,
    linewidth = 1.2
  ) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.8
  ) +
  
  facet_wrap(~ CHROM, ncol = 3, scales = "free_x") +
  
  scale_color_manual(
    values = c(
      "Salicornia sp France" = "#1F4ED8",
      "Salicornia percumbens" = "#E66100"
    )
  ) +
  
  labs(
    title = "Tajima’s D across A-subgenome chromosomes (Chr1A–Chr9A)",
    x = "Genomic position (Mb)",
    y = "Tajima’s D",
    color = "Dataset",
    caption = "LOESS-smoothed Tajima’s D profiles; dashed line indicates neutral expectation (D = 0)"
  ) +
  
  theme_minimal(base_size = 22) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 28),
    axis.title = element_text(size = 24, face = "bold"),
    axis.text = element_text(size = 20),
    strip.text = element_text(face = "bold", size = 22),
    legend.title = element_text(size = 22, face = "bold"),
    legend.text = element_text(size = 20),
    legend.position = "top",
    plot.caption = element_text(size = 18, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )

# -----------------------------
# 5. Save figure
# -----------------------------
ggsave(
  filename = "TajimaD_Chr1A_Chr9A_merged_France_vs_percumbens.flipped.orientation.pdf",
  plot = p,
  width = 22,
  height = 18,
  dpi = 600
)








##### Merged plot tajimas D procumbens and france 1C-9C chromosomes 

library(ggplot2)
library(dplyr)

# -----------------------------
# Read tetraploid europaea data
# -----------------------------
tajima_euro <- read.table(
  "F_MAF0.01_Miss20_Het20-salicornia.europae.France.genome.1-9C.merged.Tajima.D",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Tetraploid europaea"
  )

# -----------------------------
# Read Salicornia perumbens data
# -----------------------------
tajima_perc <- read.table(
  "F_MAF0.01_Miss20_Het20-salicornia.procumbens.genome.1-9C.procumbens.merged.Tajima.D",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia perumbens"
  )

# -----------------------------
# Merge datasets
# -----------------------------
tajima_all <- bind_rows(tajima_euro, tajima_perc) %>%
  
  # keep only Chr1A–Chr9A and order properly
  filter(CHROM %in% paste0("Chr", 1:9, "C")) %>%
  mutate(
    CHROM = factor(CHROM, levels = paste0("Chr", 1:9, "C"))
  )

# -----------------------------
# Plot merged Tajima's D
# -----------------------------
ggplot(tajima_all, aes(x = POS_MB, y = TajimaD, color = Dataset)) +
  
  geom_point(alpha = 0.7, size = 1.2) +
  
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = FALSE,
    linewidth = 0.8
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  facet_wrap(~ CHROM, ncol = 3, scales = "free_x") +
  
  scale_color_manual(
    values = c(
      "Tetraploid europaea" = "#1F4ED8",   # blue
      "Salicornia perumbens" = "#E66100"   # orange
    )
  ) +
  
  labs(
    title = "Tajima’s D across C-subgenome chromosomes (Chr1C–Chr9C)",
    x = "Genomic position (Mb)",
    y = "Tajima’s D",
    color = "Dataset",
    caption = "LOESS-smoothed Tajima’s D profiles; dashed line indicates neutral expectation (D = 0)"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  )

# -----------------------------
# Save figure
# -----------------------------
ggsave(
  "TajimaD_Chr1C_Chr9C_merged_euro_vs_percumbens.png",
  width = 13,
  height = 8,
  dpi = 600
)









## increase font size ## 
library(ggplot2)
library(dplyr)

# -----------------------------
# Read tetraploid europaea data
# -----------------------------
tajima_euro <- read.table(
  "F_MAF0.01_Miss20_Het20-salicornia.europae.France.genome.1-9C.merged.Tajima.D",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia sp France"
  )

# -----------------------------
# Read Salicornia perumbens data
# -----------------------------
tajima_perc <- read.table(
  "F_MAF0.01_Miss20_Het20-salicornia.procumbens.genome.1-9C.procumbens.merged.Tajima.D",
  header = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia percumbens"
  )

# -----------------------------
# Merge datasets
# -----------------------------
tajima_all <- bind_rows(tajima_euro, tajima_perc) %>%
  filter(CHROM %in% paste0("Chr", 1:9, "C")) %>%
  mutate(
    CHROM = factor(CHROM, levels = paste0("Chr", 1:9, "C"))
  )

# -----------------------------
# Plot merged Tajima's D
# -----------------------------
ggplot(tajima_all, aes(x = POS_MB, y = TajimaD, color = Dataset)) +
  
  geom_point(alpha = 0.7, size = 1.5) +
  
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = FALSE,
    linewidth = 1.2
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  
  facet_wrap(~ CHROM, ncol = 3, scales = "free_x") +
  
  scale_color_manual(
    values = c(
      "Salicornia sp France" = "#1F4ED8",
      "Salicornia percumbens" = "#E66100"
    )
  ) +
  
  labs(
    title = "Tajima’s D across C-subgenome chromosomes (Chr1C–Chr9C)",
    x = "Genomic position (Mb)",
    y = "Tajima’s D",
    color = "Dataset",
    caption = "LOESS-smoothed Tajima’s D profiles; dashed line indicates neutral expectation (D = 0)"
  ) +
  
  theme_minimal(base_size = 24) +
  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
    axis.title = element_text(size = 26, face = "bold"),
    axis.text = element_text(size = 22),
    strip.text = element_text(face = "bold", size = 24),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 22),
    legend.position = "top",
    plot.caption = element_text(size = 20, hjust = 0),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )

# -----------------------------
# Save figure
# -----------------------------
ggsave(
  "TajimaD_Chr1C_Chr9C_merged_euro_vs_percumbens.pdf",
  width = 23,
  height = 18,
  dpi = 600
)




