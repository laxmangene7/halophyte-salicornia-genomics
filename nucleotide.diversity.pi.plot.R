library(ggplot2)
library(data.table)

# Load data
df4 <- fread("nucleotide_diversity_100kb_F_MAF0.01_Miss20_Het20-salicornia.europae.France.genome.1-9C.windowed.pi",
             header = TRUE, check.names = TRUE, data.table = FALSE)

# Plot
p <- ggplot(df4, aes(x = BIN_END, y = PI)) +
  geom_point(size = 1.0, colour = "#023047") +
  geom_smooth(method = "loess", color = "orange", linewidth = 0.8, se = FALSE) +
  scale_x_continuous(
    labels = function(x) round(x / 1e6, 0),  # Mb labels
    expand = c(0.01, 0)
  ) +
  xlab("\nGenomic Position (Mb)") +
  ylab("Nucleotide Diversity (π)\n") + 
  facet_wrap(~ CHROM, nrow = 3, scales = "free_x") +
  theme(
    panel.grid.major.x = element_line(size = 0.3, colour = "gray68"),
    panel.grid.major.y = element_line(size = 0.3, colour = "gray68"),
    panel.background = element_rect(fill = "white"),
    panel.spacing = unit(0.2, "cm"),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 16),
    panel.border = element_rect(colour = "black", size = 0.4, fill = NA),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save plot
ggsave("nucleotide_diversity_100kb_F_MAF0.01_Miss20_Het20-salicornia.europae.France.genome.1-9C.windowed.pi.pdf", p,
       width = 18, height = 13, units = "in", device = cairo_pdf, limitsize = FALSE)




################################################################ merging france and procumbens 
library(ggplot2)
library(data.table)

# Load data
df_perc <- fread(
  "nucleotide_diversity_100kb_tetraploid.F_MAF0.01_Miss20_Het20-tetraploid.S.procumbens.flipped_swapped.windowed.pi",   #CHROM	BIN_START	N_VARIANTS	PI
  data.table = FALSE
)

df_france <- fread(
  "nucleotide_diversity_100kb_tetraploid.F_MAF0.01_Miss20_Het20-tetraploid.S.sp.France.flipped_swapped.windowed.pi",
  data.table = FALSE
)

# Add group labels
df_perc$Group   <- "Salicornia percumbens"
df_france$Group <- "Salicornia spp. France"

# Combine datasets
df_all <- rbind(df_perc, df_france)

# Plot
p <- ggplot(df_all, aes(x = BIN_START, y = PI, color = Group)) +
  geom_point(alpha=0.5, size = 1.0) +
  geom_smooth(
    method = "loess",
    linewidth = 0.8,
    se = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Salicornia percumbens"   = "#E66100",  # orange
      "Salicornia spp. France"  = "#1F4ED8"   # blue
    )
  ) +
  scale_x_continuous(
    labels = function(x) round(x / 1e6, 0),
    expand = c(0.01, 0)
  ) +
  xlab("\nGenomic Position (Mb)") +
  ylab("Nucleotide Diversity (π)\n") +
  facet_wrap(~ CHROM, nrow = 3, scales = "free_x") +
  theme(
    panel.grid.major.x = element_line(size = 0.3, colour = "gray68"),
    panel.grid.major.y = element_line(size = 0.3, colour = "gray68"),
    panel.background = element_rect(fill = "white"),
    panel.spacing = unit(0.2, "cm"),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 16),
    panel.border = element_rect(colour = "black", size = 0.4, fill = NA),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# Save plot
ggsave(
  "nucleotide_diversity_100kb_tetraploid_1A.9A.procumbens.and.S.France.pi.after.chromosome.swapping.and.inversion.pdf",
  p,
  width = 18,
  height = 10,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)




################################################################ merging france and procumbens ## C chromosomes 
library(ggplot2)
library(data.table)

# Load data
df_perc <- fread(
  "nucleotide_diversity_100kb_tetraploid.F_MAF0.01_Miss20_Het20-tetraploid.S.procumbens.flipped_swapped.C.windowed.pi",   #CHROM	BIN_START	N_VARIANTS	PI
  data.table = FALSE
)

df_france <- fread(
  "nucleotide_diversity_100kb_tetraploid.F_MAF0.01_Miss20_Het20-tetraploid.S.France.flipped_swapped.C.windowed.pi",
  data.table = FALSE
)

# Add group labels
df_perc$Group   <- "Salicornia percumbens"
df_france$Group <- "Salicornia spp. France"

# Combine datasets
df_all <- rbind(df_perc, df_france)

# Plot
p <- ggplot(df_all, aes(x = BIN_START, y = PI, color = Group)) +
  geom_point(alpha=0.5, size = 1.0) +
  geom_smooth(
    method = "loess",
    linewidth = 0.8,
    se = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Salicornia percumbens"   = "#E66100",  # orange
      "Salicornia spp. France"  = "#1F4ED8"   # blue
    )
  ) +
  scale_x_continuous(
    labels = function(x) round(x / 1e6, 0),
    expand = c(0.01, 0)
  ) +
  xlab("\nGenomic Position (Mb)") +
  ylab("Nucleotide Diversity (π)\n") +
  facet_wrap(~ CHROM, nrow = 3, scales = "free_x") +
  theme(
    panel.grid.major.x = element_line(size = 0.3, colour = "gray68"),
    panel.grid.major.y = element_line(size = 0.3, colour = "gray68"),
    panel.background = element_rect(fill = "white"),
    panel.spacing = unit(0.2, "cm"),
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold", size = 16),
    panel.border = element_rect(colour = "black", size = 0.4, fill = NA),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# Save plot
ggsave(
  "nucleotide_diversity_100kb_tetraploid_1C.9C.procumbens.and.S.France.pi.after.chromosome.swapping.and.inversion.pdf",
  p,
  width = 18,
  height = 10,
  units = "in",
  device = cairo_pdf,
  limitsize = FALSE
)




