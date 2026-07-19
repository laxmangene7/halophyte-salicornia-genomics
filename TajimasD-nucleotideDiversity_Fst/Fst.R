#### FST 
############################################################

fst <- read.table(
  "Fst_10kb.window.5kb.step.S_procumbens.vs.S.sp.France.weir.pop.fst",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

fst <- fst %>%
  filter(
    CHROM == "chr6",
    !is.na(MEAN_FST)
  )

fst$BIN_START <- as.numeric(fst$BIN_START)
fst$BIN_END   <- as.numeric(fst$BIN_END)
fst$MEAN_FST  <- as.numeric(fst$MEAN_FST)

fst$MID <- (
  fst$BIN_START +
    fst$BIN_END
) / 2

fst$MEAN_FST[
  fst$MEAN_FST < 0
] <- 0


############################################################
# OSCA position (bp -> Mb)
############################################################
osca_mb <- 35809855 / 1e6

############################################################
# FST plot
############################################################

p3 <- ggplot(
  fst,
  aes(x = MID / 1e6, y = MEAN_FST)
) +
  
  geom_point(
    size = 1,
    alpha = 0.6,
    color = "#444444"
  ) +
  
  geom_smooth(
    method = "loess",
    span = 0.15,
    se = FALSE,
    color = "red",
    linewidth = 1.2
  ) +
  
  geom_vline(
    xintercept = osca_mb,
    color = "blue",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  labs(
    title = "FST",
    x = "Position on chr6 (Mb)",
    y = "Mean FST"
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
    axis.text = element_text(size = 14)
  )

print(p3)

# Save plot
ggsave(
  filename = "Supplementary.Fig20_FST_plot_chr6.10kb.window.manuscript.pdf",
  plot = p3,
  width = 14,
  height = 5,
  dpi = 300
)


## Optionally save as svg

ggsave(
  filename = "Supplementary.Fig20_FST_plot_chr6.10kb.window.manuscript.svg",
  plot = p3,
  width = 14,
  height = 5,
  dpi = 300
)






