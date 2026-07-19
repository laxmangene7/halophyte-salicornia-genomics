############################################################
# R_script.plotting.R
#
# Plot normalized genome coverage across chromosomes
#
# Input:
#   *.norm.txt
#
# Output:
#   *_Normalized.reads.pdf
#
# Each input file should contain four columns:
#   1. Raw read count
#   2. Chromosome
#   3. 100-kb genomic bin position
#   4. Normalized read count
############################################################

library(data.table)
library(ggplot2)

############################################################
# Plot every normalized coverage file in the directory
############################################################

files <- list.files(
  pattern = "\\.norm\\.txt$",
  full.names = TRUE
)

for (file in files) {

  ##########################################################
  # Read normalized coverage file
  ##########################################################

  df <- fread(
    file,
    header = FALSE,
    data.table = FALSE
  )

  colnames(df) <- c(
    "Raw_Reads",
    "Chromosome",
    "Position",
    "Normalized_Reads"
  )

  ##########################################################
  # Sample name
  ##########################################################

  sample_name <- sub(
    "\\.norm\\.txt$",
    "",
    basename(file)
  )

  ##########################################################
  # Plot
  ##########################################################

  p <- ggplot(
    subset(df, Normalized_Reads < 3),
    aes(
      x = Position,
      y = Normalized_Reads
    )
  ) +

    geom_point(
      size = 1,
      colour = "blue"
    ) +

    facet_wrap(
      ~Chromosome,
      ncol = 1,
      strip.position = "right"
    ) +

    scale_x_continuous(
      labels = function(x) x / 1e6,
      breaks = seq(0, 100000000, by = 10000000),
      expand = c(0.01, 0)
    ) +

    labs(
      title = paste("Sample:", sample_name),
      x = "Genomic Position (Mb)",
      y = "Normalized Read Counts"
    ) +

    theme_bw(base_size = 13) +

    theme(

      panel.grid.major.x =
        element_line(
          colour = "grey75",
          linewidth = 0.2
        ),

      panel.grid.major.y =
        element_line(
          colour = "grey75",
          linewidth = 0.2
        ),

      panel.background =
        element_rect(fill = "white"),

      panel.border =
        element_rect(
          colour = "black",
          fill = NA,
          linewidth = 0.4
        ),

      panel.spacing =
        unit(1.2, "lines"),

      strip.text.y.right =
        element_text(
          angle = 0,
          face = "bold",
          size = 12
        ),

      axis.text =
        element_text(
          size = 10,
          face = "bold"
        ),

      axis.title =
        element_text(
          size = 13,
          face = "bold"
        ),

      plot.title =
        element_text(
          size = 15,
          face = "bold",
          hjust = 0.5
        )
    )

  ##########################################################
  # Save figure
  ##########################################################

  pdf_file <- paste0(sample_name, "_Normalized.reads.pdf")

  ggsave(
    filename = pdf_file,
    plot = p,
    width = 13.5,
    height = 45,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )

  message("Saved: ", pdf_file)
}

############################################################
# Finished
############################################################

message("All genome coverage plots generated successfully.")
