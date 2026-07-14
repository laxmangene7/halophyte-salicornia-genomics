
############################################################
# Full chromosome population genomics comparison


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
# LOAD TAJIMA'S D
############################################################

tajima_france <- fread(
  "S_sp.France.sites.with.S.France.Ref_50kb.window_All.chr.Tajima.D",
  sep = "\t",
  header = TRUE,
  data.table = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia sp France"
  )

tajima_proc <- fread(
  "S_procumbens.sites.with.S.France.Ref_50kb.window_All.chr.Tajima.D",
  sep = "\t",
  header = TRUE,
  data.table = FALSE
) %>%
  mutate(
    POS_MB = BIN_START / 1e6,
    Dataset = "Salicornia procumbens"
  )

tajima_all <- bind_rows(
  tajima_france,
  tajima_proc
)



############################################################
# LOAD PI
############################################################

pi_france <- fread(
  "S_sp.France.sites.with.S.France.Ref_50kb.window_All.chr.pi",
  sep = "\t",
  header = TRUE,
  data.table = FALSE
) %>%
  mutate(
    POS_MB = BIN_END / 1e6,
    Dataset = "Salicornia sp France"
  )

pi_proc <- fread(
  "S_procumbens.sites.with.S.France.Ref_50kb.window_All.chr.pi",
  sep = "\t",
  header = TRUE,
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
    !is.na(PI)
  )

############################################################
# CHROMOSOME GROUPS
############################################################

chrA <- paste0(
  "Chr",
  1:9,
  "A"
)

chrC <- paste0(
  "Chr",
  1:9,
  "C"
)

############################################################
# FUNCTION TO GENERATE PANELS
############################################################

make_population_plot <- function(chr_set, output_prefix){
  
  ##########################################################
  # Filter chromosomes
  ##########################################################
  
  tajima_subset <- tajima_all %>%
    filter(CHROM %in% chr_set)
  
  pi_subset <- pi_all %>%
    filter(CHROM %in% chr_set)
  
  ##########################################################
  # Tajima's D plot
  ##########################################################
  
  p1 <- ggplot(
    tajima_subset,
    aes(
      x = POS_MB,
      y = TajimaD,
      color = Dataset
    )
  ) +
    
    geom_point(
      alpha = 0.5,
      size = 0.8
    ) +
    
    geom_smooth(
      method = "loess",
      span = 0.20,
      se = FALSE,
      linewidth = 0.8
    ) +
    
    scale_color_manual(
      values = plot_colors
    ) +
    
    facet_wrap(
      ~CHROM,
      scales = "free_x",
      ncol = 3
    ) +
    
    labs(
      title = "Tajima's D",
      x = "",
      y = "Tajima's D"
    ) +
    
    theme_bw(base_size = 16) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 20,
        hjust = 0.5
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = 16
      ),
      
      axis.text = element_text(
        size = 10
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      
      legend.position = "top",
      
      legend.title = element_blank(),
      
      legend.text = element_text(
        size = 12
      )
    )
  
  ##########################################################
  # PI plot
  ##########################################################
  
  p2 <- ggplot(
    pi_subset,
    aes(
      x = POS_MB,
      y = PI,
      color = Dataset
    )
  ) +
    
    geom_point(
      alpha = 0.5,
      size = 0.8
    ) +
    
    geom_smooth(
      method = "loess",
      span = 0.20,
      se = FALSE,
      linewidth = 0.8
    ) +
    
    scale_color_manual(
      values = plot_colors
    ) +
    
    facet_wrap(
      ~CHROM,
      scales = "free_x",
      ncol = 3
    ) +
    
    labs(
      title = expression(
        Nucleotide~diversity~(pi)
      ),
      x = "Position (Mb)",
      y = expression(pi)
    ) +
    
    theme_bw(base_size = 16) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 20,
        hjust = 0.5
      ),
      
      axis.title = element_text(
        face = "bold",
        size = 16
      ),
      
      axis.text = element_text(
        size = 10
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      
      legend.position = "none"
    )
  
  ##########################################################
  # Combine panels
  ##########################################################
  
  combined_plot <-
    p1 /
    p2 +
    plot_layout(
      heights = c(1,1)
    )
  
  ##########################################################
  # Save PDF
  ##########################################################
  
  ggsave(
    paste0(output_prefix, ".pdf"),
    combined_plot,
    width = 18,
    height = 16,
    dpi = 600
  )
  
  ##########################################################
  # Save SVG
  ##########################################################
  
  ggsave(
    paste0(output_prefix, ".svg"),
    combined_plot,
    width = 18,
    height = 16
  )
  
  ##########################################################
  # Save PNG
  ##########################################################
  
  ggsave(
    paste0(output_prefix, ".png"),
    combined_plot,
    width = 18,
    height = 16,
    dpi = 600
  )
  
  ##########################################################
  # Print
  ##########################################################
  
  print(combined_plot)
}

############################################################
# GENERATE FIGURE 1 (A genomes)
############################################################

make_population_plot(
  chrA,
  "S_sp.France.and_S.procumbens_50kb.window_A-subgenome.Chr1A-Chr9A.Tajima.D_and.Pi"
)

############################################################
# GENERATE FIGURE 2 (C genomes)
############################################################

make_population_plot(
  chrC,
  "S_sp.France.and_S.procumbens_50kb.window_C-subgenome.Chr1C-Chr9C.Tajima.D_and.Pi"
)









