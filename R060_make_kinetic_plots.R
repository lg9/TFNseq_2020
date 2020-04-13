## For plotting gene depletions


library(tidyverse)
library(ggrepel)

source('R001_common.R')



## FUNCTIONS -------------------------------------------------------------------
make_panel <- function(panel_def_file, pdat = panel_plot_dat, trial = "T1",
                       show_labels = TRUE, save_file = NULL,
                       plot_params = parent.frame()$plot_params, seed = NULL) {
  if (trial == "T1") {
    trial_samples <- c('MRg1a', 'T1_2h', 'T1_4h', 'T1_6h', 'T1_8h', 'AyL4')
  } else if (trial == "T2") {
    trial_samples <- c('MRg1a', 'T2_2h', 'T2_4h', 'T2_6h', 'T2_8h', 'AyL3')
  } else {
    return("Unable to define trial samples / times")
  }
  if (!is.null(seed)) { set.seed(seed) }
  # get panel data
  panel_def_dat <- read_tsv(paste(plot_params$panel_files_path, panel_def_file, sep = ""))
  # make panel plot data tibble - w required transformed data for plotting
  ppdat <- pdat %>%
    filter(
      Locus %in% panel_def_dat$Locus,
      Sample %in% trial_samples
    ) %>%
    left_join(
      panel_def_dat,
      by = "Locus"
    ) %>%
    group_by(
      Group, Sample
    ) %>%
    mutate(
      n = n(),
      avg_lrpk = mean(lrpk)
    ) %>%
    select(
      -Locus, -lrpk
    ) %>%
    distinct(
      Group, Sample,
      .keep_all = TRUE
    ) %>%
    arrange(
      Group
    ) %>%
    group_by(
      Group
    ) %>%
    # assign label points (the last one still in the plot)
    mutate(
      in_plot = avg_lrpk >= plot_params$ylim[1] & 
        time <= plot_params$xlim[2]
    ) %>%
    mutate(
      label_in_plot = 
        ifelse(in_plot & (cumsum(in_plot) == sum(in_plot)), Label, "")
    ) %>%
    select(
      -in_plot
    )
  
  # set up named color list
  col_tbl <- ppdat %>%
    group_by(
      Group
    ) %>%
    summarise(
      col = first(Color)
    )
  group_colors <- col_tbl$col
  names(group_colors) <- col_tbl$Group
  
  # set up line size list
  sz_tbl <- ppdat %>%
    group_by(
      Group
    ) %>%
    summarize(
      sz = ifelse(first(n) == 1,
                  plot_params$line_sz_gene, plot_params$line_sz_group)
    )
  group_sizes <- sz_tbl$sz
  names(group_sizes) <- sz_tbl$Group
  
  # plot the data
  p <- ppdat %>%
    filter(
      !(is.infinite(avg_lrpk))
    ) %>%
    ggplot(
      aes(time, avg_lrpk, color = factor(Group))
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "Black", fill = NA)
    ) +
    { if (!is.null(save_file)) {
      theme(
        axis.title = element_text(size = plot_params$axis_title_size),
        axis.text = element_text(size = plot_params$axis_text_size)
      ) } 
    } +
    coord_cartesian(
      xlim = plot_params$xlim,
      ylim = plot_params$ylim
    ) +
    scale_x_continuous(
      breaks = c(0, 2, 4, 6, 8, 10, 12)
    ) +
    labs(
      x = "hours",
      y = "relative recovery (log2)"  #y = expression(Delta~"log reads per gene length"))
    ) +
    geom_line(
      aes(size = factor(Group)),
      alpha = plot_params$alpha,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = group_colors
    ) +
    scale_size_manual(
      values = group_sizes
    )
  
  if (show_labels) {
    p <- p +
      geom_text_repel(
        aes(label = label_in_plot),
        size = ifelse(!is.null(save_file), plot_params$label_text_size, NA),
        show.legend = FALSE)
  }
    
  if (!is.null(save_file)) {
    ggsave(
      save_file,
      plot = p,
      path = plot_params$save_plots_path,
      width = plot_params$save_width,
      height = plot_params$save_height,
      units = plot_params$save_units )
  }
  p
}

## LOAD DATA ------------------------------------------------------------------

TBGpkNN <- read_tsv(files$TBGpkNmNp)

## MAKE data table for plotting -----------------------------------------------

panel_plot_dat <- TBGpkNN %>%
  left_join(
    outgr_times,
    by = "Sample"
  ) %>%
  rename(
    rpk = starts_with("rpk")
  ) %>%
  mutate(
    lrpk = log2(rpk)
  ) %>%
  select(
    Locus, Sample, time, lrpk
  )

## Make plots of depletion kinetic (for Fig 3) --------------------------------

# object of parameters for the plotting function; these are for export for printing
plot_params <- list(
  trial = 'T1',
  panel_files_path = '070_plots/plot_def_files/',
  xlim = c(0, 8),
  ylim = c(-1.2, 0.4),
  line_sz_gene = 0.5,
  line_sz_group = 1.8,
  alpha = 0.85,
  save_plots_path = '070_plots/',
  save_width = 90,  #86
  save_height = 56,  #56
  save_units = "mm",
  label_text_size = 7 * 25.4 / 72.27,     # e.g., 7 points expressed in mm
  axis_title_size = 9,
  axis_text_size = 8
)



#make_panel('General_processes.xls', seed = 122)
make_panel('General_processes.xls', show_labels = FALSE, save_file = "Kinetic_plot_general_processes.pdf")

#make_panel('Nucleotide_synthesis.xls')
make_panel('Nucleotide_synthesis.xls', show_labels = FALSE, save_file = "Kinetic_plot_nucleotide_synthesis.pdf")

#make_panel('DNA_replication.xls')
make_panel('DNA_replication.xls', show_labels = FALSE, save_file = "Kinetic_plot_DNA_replication.pdf")

#make_panel('Trxl_regulation.xls')
make_panel('Trxl_regulation.xls', show_labels = FALSE, save_file = "Kinetic_plot_trxl_regulation.pdf")



