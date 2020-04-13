## For plotting hits and reads between given coorrdinates

library(tidyverse)

source('R001_common.R')


## PARAMETERS -----------------------------------------------------------------

## Source data
rcmp_file <- files$rcmp_Nm

# object of parameters for plotting
hit_plot_params <- list(
  bp_per_cm = 482,                 # scale bp per cm
  max_scale = 10,                  # maximum hight of graphed bars in adjusted log2 units
  log2_units_per_cm = 24,          # log2 units per cm (vertical scale for graphed bars) was 18
  sample_spacing_log2_units = 16,  # spacing between graphed samples (in log2 units)     was 20
  cm_for_arrows = 0.5,             # space for the gene arrows
  hit_lwd = 0.5                    # lwd (relative line width) for the hit lines
)

hit_plot_sample_params <- list(
  samples = c('MRg1a',
              'T1_2h',
              'T1_4h',
              'T1_6h',
              'T1_8h',
              'AyL4'),
  sample_cols = c('#1D91C0',
                  '#ADDD8E', '#78C679', '#41AB5D', '#238443',
                  '#000000')
)


## FUNCTIONS ------------------------------------------------------------------

plot_hits <- function(nuc_fr, nuc_to, reverse_dir = FALSE, # replicon = 'ADP1',
                      save_pdf_or_svg_file = NULL,
                      plot_params = parent.frame()$hit_plot_params,
                      sample_params = parent.frame()$hit_plot_sample_params) {
  
  # Initial setup
  n_samples <- length(sample_params$samples)
  
  plot_cm_vert <- ((n_samples - 1) * plot_params$sample_spacing_log2_units +
                     2 * plot_params$max_scale) / plot_params$log2_units_per_cm +
    plot_params$cm_for_arrows
  
  plot_cm_horiz <- (nuc_to - nuc_fr + 1) / plot_params$bp_per_cm
  
  line_y_vals <- plot_cm_vert - plot_params$max_scale / plot_params$log2_units_per_cm -
    ((seq(n_samples) - 1) * plot_params$sample_spacing_log2_units / plot_params$log2_units_per_cm)
  
  plot_width_in <- (plot_cm_horiz * 1.08 / 2.54) + par()$mai[2] + par()$mai[4] + par()$omi[2] + par()$omi[4]
  plot_height_in <- (plot_cm_vert * 1.08 / 2.54) + par()$mai[1] + par()$mai[3] + par()$omi[1] + par()$omi[3]
  
  if (!is.null(save_pdf_or_svg_file)) {
    if (endsWith(save_pdf_or_svg_file, ".pdf")) {
      pdf(file = save_pdf_or_svg_file,
          width = plot_width_in,
          height = plot_height_in,
          pointsize = 12, bg = "transparent")
    } else if (endsWith(save_pdf_or_svg_file, ".svg")) {
      svg(filename = save_pdf_or_svg_file,
          width = plot_width_in,
          height = plot_height_in,
          pointsize = 12, bg = "transparent")
    } else {
      return("Ambiguous save file suffix")
    }
  }
  
  # Make subset table of the hits in the region in relative cm measurements
  sub_rcmp_rel_cms <- rcmp_al2 %>%
    filter(
      # Replicon == replicon,
      Sample %in% c(sample_params$samples),
      EffPos >= nuc_fr,
      EffPos <= nuc_to
    ) %>%
    mutate(
      rel_ep = if(reverse_dir) {
        nuc_to - EffPos
      } else {
        EffPos - nuc_fr
      },
      rel_cm = rel_ep / plot_params$bp_per_cm,
      al2_cm = al2r / plot_params$log2_units_per_cm *
        if (reverse_dir) -1 else 1
    ) %>%
    select(
      Sample, rel_cm, al2_cm
    )
  
  # Create the plot
  plot(-1, -1, xlim = c(0, plot_cm_horiz), ylim = c(0, plot_cm_vert),
       bg = "transparent", xaxt = "n", yaxt = "n", ann = FALSE, bty = 'n')
  
  # Plot horiz line and hits for each sample
  for (sn in 1:n_samples) {
    # Horizontal lines
    lines(c(0, plot_cm_horiz), rep(line_y_vals[sn], 2), col = sample_params$sample_cols[sn],
          lwd = 0.75)
    
    # Make table of the relative cm measurements for the sample
    sample_hits_rel_cms <- sub_rcmp_rel_cms %>%
      filter(
        Sample %in% c(sample_params$samples[sn])
      )
    
    # Plot the hits (vertical bars)
    for (hn in 1:nrow(sample_hits_rel_cms)) {
      lines(rep(sample_hits_rel_cms$rel_cm[hn], 2),
            line_y_vals[sn] + c(0, sample_hits_rel_cms$al2_cm[hn]),
            col = sample_params$sample_cols[sn],
            lwd = plot_params$hit_lwd)
    }
  }
  
  # Plot annotation arrows
  #  make sub-annotation table of the region in relative cm measurements
  reg_locs_rel_cms <- annot %>%
    filter(
      # Include all loci that even partially cover (or completely span) the region
      # Replicon == replicon &
      (From <= nuc_to & To >= nuc_fr) |
        (From < nuc_fr & To > nuc_to)
    ) %>%
    mutate(
      fr_bp = if (reverse_dir) {
        nuc_to - ifelse(Strand == "+", From, To)
      } else {
        ifelse(Strand == "+", From, To) - nuc_fr
      },
      fr_cm = fr_bp / plot_params$bp_per_cm,
      to_bp = if (reverse_dir) {
        nuc_to - ifelse(Strand == "+", To, From)
      } else {
        ifelse(Strand == "+", To, From) - nuc_fr
      },
      to_cm = to_bp / plot_params$bp_per_cm
    ) %>%
    select(
      Label, fr_cm, to_cm
    )
  
  arrow_y_val <- plot_params$cm_for_arrows / 2
  
  for (an in 1:nrow(reg_locs_rel_cms)) {
    # draw arrow showing gene coordinates
    arrows(reg_locs_rel_cms$fr_cm[an], arrow_y_val, x1 = reg_locs_rel_cms$to_cm[an],
           length = 0.1, lwd = 1)
    # add label to arrow
    label_center_cm = mean(c(max(0, min(reg_locs_rel_cms$fr_cm[an],
                                        reg_locs_rel_cms$to_cm[an])),
                             min(plot_cm_horiz, max(reg_locs_rel_cms$fr_cm[an],
                                                    reg_locs_rel_cms$to_cm[an]))))  # places within plot area
    text(label_center_cm, 0, reg_locs_rel_cms$Label[an],
         cex = 0.75)
  }
  
  if (!is.null(save_pdf_or_svg_file)) {
    dev.off()
  }
}


## GET DATA -------------------------------------------------------------------

rcmp_tidy <- read_tsv(rcmp_file)




## PROGRAM --------------------------------------------------------------------

## Make new table of adjusted Log2 read values, signed by dir

colnames(rcmp_tidy)[startsWith(colnames(rcmp_tidy), "reads")] <- "reads"

minnz <- min(rcmp_tidy$reads[rcmp_tidy$reads > 0])

rcmp_al2 <- mutate(
  rcmp_tidy,
  al2r = log2(reads * 2 / minnz) * ifelse(Dir == 'R', -1, 1)
) %>% select(
  -reads,
  -Dir
)
rm(minnz, rcmp_tidy)

## Plots the hits at each locus of interest

plot_hits(705395, 710824, reverse_dir = TRUE, save_pdf_or_svg_file = '070_plots/locus_hits_nrdAB.pdf')
plot_hits(415991, 418524, reverse_dir = TRUE, save_pdf_or_svg_file = '070_plots/locus_hits_lnt.pdf')

