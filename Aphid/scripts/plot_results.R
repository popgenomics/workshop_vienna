library(tidyverse)
library(argparse)

## function ##
plot_aphid <- function(aphid_results, plot_out){
  df = as_tibble(read.csv(aphid_results))
  p_noevent = df %>% filter(topology == "((AB)C)") %>% pull(contribution_noevent)
  p_ILS_ABC = df %>% filter(topology == "((AB)C)") %>% pull(contribution_ILS)
  p_ILS_ACB = df %>% filter(topology == "((AC)B)") %>% pull(contribution_ILS)
  p_ILS_BCA = df %>% filter(topology == "((BC)A)") %>% pull(contribution_ILS)
  p_GF_ABC = df %>% filter(topology == "((AB)C)") %>% pull(contribution_GF)
  p_GF_ACB = df %>% filter(topology == "((AC)B)") %>% pull(contribution_GF)
  p_GF_BCA = df %>% filter(topology == "((BC)A)") %>% pull(contribution_GF)
  
  df_values = as_tibble(data.frame(
    p_noevent = p_noevent,
    p_ILS_ABC = p_ILS_ABC,
    p_ILS_ACB = p_ILS_ACB,
    p_ILS_BCA = p_ILS_BCA,
    p_GF_ABC = p_GF_ABC,
    p_GF_ACB = p_GF_ACB,
    p_GF_BCA = p_GF_BCA
  ))
  df_plot  = df_values %>% pivot_longer(cols = c(p_noevent, p_ILS_ABC, p_ILS_ACB, p_ILS_BCA, p_GF_ABC, p_GF_ACB, p_GF_BCA),
                                              names_to = "variable",
                                              values_to = "contribution")
  max_y = max(df_plot$contribution)
  limit_y = ceiling(max_y * 10) / 10
  if (limit_y <= 0.4){
    break_y = 0.05
  } else {
    break_y = 0.1
  }
  
  
  fig = ggplot(df_plot, aes(x = variable, y = contribution, color = variable)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) + 
    scale_color_manual(
      values = c("p_noevent" = "#fc8d59", 'p_GF_ABC' = "#00bfdd", 'p_GF_ACB' = "#55bfdb", 'p_GF_BCA' = "#99bfdb", "p_ILS_ABC" = "#00ef8b", "p_ILS_ACB" = "#55ef8b", "p_ILS_BCA" = "#99ef8b"),
      labels = c('Gene flow (AB)C)', 'Gene flow (AC)B)', 'Gene flow (BC)A)', 'ILS (AB)C)', 'ILS (AC)B)', 'ILS (BC)A)', 'No event')
    ) + 
    labs(
      x = "",
      y = "Contribution",
      color = "") +
    scale_x_discrete(
      labels = c('Gene flow (AB)C)', 'Gene flow (AC)B)', 'Gene flow (BC)A)', 'ILS (AB)C)', 'ILS (AC)B)', 'ILS (BC)A)', 'No event')
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = -45, hjust = 0)) + scale_y_continuous(limits = c(0, limit_y), breaks = seq(0, limit_y, break_y))
  
  ggsave(filename = plot_out, plot = fig, width = 10, height = 5)
}

parser <- ArgumentParser(description = "Processing results and generating graphs")
parser$add_argument(
  "-r", "--results",
  type = "character",
  required = TRUE,
  help = "Path to the processed Aphid file"
)

parser$add_argument(
  "-p", "--plot_out",
  type = "character",
  help = "Path where the plot will be save"
)
args <- parser$parse_args()

processed = args$results
graph = args$plot_out

plot_aphid(aphid_results = processed, plot_out = graph)