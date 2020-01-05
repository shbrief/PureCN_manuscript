#!/usr/bin/env Rscript

# figure styling
source("~/data2/PureCN_manuscript/Figures/manuscript_figure_style.R")

if (plot_for == "purity") {
  customPlot = c(customPlot, 
                 scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2)),
                 scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.2)))
} else if (plot_for == "ploidy") {
  customPlot = c(customPlot,
                 scale_x_continuous(limits = c(0, 7), breaks = seq(0,7,1)),
                 scale_y_continuous(limits = c(0, 7), breaks = seq(0,7,1)))  
} else {
  customPlot = customPlot
}

# plotting
res_plot = ggplot(df) +
  geom_point(aes(x = df[,x], y = df[,y]),
             col = "black",
             alpha = 1/3) + 
  customLabs +
  geom_abline(col = "grey", lwd = 0.5) + 
  customPlot

# add concordance
cor = round(cor(df[,x], df[,y], use = "complete.obs"), 2)
final_plot = ggdraw(add_sub(res_plot, paste("Concordance\ncorrelation =", cor), 
                          vpadding = grid::unit(0, "lines"),
                          y = 11, x = 0.05,
                          hjust = 0, 
                          fontface = "italic",
                          size = 12))

# print(final_plot)