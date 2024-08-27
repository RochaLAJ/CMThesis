require(tidyverse)
require(dplyr)
require(ggplot2)
require(cowplot)
require(RColorBrewer)

list_of_files <- lapply(list.files(pattern=".dat"), FUN=read.table)
titles <- c("Wildtype", "p.Met918Val", "p.Met918Thr")

for (i in seq_along(list_of_files)) {
  colnames(list_of_files[[i]]) <- c("Frames", "Phi215", "Psi215", "Phi216", "Psi216")
}

create_plot <- function(file, colour, title) {
  ggplot() +
    geom_point(aes(x = Phi215, y = Psi215, colour=Psi215), shape=21, size = 2, stroke=1, fill=colour,data = file) +
    geom_point(aes(x = Phi216, y = Psi216, colour=Psi216), shape=21, size = 2, stroke=1, fill=colour, data = file) + xlab(expression(phi)) + ylab(expression(psi))  +
    scale_colour_gradient(low = colour, high = "black") + 
    xlim(-180, 180) + ylim(-180, 180) + theme_bw() + 
    ggtitle(title) +
    labs(color='Angles')
}

plots_list <- mapply(create_plot, file = list_of_files, colour = c("darkgreen", "orange", "red"), title = titles, SIMPLIFY = FALSE)

combined_plots <- plot_grid(plotlist = plots_list, nrow = 1, ncol = 3, rel_heights = c(1,1))

concise_title <- ggdraw() +
  draw_label("PhiPsi Relation Between Residues 913 and 914", size = 21)


final_plot <- plot_grid(concise_title, combined_plots, nrow = 2, rel_heights = c(0.1, 1))

print(final_plot)
ggsave("PhiPsi918v2.png", dpi=600, width=10, height=7, bg="white")
