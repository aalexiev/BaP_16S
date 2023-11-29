# my_cowplot_tweaks.R

my_theme <- theme_update(
  legend.position = "top",
  legend.box = "vertical",
  legend.box.just = "left",
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 9),
  strip.text = element_text(size = 10),
  plot.caption = element_text(hjust = 0, size = 10),
  axis.text = element_text(size = 8)
)
