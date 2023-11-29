# my_cowplot_tweaks.R

my_theme <- theme_update(
  legend.position = "bottom",
  legend.box = "vertical",
  legend.box.just = "left",
  legend.title = element_text(size = 18),
  legend.text = element_text(size = 16),
  strip.text = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.title = element_text(size = 16)
)
