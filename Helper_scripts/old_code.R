# old_code.R

for (f in list.files(path = "Chunks", pattern = "^2[346789]|^3[0"))

### gutCheck-alpha-by-plate
alpha.plate.dt <- merge(
  alpha.base,
  sample.dt[, .(Sequencing.ID, Sample, Expt_plateID, Extract_batch)],
  by = "Sequencing.ID"
) %>%
  melt(
    measure.vars = alpha.methods,
    variable.name = "Alpha",
    value.name = "Score"
  )
{
  ggplot(alpha.plate.dt, aes(x = Expt_plateID, y = Score)) +
    geom_quasirandom() +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", color = "red") +
    facet_wrap(~ Alpha, scales = "free_y")
  } %>%
  ggsave(
    filename = file.path(tomDir, "beeswarms_alphaDiv_by_plate.pdf"),
    width = 12,
    height = 12
  )

res.dt <- NULL
for (plate in unique(sample.dt$Expt_plateID)) {
  n.taxa <- subset_samples(ps.decontam, Expt_plateID == plate) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    ntaxa()
  res.dt <- rbind(
    res.dt,
    data.table(
      Plate = plate,
      nTaxa = n.taxa
    )
  )
}
res.dt

chao1gt100.dt <- sample.dt[Sample %in% alpha.plate.dt[Alpha == "Chao1" & Score > 100]$Sample]
chao1gt100.dt %>%
  saveRDS(file = file.path(tomDir, "dt_chao1_gt_100_samples.rds"))

for (col in names(chao1gt100.dt)) {
  print(chao1gt100.dt[, .(.N), keyby = col])
}

alpha.plate.dt <- merge(
  alpha.base,
  sample.dt,
  by = "Sequencing.ID"
) %>%
  melt(
    measure.vars = alpha.methods,
    variable.name = "Alpha",
    value.name = "Score"
  )

{
  ggplot(alpha.plate.dt, aes(x = as.numeric(Extract_batch), y = Score)) +
    geom_quasirandom() +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", color = "red") +
    facet_wrap(~ Alpha, scales = "free_y")
  } %>%
  ggsave(
    filename = file.path(tomDir, "beeswarms_alphaDiv_by_extractionBatch.pdf"),
    width = 14,
    height = 12
  )
######

# beta-lpr-stats-bapOnly-for-Tom
beta.lpr.bapOnly.res.file <- file.path(saveDir, "list_beta_div_lpr_bapOnly_results.rds")
beta.lpr.bapOnly.res <- redo.if("beta.stats", beta.lpr.bapOnly.res.file, {
  lapply(beta.methods, function(beta) {
    lapply(lpr.cycles, function(cycle) {
      lapply(lpr.epochs, function(epoch) {
        dist.mat <- dist.list[[beta]] %>%
          dist_subset(rownames(beta.lpr.dfs[[cycle]][[epoch]]))
        dbrda.full <- capscale(
          dist.mat ~ BaP_cat,
          data = beta.lpr.bapCat.dfs[[cycle]][[epoch]]
        )
        return(dbrda.full)
      })
    })
  })
})
beta.lpr.bapOnly.aovs.file <- file.path(
  saveDir,
  "list_beta_div_lpr_bapOnly_permanovas.rds"
)
beta.lpr.bapOnly.aovs <- redo.if("beta.stats", beta.lpr.bapOnly.aovs.file, {
  lapply(beta.methods, function(beta) {
    lapply(lpr.cycles, function(cycle) {
      lapply(lpr.epochs, function(epoch) {
        anova(beta.lpr.bapOnly.res[[beta]][[cycle]][[epoch]], by = "term") %>%
          tidy() %>%
          as.data.table() %>%
          return()
      })
    })
  })
})

sig.lpr.aov.bapOnly.id.dt <- lapply(beta.methods, function(beta) {
  lapply(lpr.cycles, function(cycle) {
    lapply(lpr.epochs, function(epoch) {
      aov <- beta.lpr.bapOnly.aovs[[beta]][[cycle]][[epoch]]
      if (any(aov$p.value <= 0.05, na.rm = T)) {
        dt <- data.table(Method = beta, Cycle = cycle, Epoch = epoch)
      }
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()

print(sig.lpr.aov.bapOnly.id.dt)
######

# beta-lpr-bapCat-plot
merge.cols <- c("Variable", "Method", "Cycle", "Epoch", "Sig")
biplot.data.bapCat.lists.file <- file.path(
  saveDir,
  "list_lpr_beta_biplot_bapCat_data.rds"
)
biplot.data.bapCat.lists <- redo.if("beta.stats", biplot.data.bapCat.lists.file, {
  lapply(1:nrow(sig.lpr.aov.id.bapCat.dt), function(i) {
    ids <- sig.lpr.aov.id.bapCat.dt[i] %>% unlist()
    dbrda <- beta.lpr.bapCat.res[[ids[1]]][[ids[2]]][[ids[3]]]$Step
    biplot.data <- get.biplot.data(ps = ps.decontam, ord = dbrda)
    smpls.dt <- copy(biplot.data$sample.coords)
    smpls.dt[, names(ids) := as.list(ids)]
    names(smpls.dt)[2:3] <- c("Axis1", "Axis2")
    vctrs.dt0 <- biplot.data$vector.coords
    vctrs.dt0[, 2:3] <- vctrs.dt0[, 2:3] * biplot.data$coord.scale
    vctrs.dt0[, names(ids) := as.list(ids)]
    names(vctrs.dt0)[2:3] <- c("Axis1", "Axis2")
    vctrs.dt <- sig.lpr.aov.bapCat.mod.terms.dt[, ..merge.cols][
      vctrs.dt0, on = merge.cols[-5]
    ]
    if (any(is.na(vctrs.dt$Sig))) {
      aov.dt <- sig.lpr.aov.bapCat.mod.terms.dt[
        Method == ids[1] & Cycle == ids[2] & Epoch == ids[3]
      ]
      if (aov.dt[str_detect(Variable, ":")]$Sig) {
        vctrs.dt[str_detect(Variable, ":")]$Sig <- TRUE
      }  else {
        vctrs.dt[str_detect(Variable, ":")]$Sig <- FALSE
      }
      if (aov.dt[str_detect(Variable, "^BaP_cat")]$Sig) {
        vctrs.dt[str_detect(Variable, "^BaP_cat")]$Sig <- TRUE
      } else {
        vctrs.dt[str_detect(Variable, "^BaP_cat")]$Sig <- FALSE
      }
    }
    vctrs.dt[
      , Sig.lab := factor(
        ifelse(Sig, "p <= 0.05", "p > 0.05"),
        levels = c("p <= 0.05", "p > 0.05")
      )
    ]
    axes.labs.dt <- data.table(
      Axis1 = biplot.data$axes.labs[1],
      Axis2 = biplot.data$axes.labs[2]
    )
    axes.labs.dt[, names(ids) := as.list(ids)]
    return.list <- list(
      Samples = smpls.dt,
      Vectors = vctrs.dt,
      Axes.labs = axes.labs.dt,
      Coord.scale = biplot.data$coord.scale
    )
    return(return.list)
  })
})

plot.env <- new.env()
plot.list <- lapply(biplot.data.bapCat.lists, function(x) {
  smpls.dt <- copy(x$Samples)
  cat(smpls.dt$Method[1], sep = "\n")
  smpls.dt <- smpls.dt[
    lpr.auc.dt[, .(Sample, Cycle, Epoch, AUC.trans)],
    on = c("Origin==Sample", "Cycle", "Epoch"),
    nomatch = 0
  ]
  vctrs.dt <- copy(x$Vectors)
  axes.labs.dt <- copy(x$Axes.labs)
  coord.scale <- x$Coord.scale
  plot <- ggplot(smpls.dt, aes(x = Axis1, y = Axis2))
  if (any(str_detect(vctrs.dt$Variable, ":")) | nrow(vctrs.dt) > 1 & all(vctrs.dt$Sig)) {
    plot <- plot +
      geom_point(size = 4, aes(color = as.factor(BaP_uM))) +
      bap.color.scale +
      geom_point(size = 2.5, shape = 21, color = "white", aes(fill = AUC.trans)) +
      scale_fill_gradient(
        name = "AUC (Tukey-transformed)",
        low = "white",
        high = "black"
      ) +
      new_scale_color() +
      geom_segment(
        data = vctrs.dt,
        x = 0, y = 0,
        size = 1.5,
        aes(xend = Axis1 / coord.scale, yend = Axis2 / coord.scale, color = Sig.lab),
        arrow = arrow(length = unit(0.02, "npc"))
      ) +
      geom_label_repel(
        data = vctrs.dt,
        aes(
          x = Axis1 / coord.scale,
          y = Axis2 / coord.scale,
          label = Variable,
          color = Sig.lab
        ),
        fill = rgb(1, 1, 1, 0.4)
      ) +
      scale_color_manual("Term vector sig.", values = sig.colors)

    plot.env$full.legend <- get_legend(plot)
  } else {
    if (!("AUC.trans" %in% vctrs.dt$Variable)) {
      auc.term.check <- FALSE
    } else {
      auc.term.check <- vctrs.dt[Variable == "AUC.trans"]$Sig
    }
    if (auc.term.check) {
      plot <- plot +
        geom_point(size = 2.5, shape = 21, color = "black", aes(fill = AUC.trans)) +
        scale_fill_gradient(
          name = "AUC (Tukey-transformed)",
          low = "white",
          high = "black"
        )
    } else {
      plot <- plot + geom_point(size = 2.5, aes(color = as.factor(BaP_uM))) +
        bap.color.scale
    }
  }
  sig.terms <- vctrs.dt[(Sig)]$Variable
  nonsig.terms <- vctrs.dt[!(Sig)]$Variable
  if (length(nonsig.terms) == 0) {
    caption <- paste0(
      "All terms (", paste(sig.terms, collapse = ", "),
      ") were significant (p <= 0.05) in the optimized model"
    )
  } else {
    caption <- paste(
      paste(sig.terms, collapse = " and "),
      ifelse(length(sig.terms) > 1, "were", "was"),
      "signifcant (p <= 0.05) in the optimized model.",
      paste(nonsig.terms, collapse = " and "),
      ifelse(length(nonsig.terms) > 1, "were", "was"),
      "included in the optimized model, but not significant"
    )
  }
  plot <- plot +
    labs(
      title = smpls.dt$Method[1],
      subtitle = paste0(
        tools::toTitleCase(smpls.dt$Cycle[1]), " cycle (",
        str_replace(smpls.dt$Epoch[1], "E", "Epoch "), ")"
      ),
      x = axes.labs.dt$Axis1,
      y = axes.labs.dt$Axis2,
      caption = wrap.caption(caption, w = 80)
    ) +
    theme(legend.position = "none")
  return(plot)
})

plot.list[[length(plot.list) + 1]] <- plot.env$full.legend
plot_grid(plotlist = plot.list, nrow = 4, align = "h", axis = "tb") %>%
  ggsave(
    filename = file.path(plotDir, "ordinations_bapCat_lpr.png"),
    width = img.width * 5,
    height = img.height * 3,
    dpi = img.dpi
  )
######


# beta-lprs-for-Tom
{ plot.list[[11]] + theme(legend.position = "right") } %>%
  ggsave(
    filename = file.path(tomDir, "ord_LPR_Aitchison_E1_darkCycle.pdf"),
    width = 11,
    height = 10
  )
plot.list <- NULL
smpls.dt <- copy(biplot.data.lists[[13]]$Samples)
smpls.dt <- smpls.dt[
  lpr.auc.dt[, .(Sample, Cycle, Epoch, AUC)],
  on = c("Origin==Sample", "Cycle", "Epoch"),
  nomatch = 0
]
vctrs.dt <- copy(biplot.data.lists[[13]]$Vectors)
axes.labs.dt <- copy(biplot.data.lists[[13]]$Axes.labs)
base.plot <- ggplot(smpls.dt, aes(x = Axis1, y = Axis2))

plot.list[[1]] <- base.plot +
  geom_point(size = 2, aes(color = as.factor(BaP_uM))) +
  bap.color.scale +
  labs(
    title = smpls.dt$Method[1],
    subtitle = paste0(
      tools::toTitleCase(smpls.dt$Cycle[1]), " cycle (",
      str_replace(smpls.dt$Epoch[1], "E", "Epoch "), ")"
    ),
    x = axes.labs.dt$Axis1,
    y = axes.labs.dt$Axis2
  )

plot.list[[2]] <- base.plot +
  geom_point(size = 2, aes(color = AUC)) +
  scale_color_gradientn(
    name = "AUC (µm × s)",
    colors = auc.colors
  ) +
  theme(
    legend.text = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_text(vjust = 0.9)
  ) +
  labs(x = axes.labs.dt$Axis1, y = axes.labs.dt$Axis2)

plot_grid(plotlist = plot.list, nrow = 1, align = "h", axis = "l") %>%
  ggsave(
    filename = file.path(tomDir, "ords_LPR_WUnifrac_E1_darkCycle_byBaP_byAUC.pdf"),
    width = 22,
    height = 10
  )

# {
#   base.plot +
#     geom_point(size = 4, aes(color = as.factor(BaP_uM))) +
#     bap.color.scale +
#     geom_point(size = 2, shape = 21, color = "white", aes(fill = AUC)) +
#     scale_fill_gradient(
#       name = "AUC (µm × s)",
#       low = "white",
#       high = "black"
#     ) +
#     geom_segment(
#       data = vctrs.dt[(Sig)],
#       x = 0, y = 0,
#       aes(xend = Axis1 * 2, yend = Axis2 * 2),
#       color = "black",
#       arrow = arrow(length = unit(0.01, "npc"))
#     ) +
#     labs(
#       title = smpls.dt$Method[1],
#       subtitle = paste0(
#         tools::toTitleCase(smpls.dt$Cycle[1]), " cycle (",
#         str_replace(smpls.dt$Epoch[1], "E", "Epoch "), ")"
#       ),
#       x = axes.labs.dt$Axis1,
#       y = axes.labs.dt$Axis2
#     ) +
#     theme(legend.box = "horizontal", legend.position = "bottom")
#   }


smpls.dt <- copy(biplot.data.bapCat.lists[[15]]$Samples)
smpls.dt <- smpls.dt[
  lpr.auc.dt[, .(Sample, Cycle, Epoch, AUC)],
  on = c("Origin==Sample", "Cycle", "Epoch"),
  nomatch = 0
]
vctrs.dt <- copy(biplot.data.bapCat.lists[[15]]$Vectors)
vctrs.dt$Label <- str_replace_all(vctrs.dt$Variable, "BaP_cat", "BaP[") %>%
  str_replace_all("u", "µ") %>%
  paste0("]")
vctrs.dt[Variable == "AUC"]$Label <- "AUC:BaP[0µM]"
axes.labs.dt <- copy(biplot.data.bapCat.lists[[15]]$Axes.labs)
{
  ggplot(smpls.dt, aes(x = Axis1, y = Axis2)) +
    geom_point(size = 4, aes(color = as.factor(BaP_uM))) +
    bap.color.scale +
    geom_point(size = 2, shape = 21, color = "white", aes(fill = AUC)) +
    scale_fill_gradient(
      name = "AUC (µm × s)",
      low = "white",
      high = "black"
    ) +
    geom_segment(
      data = vctrs.dt[str_detect(Variable, "AUC")],
      x = 0, y = 0,
      aes(xend = Axis1 * 3, yend = Axis2 * 3),
      color = "black",
      arrow = arrow(length = unit(0.01, "npc"))
    ) +
    geom_label_repel(
      data = vctrs.dt[str_detect(Variable, "AUC")],
      aes(x = Axis1 * 3, y = Axis2 * 3, label = Label),
      color = "black",
      fill = rgb(1, 1, 1, 0.6)
    ) +
    labs(
      title = smpls.dt$Method[1],
      subtitle = paste0(
        tools::toTitleCase(smpls.dt$Cycle[1]), " cycle (",
        str_replace(smpls.dt$Epoch[1], "E", "Epoch "), ")"
      ),
      x = axes.labs.dt$Axis1,
      y = axes.labs.dt$Axis2
    ) +
    theme(legend.box = "horizontal", legend.position = "bottom")
} %>%
  ggsave(
    filename = file.path(tomDir, "ord_LPR_UUnifrac_E2_darkCycle_BaPxAUC.pdf"),
    width = 11,
    height = 10
  )

base.plot <- ggplot(smpls.dt, aes(x = Axis1, y = Axis2))
plot.list <- NULL
plot.list[[1]] <- base.plot +
  geom_point(size = 2, aes(color = as.factor(BaP_uM))) +
  bap.color.scale +
  labs(
    title = smpls.dt$Method[1],
    subtitle = paste0(
      tools::toTitleCase(smpls.dt$Cycle[1]), " cycle (",
      str_replace(smpls.dt$Epoch[1], "E", "Epoch "), ")"
    ),
    x = axes.labs.dt$Axis1,
    y = axes.labs.dt$Axis2
  )

plot.list[[2]] <- base.plot +
  geom_point(size = 2, aes(color = AUC)) +
  scale_color_gradientn(
    name = "AUC (µm × s)",
    colors = auc.colors
  ) +
  theme(
    legend.text = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_text(vjust = 0.9)
  ) +
  labs(x = axes.labs.dt$Axis1, y = axes.labs.dt$Axis2)

plot_grid(plotlist = plot.list, nrow = 1, align = "h", axis = "l") %>%
  ggsave(
    filename = file.path(tomDir, "ords_LPR_UUnifrac_E2_darkCycle_byBaP_byAUC.pdf"),
    width = 22,
    height = 10
  )
######

# lpr-lms-sigIntrxns-dark-E1-plots-for-Tom, include=FALSE}
focal.taxa <- c("ASV00063", "ASV00116")
dark.e1.reg.dt <- copy(intrxn.reg.data)[
  Cycle == "dark" & Epoch == "E1" & Taxon %in% focal.taxa
]
impt.dt <- dark.e1.reg.dt[, .(Importance = Importance[1]), by = Taxon]
plot.taxa <- impt.dt[order(Importance, decreasing = T)]$Taxon
taxa.labs.dt <- taxa.dts.list$rarefied[
  Taxon %in% plot.taxa[str_detect(plot.taxa, "ASV")]
]
taxa.labs.dt[, Kingdom := NULL]
taxa.labs.dt[, Lab := paste(Family, Genus, Taxon, sep = "|")]
taxa.labs.dt[, Orig.ID := Taxon]
non.asv.taxa <- plot.taxa[!str_detect(plot.taxa, "ASV")]
if (length(non.asv.taxa) > 0) {
  non.asv.taxa.dt <- str_split(non.asv.taxa, "\\.") %>%
    lapply(function(x) { x <- c(x, rep(NA, 5 - length(x))) }) %>%
    lapply(function(x) { str_remove(x, "^[a-z]__") }) %>%
    do.call(rbind, .) %>%
    as.data.table() %>%
    set_colnames(c("Phylum", "Class", "Order", "Family", "Genus"))
  non.asv.taxa.dt[, Taxon := NA]
  non.asv.taxa.dt$Lab <- str_split(non.asv.taxa, "\\.") %>%
    lapply(tail, 2) %>%
    lapply(function(x) { str_remove(x, "^[a-z]__") }) %>%
    sapply(paste, collapse = "|")
  non.asv.taxa.dt[, Orig.ID := non.asv.taxa]
  taxa.labs.dt <- rbind(taxa.labs.dt, non.asv.taxa.dt)
}
setkey(taxa.labs.dt, Orig.ID)

plot.taxa.dt <- merge(
  dark.e1.reg.dt[Taxon %in% plot.taxa],
  taxa.labs.dt[, .(Orig.ID, Lab)],
  by.x = "Taxon",
  by.y = "Orig.ID"
)

keep.cols <- c("BaP_uM", "AUC", plot.taxa)
lpr.taxa.set <- best.lpr.rf.datasets[LPR.cycle == "dark" & LPR.epoch == "E1"] %>%
  use_series(Taxa.set)
lpr.abund.set <- best.lpr.rf.datasets[LPR.cycle == "dark" & LPR.epoch == "E1"] %>%
  use_series(Abundance.type)
plot.data <- rbind(
  copy(lpr.rf.data.list[[lpr.taxa.set]][[lpr.abund.set]]$dark$E1$Train),
  copy(lpr.rf.data.list[[lpr.taxa.set]][[lpr.abund.set]]$dark$E1$Test)
)[, ..keep.cols] %>%
  melt(
    id.vars = c("BaP_uM", "AUC"),
    variable.name = "Orig.ID",
    value.name = "Abund"
  ) %>%
  merge(taxa.labs.dt[, .(Orig.ID, Lab)], by = "Orig.ID")
plot.data[, Lab := factor(Lab, levels = taxa.labs.dt[plot.taxa]$Lab)]

plot.file <- file.path(
  tomDir,
  paste0(
    "scatterPlots_darkCycle_E1_",
    paste(focal.taxa, collapse = "-"),
    "_wRegressions.pdf"
  )
)
{
  ggplot(plot.data, aes(x = AUC, y = sqrt(Abund))) +
    geom_point(aes(color = as.factor(BaP_uM))) +
    geom_abline(
      data = plot.taxa.dt,
      aes(intercept = Intercept, slope = Slope, color = as.factor(BaP_uM))
    ) +
    facet_wrap(~ Lab, scales = "free_y", ncol = 4) +
    bap.color.scale +
    labs(
      caption = ,
      x = "AUC (µm × s)",
      y = "sqrt(Abundance)"
    )
} %>%
  ggsave(
    filename = plot.file,
    width = 16,
    height = 8
  )

plot.file <- file.path(
  tomDir,
  paste0(
    "beeswarm_darkCycle_E1_",
    paste(focal.taxa, collapse = "-"),
    ".pdf"
  )
)
{
  ggplot(plot.data, aes(x = BaP_uM, y = sqrt(Abund))) +
    geom_quasirandom(aes(color = as.factor(BaP_uM))) +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", color = "black") +
    bap.color.scale +
    facet_wrap(~ Lab, scales = "free_y", ncol = 4) +
    labs(
      caption = ,
      x = "BaP exposure (µM)",
      y = "sqrt(Abundance)"
    ) +
    theme(legend.position = "none")
} %>%
  ggsave(
    filename = plot.file,
    width = 16,
    height = 8
  )

######

#### Significant __interactions__ for LPR dark cycle (shared by epochs 1 & 2)

# lpr-lms-sigIntrxns-dark-E1&2-plots, include=FALSE}
dark.e1e2.reg.dt <- copy(intrxn.reg.data)[Cycle == "dark" & Epoch %in% c("E1", "E2")]
impt.dt <- dark.e1e2.reg.dt[
  , .(Importance = Importance[1]),
  by = .(Taxon, Epoch)
] %>%
  .[
    Taxon %in% names(table(.$Taxon)[table(.$Taxon) > 1]),
    .(Importance = mean(Importance)),
    by = Taxon
  ]
max.facets <- 16
for (i in seq(1, nrow(impt.dt), by = max.facets)) {
  cat(i, sep = "\n")
  end <- min(c(i + max.facets - 1, nrow(impt.dt)))
  plot.taxa <- impt.dt[order(Importance, decreasing = T)]$Taxon[i:end]
  taxa.labs.dt <- taxa.dts.list$rarefied[
    Taxon %in% plot.taxa[str_detect(plot.taxa, "ASV")]
  ]
  taxa.labs.dt[, Kingdom := NULL]
  taxa.labs.dt[, Lab := paste(Family, Genus, Taxon, sep = "|")]
  taxa.labs.dt[, Orig.ID := Taxon]
  non.asv.taxa <- plot.taxa[!str_detect(plot.taxa, "ASV")]
  if (length(non.asv.taxa) > 0) {
    non.asv.taxa.dt <- str_split(non.asv.taxa, "\\.") %>%
      lapply(function(x) { x <- c(x, rep(NA, 5 - length(x))) }) %>%
      lapply(function(x) { str_remove(x, "^[a-z]__") }) %>%
      do.call(rbind, .) %>%
      as.data.table() %>%
      set_colnames(c("Phylum", "Class", "Order", "Family", "Genus"))
    non.asv.taxa.dt[, Taxon := NA]
    non.asv.taxa.dt$Lab <- str_split(non.asv.taxa, "\\.") %>%
      lapply(tail, 2) %>%
      lapply(function(x) { str_remove(x, "^[a-z]__") }) %>%
      sapply(paste, collapse = "|")
    non.asv.taxa.dt[, Orig.ID := non.asv.taxa]
    taxa.labs.dt <- rbind(taxa.labs.dt, non.asv.taxa.dt)
  }
  setkey(taxa.labs.dt, Orig.ID)

  plot.taxa.dt <- merge(
    dark.e1.reg.dt[Taxon %in% plot.taxa],
    taxa.labs.dt[, .(Orig.ID, Lab)],
    by.x = "Taxon",
    by.y = "Orig.ID"
  )

  keep.cols <- c("BaP_uM", "AUC", plot.taxa)
  lpr.taxa.set <- best.lpr.rf.datasets[LPR.cycle == "dark" & LPR.epoch == "E1"] %>%
    use_series(Taxa.set)
  lpr.abund.set <- best.lpr.rf.datasets[LPR.cycle == "dark" & LPR.epoch == "E1"] %>%
    use_series(Abundance.type)
  plot.data <- rbind(
    copy(lpr.rf.data.list[[lpr.taxa.set]][[lpr.abund.set]]$dark$E1$Train),
    copy(lpr.rf.data.list[[lpr.taxa.set]][[lpr.abund.set]]$dark$E1$Test)
  )[, ..keep.cols] %>%
    melt(
      id.vars = c("BaP_uM", "AUC"),
      variable.name = "Orig.ID",
      value.name = "Abund"
    ) %>%
    merge(taxa.labs.dt[, .(Orig.ID, Lab)], by = "Orig.ID")
  plot.data[, Lab := factor(Lab, levels = taxa.labs.dt[plot.taxa]$Lab)]
  cap <- paste(
    "Important (by random forest) taxa", i, "-", end,
    "that had significant relationships with an interaction between",
    "LPR dark cycle epoch 1 AUC and BaP treatment"
  )
  p <- ggplot(plot.data, aes(x = AUC, y = sqrt(Abund))) +
    geom_point(aes(color = as.factor(BaP_uM))) +
    geom_abline(
      data = plot.taxa.dt,
      aes(intercept = Intercept, slope = Slope, color = as.factor(BaP_uM))
    ) +
    facet_wrap(~ Lab, scales = "free_y", ncol = 4) +
    bap.color.scale +
    labs(
      caption = ,
      x = "AUC (µm × s)",
      y = "sqrt(Abundance)"
    )
  plot.file <- file.path(
    plotDir,
    paste0(
      "scatterPlots_darkCycle_E1-E2shared_sigImptTaxa", i, "-", end,"_wRegressions.png"
    )
  )
  ggsave(
    p,
    filename = plot.file,
    width = 22,
    height = 11,
    dpi = img.dpi
  )
}
######

###
both.auc.dt <- merge(
  epr.auc.dt[, .(Sample, BaP_uM, AUC)],
  lpr.auc.dt[, .(Sample, Cycle, Epoch, AUC)],
  by = "Sample",
  all = TRUE,
  suffixes = c(".EPR", ".LPR")
)
ggplot(
  both.auc.dt[complete.cases(both.auc.dt)],
  aes(x = AUC.EPR, y = AUC.LPR, color = as.factor(BaP_uM))
) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, se = F) +
  bap.color.scale +
  facet_grid(Cycle ~ Epoch, scales = "free")

behaviorDir <- "Input/Behavior_results/GUI files/RemovedAffected"
epr.file <- file.path(
  behaviorDir,
  "behavior_24h_FOR_GUI_HMAT_With_Removal_10-23-2020 112541.csv"
)
epr.raw.dt <- read.csv(epr.file) %>% as.data.table()
epr.raw.dt[, remove := ifelse(is.na(remove), 1, remove)]
epr.removed <- epr.raw.dt[, .(EPR.removed = sum(remove)), by = .(conc)]

lpr.file <- file.path(behaviorDir, "behavior_5d_FOR_GUI_LPR_10-23-2020 112541.csv")
lpr.raw.dt <- read.csv(lpr.file) %>% as.data.table()
lpr.raw.dt[, remove := ifelse(is.na(t1), 1, 0)]
lpr.removed <- lpr.raw.dt[, .(LPR.removed = sum(remove)), by = .(conc)]

epr.removed[lpr.removed, on = "conc"]

ggplot(lpr.auc.dt, aes(x = AUC)) +
  geom_histogram(aes(fill = as.factor(BaP_uM)), position = "dodge", binwidth = 100) +
  bap.fill.scale +
  facet_grid(Cycle ~ Epoch, scales = "free")
######
<<<<<<< HEAD

### trying-stuff, eval=FALSE, echo=FALSE}
profile <- tweedie.profile(
  AUC ~ BaP_uM,
  data = lpr.auc.dt, p.vec = seq(1.1, 1.9, 0.1)
)

qqplot(
  x = qtweedie(
    p = ppoints(500),
    xi = profile$p.max,
    mu = mean(lpr.auc.dt$AUC),
    phi = profile$phi.max
  ),
  y = quantile(lpr.auc.dt$AUC, probs = ppoints(500))
)
qqline(
  y = quantile(lpr.auc.dt$AUC, probs = ppoints(500)),
  distribution = function(p) {
    qtweedie(p, xi = profile$p.max, mu = mean(lpr.auc.dt$AUC), phi = profile$phi.max)
  },
  probs = c(0.1, 0.6),
  col = 2
)

lpr.auc.dt$AUC.trans <- transformTukey(lpr.auc.dt$AUC, plotit = F, quiet = T)
profile <- tweedie.profile(
  AUC.trans ~ BaP_uM,
  data = lpr.auc.dt, p.vec = seq(1.1, 1.9, 0.1)
)

qqplot(
  x = qtweedie(
    p = ppoints(500),
    xi = profile$p.max,
    mu = mean(lpr.auc.dt$AUC.trans),
    phi = profile$phi.max
  ),
  y = quantile(lpr.auc.dt$AUC.trans, probs = ppoints(500))
)
qqline(
  y = quantile(lpr.auc.dt$AUC.trans, probs = ppoints(500)),
  distribution = function(p) {
    qtweedie(
      p,
      xi = profile$p.max,
      mu = mean(lpr.auc.dt$AUC.trans),
      phi = profile$phi.max
    )
  },
  probs = c(0.1, 0.6),
  col = 2
)
######
# Spearman correlations

### spearman-ASV-cors, include=FALSE, eval=FALSE}
taxon.pairs.file <- file.path(saveDir, "list_all_taxon_pairs_by_abundType.rds")
taxon.pairs.list <- redo.if("spearman.cors", taxon.pairs.file, {
  cl <- makeCluster(length(abundance.types), type = "FORK", outfile = "")
  registerDoParallel(cl, length(abundance.types))
  taxon.pairs <- foreach(
    abund.type = abundance.types,
    .final = abundType.setNames
  ) %dopar% {
    taxa <- taxa_names(ps.list[[abund.type]])
    pair.mat <- matrix(
      FALSE,
      nrow = length(taxa),
      ncol = length(taxa),
      dimnames = list(taxa, taxa)
    )
    pair.mat[upper.tri(pair.mat)] <- TRUE
    pair.df <- reshape2::melt(
      pair.mat,
      varnames = c("Taxon1", "Taxon2"),
      value.name = "Keep"
    )
    pair.dt <- as.data.table(pair.df[pair.df$Keep, ])
    pair.dt[, Keep := NULL]
    return(pair.dt)
  }
  stopCluster(cl)
  taxon.pairs
})

### the following labelled everythign "rarefied". Must be a scoping issue with foreach...
spearman.cors.file <- file.path(saveDir, "list_spearman_cors_by_abundType_and_bap.rds")
spearman.cors <- redo.if("spearman.cors", spearman.cors.file, {
  lapply(abundance.types, function(abund.type) {
    lapply(bap.concs, function(conc) {
      ps.sub <- ps.list[[abund.type]] %>%
        subset_samples(BaP_uM == conc) %>%
        prune_taxa(taxa_sums(.) > 0, .)
      sub.taxon.pairs <- taxon.pairs.list[[abund.type]][
        Taxon1 %in% taxa_names(ps.sub) & Taxon2 %in% taxa_names(ps.sub)
      ]
      asv.tbl <- otu.matrix(ps.sub)
      cat(paste(abund.type, conc, sep = "-"), sep = "\n")
      cl <- makeCluster(nCores, type = "FORK")
      registerDoParallel(cl, nCores)
      cor.res.dt <- foreach(
        i = 1:nrow(sub.taxon.pairs),
        .final = rbindlist,
        .verbose = TRUE
      ) %dopar% {
        t1 <- as.character(sub.taxon.pairs[i]$Taxon1)
        t2 <- as.character(sub.taxon.pairs[i]$Taxon2)
        res <- cor.test(asv.tbl[, t1], asv.tbl[, t2], method = "spearman") %>%
          suppressWarnings()
        res.dt <- data.table(
          Abund.type = abund.type,
          BaP_uM = conc,
          Taxon1 = t1,
          Taxon2 = t2,
          Rho = res$estimate,
          Pval = res$p.value
        )
        return(res.dt)
      }
      stopCluster(cl)
      return(cor.res.dt)
    }) %>% rbindlist()
  }) %>% rbindlist()
})
######

### spearman-genus-cors, include=FALSE, eval=FALSE}
genus.ps.file <- file.path(saveDir, "physeq_genus_glom.rds")
genus.ps <- redo.if("spearman.cors", genus.ps.file, {
  asv.taxa.dt <- taxa.data.table(ps.list$rarefied) %>% setkeyv("Taxon")
  tax_glom(ps.list$rarefied, taxrank = "Genus", NArm = T)
  taxa_names(genus.ps) <- asv.taxa.dt[taxa_names(genus.ps)]$Genus
  genus.ps
})

pct <- 0.50
top.genera <- prune_taxa(
  taxa_sums(genus.ps) >= quantile(taxa_sums(genus.ps), probs = pct),
  genus.ps
) %>%
  taxa_names()

genus.pairs.file <- file.path(saveDir, "dt_all_genus_pairs_by_abundType.rds")
genus.pairs.dt <- redo.if("spearman.cors", genus.pairs.file, {
  taxa <- taxa_names(genus.ps)
  pair.mat <- matrix(
    FALSE,
    nrow = length(taxa),
    ncol = length(taxa),
    dimnames = list(taxa, taxa)
  )
  pair.mat[upper.tri(pair.mat)] <- TRUE
  pair.df <- reshape2::melt(
    pair.mat,
    varnames = c("Taxon1", "Taxon2"),
    value.name = "Keep"
  )
  pair.dt <- as.data.table(pair.df[pair.df$Keep, ])
  pair.dt[, Keep := NULL]
  pair.dt
})


genus.cors.file <- file.path(saveDir, "dt_genus_spearman_cors_by_bap.rds")
genus.cors <- redo.if("spearman.cors", genus.cors.file, {
  lapply(bap.concs, function(conc) {
    cat(conc, sep = "\n")
    subset.expr <<- paste("BaP_uM ==", conc)
    ps.sub <- subset_samples(genus.ps, eval(parse(text = subset.expr))) %>%
      prune_taxa(taxa_sums(.) > 0, .)
    sub.taxon.pairs <- genus.pairs.dt[
      Taxon1 %in% taxa_names(ps.sub) & Taxon2 %in% taxa_names(ps.sub)
    ]
    asv.tbl <- otu.matrix(ps.sub)
    cl <- makeCluster(nCores, type = "FORK")
    registerDoParallel(cl, nCores)
    cor.res.dt <- foreach(
      i = 1:nrow(sub.taxon.pairs),
      .final = rbindlist,
      .verbose = FALSE
    ) %dopar% {
      t1 <- as.character(sub.taxon.pairs[i]$Taxon1)
      t2 <- as.character(sub.taxon.pairs[i]$Taxon2)
      res <- cor.test(asv.tbl[, t1], asv.tbl[, t2], method = "spearman") %>%
        suppressWarnings()
      res.dt <- data.table(
        BaP_uM = conc,
        Taxon1 = t1,
        Taxon2 = t2,
        Rho = res$estimate,
        Pval = res$p.value
      )
      return(res.dt)
    }
    stopCluster(cl)
    return(cor.res.dt)
  }) %>% rbindlist()
})
genus.cors[, Adj.pval := p.adjust(Pval, method = "BY")]
genus.cors.sig <- genus.cors[Adj.pval <= 0.01]
genus.cors.sig[, .(Count = .N), by = BaP_uM]
genus.cors.sig[Rho < 0]

{
  ggplot(genus.cors.sig, aes(x = Rho)) +
    geom_histogram(binwidth = 0.01)
} %>%
  ggsave(filename = file.path(tomDir, "histogram_genus_spearman_rhos.pdf"))

{
  ggplot(genus.cors.sig, aes(x = Rho)) +
    geom_histogram(binwidth = 0.01) +
    scale_y_log10() +
    labs(y = "log10(count)")
  } %>%
  ggsave(filename = file.path(tomDir, "histogram_genus_spearman_rhos_logCounts.pdf"))
######

### spearman-genus-plots, include=FALSE, eval=FALSE}
genus.taxonomy <- taxa.data.table(genus.ps)
genus.taxonomy[, Taxon := Genus]
genus.taxonomy[, Genus := NULL]
set.redo.false("pick.random.colors")
v.colors.dt.file <- file.path(saveDir, "dt_genus_vertex_phylum_colors.rds")
v.colors.dt <- redo.if("pick.random.colors", v.colors.dt.file, {
  data.table(
    Phylum = sort(unique(genus.taxonomy$Phylum)),
    Color = randomColor(
      count = length(unique(genus.taxonomy$Phylum)),
      luminosity = "bright"
    )
  ) %>% setkeyv("Phylum")
})

set2.colors <- brewer.pal(8, "Set2")
n.bins <- 30
rho.bins <- bin(genus.cors.sig$Rho, nbins = n.bins)
e.colors.dt <- data.table(
  Rho.bin = levels(rho.bins),
  Color = colorRampPalette(c(set2.colors[2], "grey80", set2.colors[1]))(n.bins)
) %>% setkeyv("Rho.bin")
e.colors.dt$Start <- str_remove_all(e.colors.dt$Rho.bin, "\\(|\\]") %>%
  str_split(",") %>%
  sapply(`[`, 1) %>%
  as.numeric()
e.colors.dt$End <- str_remove_all(e.colors.dt$Rho.bin, "\\(|\\]") %>%
  str_split(",") %>%
  sapply(tail, 1) %>%
  as.numeric()
# shared.taxa <- sapply(
#   unique(c(genus.cors.sig$Taxon1, genus.cors.sig$Taxon2)),
#   function(taxon) {
#     pres.count <- sapply(bap.concs, function(conc) {
#       subset.dt <- genus.cors.sig[BaP_uM == conc]
#       taxon %in% unique(c(subset.dt$Taxon1, subset.dt$Taxon2))
#     }) %>% sum()
#     if (pres.count == 4) { return(taxon) }
# }) # this does not reduce the number of taxa plotted

network.stats.dt <- data.table()
par(mfrow = c(2,2), mar = c(1,1,1,1))
for (conc in bap.concs) {
  subset.dt <- genus.cors.sig[BaP_uM == conc & abs(Rho) > 0.3]
  edges.dt <- subset.dt[, .(Taxon1, Taxon2, Rho)]
  vertices.dt <- data.table(
    Taxon = unique(c(subset.dt$Taxon1, subset.dt$Taxon2))
  ) %>%
    merge(genus.taxonomy, by = "Taxon")
  net <- graph_from_data_frame(d = edges.dt, vertices = vertices.dt, directed = F)
  hub.scores <- hub_score(net, weights = NA)$vector
  V(net)$Color <- v.colors.dt[V(net)$Phylum]$Color
  V(net)$Size <- (sqrt(hub.scores / max(hub.scores)) * 3) + 0.75
  E(net)$Color <- e.colors.dt[
    sapply(E(net)$Rho, function (rho) {
      which(e.colors.dt$Start < rho & e.colors.dt$End >= rho)
    })
  ]$Color
  set.seed(user.seed)
  l <- layout_with_kk(net)

  plot(
    net,
    vertex.label = NA,
    vertex.size = V(net)$Size,
    vertex.color = V(net)$Color,
    edge.color = E(net)$Color,
    edge.curved=.1,
    layout = l
  )
  title(paste(conc, "µM BaP"))
  if (conc == 0) {
    legend(
      title = "Phylum",
      x = 1.1,
      y = 1,
      v.colors.dt$Phylum,
      pch = 21,
      col = "#777777",
      pt.bg = v.colors.dt$Color,
      pt.cex = 2,
      cex = 0.8,
      bty = "n",
      ncol = 2,
      x.intersp = 0.5,
      title.adj = 0
    )
  }
  node.degree.all <- degree(net, mode = "all")
  network.stats.dt <- rbind(
    network.stats.dt,
    data.table(
      BaP_uM = conc,
      Edge.density = edge_density(net, loops = F),
      Mean.node.degree = mean(node.degree.all),
      Median.node.degree = median(node.degree.all),
      Max.node.degree = max(node.degree.all),
      Min.node.degree = min(node.degree.all),
      Mean.distance = mean_distance(net, directed = F)
    )
  )
}
saveRDS(
  network.stats.dt,
  file = file.path(saveDir, "dt_genus_spearman_network_statistics.rds")
)
saveRDS(
  genus.cors.sig,
  file = file.path(saveDir, "dt_genus_spearman_correlations_sig0_01.rds")
)
######

### spearman-top50pct-genus-plots, include=FALSE, eval=FALSE}
genus.taxonomy <- prune_taxa(top.genera, genus.ps) %>% taxa.data.table()
genus.taxonomy[, Taxon := Genus]
genus.taxonomy[, Genus := NULL]
set.redo.false("pick.random.colors")
v.colors.dt.file <- file.path(saveDir, "dt_genus_vertex_phylum_colors.rds")
v.colors.dt <- redo.if("pick.random.colors", v.colors.dt.file, {
  data.table(
    Phylum = sort(unique(genus.taxonomy$Phylum)),
    Color = randomColor(
      count = length(unique(genus.taxonomy$Phylum)),
      luminosity = "bright"
    )
  ) %>% setkeyv("Phylum")
})

top.genus.cors.sig <- genus.cors.sig[Taxon1 %in% top.genera & Taxon2 %in% top.genera]

set2.colors <- brewer.pal(8, "Set2")
n.bins <- 30
rho.bins <- bin(top.genus.cors.sig$Rho, nbins = n.bins)
e.colors.dt <- data.table(
  Rho.bin = levels(rho.bins),
  Color = colorRampPalette(c(set2.colors[2], "grey80", set2.colors[1]))(n.bins)
) %>% setkeyv("Rho.bin")
e.colors.dt$Start <- str_remove_all(e.colors.dt$Rho.bin, "\\(|\\]") %>%
  str_split(",") %>%
  sapply(`[`, 1) %>%
  as.numeric()
e.colors.dt$End <- str_remove_all(e.colors.dt$Rho.bin, "\\(|\\]") %>%
  str_split(",") %>%
  sapply(tail, 1) %>%
  as.numeric()
# shared.taxa <- sapply(
#   unique(c(top.genus.cors.sig$Taxon1, top.genus.cors.sig$Taxon2)),
#   function(taxon) {
#     pres.count <- sapply(bap.concs, function(conc) {
#       subset.dt <- top.genus.cors.sig[BaP_uM == conc]
#       taxon %in% unique(c(subset.dt$Taxon1, subset.dt$Taxon2))
#     }) %>% sum()
#     if (pres.count == 4) { return(taxon) }
# }) # this does not reduce the number of taxa plotted

network.stats.dt <- data.table()
par(mfrow = c(2,2), mar = c(1,1,1,1))
for (conc in bap.concs) {
  subset.dt <- top.genus.cors.sig[BaP_uM == conc & abs(Rho) > 0.3]
  edges.dt <- subset.dt[, .(Taxon1, Taxon2, Rho)]
  vertices.dt <- data.table(
    Taxon = unique(c(subset.dt$Taxon1, subset.dt$Taxon2))
  ) %>%
    merge(genus.taxonomy, by = "Taxon")
  net <- graph_from_data_frame(d = edges.dt, vertices = vertices.dt, directed = F)
  hub.scores <- hub_score(net, weights = NA)$vector
  V(net)$Color <- v.colors.dt[V(net)$Phylum]$Color
  V(net)$Size <- (sqrt(hub.scores / max(hub.scores)) * 2.5) + 1.5
  E(net)$Color <- e.colors.dt[
    sapply(E(net)$Rho, function (rho) {
      which(e.colors.dt$Start < rho & e.colors.dt$End >= rho)
    })
  ]$Color
  set.seed(user.seed)
  l <- layout_with_kk(net)

  plot(
    net,
    vertex.label = NA,
    vertex.size = V(net)$Size,
    vertex.color = V(net)$Color,
    edge.color = E(net)$Color,
    edge.curved = .1,
    layout = l
  )
  title(paste(conc, "µM BaP"))
  if (conc == 0) {
    legend(
      title = "Phylum",
      x = 1,
      y = 1,
      v.colors.dt$Phylum,
      pch = 21,
      col = "#777777",
      pt.bg = v.colors.dt$Color,
      pt.cex = 2,
      cex = 0.8,
      bty = "n",
      ncol = 2,
      x.intersp = 0.5,
      title.adj = 0
    )
  }
  node.degree.all <- degree(net, mode = "all")
  network.stats.dt <- rbind(
    network.stats.dt,
    data.table(
      BaP_uM = conc,
      Edge.density = edge_density(net, loops = F),
      Mean.node.degree = mean(node.degree.all),
      Median.node.degree = median(node.degree.all),
      Max.node.degree = max(node.degree.all),
      Min.node.degree = min(node.degree.all),
      Mean.distance = mean_distance(net, directed = F)
    )
  )
}
saveRDS(
  network.stats.dt,
  file = file.path(saveDir, "dt_genus_spearman_network_statistics.rds")
)
saveRDS(
  top.genus.cors.sig,
  file = file.path(saveDir, "dt_genus_spearman_correlations_sig0_01.rds")
)
#####

morph.dt <- read.csv("~/Downloads/Untitled spreadsheet - Sheet1.csv") %>% as.data.table()

morph.dt[, value := ifelse(is.na(value), 0, value)]
names(morph.dt)[str_which(names(morph.dt), "value")] <- "BaP.uM"
id.cols <- c("plate_id", "row_id", "col_id", "BaP.uM")
var.test <- names(morph.dt) %>%
  { str_detect(., "_$|[0-9]") | str_detect(., "^[A-Z]+$") }
var.cols <- names(morph.dt)[var.test]
morph.melt.dt <- melt(
  morph.dt,
  id.vars = id.cols,
  measure.vars = var.cols,
  variable.name = "Issue"
)
morph.melt.dt[, .(Count = sum(value, na.rm = T)), by = .(BaP.uM, Issue)] %>%
  dcast(BaP.uM ~ Issue)
