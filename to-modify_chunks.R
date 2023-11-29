## Alpha-diversity

```{r alpha-bapOnly-stats}
alpha.bap.res.file <- file.path(dirs$plots, "alphaDiv_bapOnly_stats_table.png")
flextable::flextable(alpha.bap.res) %>%
  set_caption(caption = "Alpha-diversity by BaP concentration") %>%
  theme_box() %>%
  align(j = 1:3, align = "right") %>%
  colformat_double(j = 5:7, digits = 2) %>%
  colformat_double(j = 8, digits = 3) %>%
  merge_v(j = c(1:3)) %>%
  autofit() %>%
  save_as_image(path = alpha.bap.res.file)
```

![](Plots/alphaDiv_bapOnly_stats_table.png)

### EPR


```{r alpha-epr-table}
to.print <- lapply(alpha.methods, function(alpha) {
  cat(paste0(alpha, ":"), sep = "\n")
  print(summary(alpha.epr.res[[alpha]]))
  cat('################', sep = "\n")
})
### MOVE THIS TO ALL ONE FLEXTABLE
```

```{r alpha-epr-plot, echo=FALSE, fig.height=10}
alpha.epr.plot.dt <- melt(
  alpha.epr.dt,
  measure.vars = alpha.methods,
  variable.name = "Alpha",
  value.name = "Score"
)
fig.cap <- "Plot of transformed area under the curve (AUC) values calculated from movement during the EPR test against alpha-diversity metric scores. EPR data was acquired 1 dpf, alpha-diversity metrics were assessed by sequencing gut contents 9 dpf."
{
  ggplot(
    alpha.epr.plot.dt[Alpha != "Observed"],
    aes(x = AUC.trans, y = Score, color = as.factor(BaP_uM))
  ) +
    geom_point() +
    stat_smooth(method = "lm", formula = y ~ x, se = F) +
    facet_wrap(~ Alpha, scales = "free_y") +
    bap.color.scale +
    labs(title = "EPR", x = "AUC (Tukey-transformed)", y = "Alpha Metric Score") +
    gg_figure_caption(fig.cap)
} %>%
  ggsave(
    filename = file.path(dirs$plots, "scatter_EPR_alphaDiv.png"),
    width = img.width,
    height = img.height,
    dpi = img.dpi
  )
```

![](Plots/scatter_EPR_alphaDiv.png)

### LPR

```{r alpha-lpr-tables, echo=FALSE, fig.height=4}
to.print <- lapply(alpha.methods, function(alpha) {
  lapply(lpr.cycles, function(cycle) {
    lapply(lpr.epochs, function(epoch) {
      mod <- alpha.lpr.res[[cycle]][[epoch]][[alpha]]
      cap <- paste0(
        alpha, ": ", toTitleCase(cycle),
        " cycle (", str_replace(epoch, "E", "Epoch "), ")"
      )
      cat(cap, sep = "\n")
      if (alpha.trans.dt[alpha]$Transformed) {
        summary(mod) %>%
          print()
      } else {
        summary(mod) %>%
          print()
      }
      cat.sep()
    })
  })
})

### MOVE THIS TO ALL ONE FLEXTABLE
```

```{r alpha-lpr-plot, echo=FALSE}
alpha.lpr.plot.dt <- melt(
  alpha.lpr.dt,
  measure.vars = alpha.methods,
  variable.name = "Alpha",
  value.name = "Score"
)
fig.cap <- "Plot of transformed area under the curve (AUC) values calculated from movement during the LPR test against alpha-diversity metric scores. LPR data was acquired 5 dpf, alpha-diversity metrics were assessed by sequencing gut contents at 9 dpf. LPR data was acquired 5 dpf and consists of three epochs of light and dark cycles."
{
  ggplot(
    alpha.lpr.plot.dt[Alpha != "Observed"],
    aes(x = AUC.trans, y = Score, color = as.factor(BaP_uM))
  ) +
    geom_point(alpha = 0.8) +
    stat_smooth(method = "lm", formula = y ~ x, se = F) +
    facet_grid(Alpha ~ Cycle + Epoch, scales = "free") +
    bap.color.scale +
    labs(
      title = "LPR",
      x = "AUC (Tukey-transformed)",
      y = "Alpha Metric Score"
    ) +
    gg_figure_caption(fig.cap, caption.width = cap.width * 3)
} %>%
  ggsave(
    filename = file.path(
      dirs$plots,
      "scatterPlot_lpr_alpha_by_AUC_BaP_Cycle_Epoch.png"
    ),
    width = img.width * 3,
    height = img.height * 2,
    dpi = img.dpi
  )
```

![figure](Plots/scatterPlot_lpr_alpha_by_AUC_BaP_Cycle_Epoch.png)

take the Simpson lines with a grain of salt, the actual models run with beta glms, not normal glms, and the link function was "logit", so those lines should look very different!

  ## Beta-diversity


  ### EPR

  ```{r beta-epr-stats}

beta.epr.res.file <- file.path(dirs$save, "list_beta_div_EPR_results.rds")
beta.epr.res <- redo.if("beta.stats", beta.epr.res.file, {
  lapply(beta.methods, function(beta) {
    cat(beta, sep = "\n")

    if (beta == "Aitchison") {
      beta.epr.df <- beta.epr.dfs.list$clr.trans
    } else {
      beta.epr.df <- beta.epr.dfs.list$rarefied
    }

    dist.mat <- dist.list[[beta]] %>% dist_subset(rownames(beta.epr.df))
    dbrda.full <- capscale(dist.mat ~ AUC.trans * BaP_uM, data = beta.epr.df)
    dbrda.step <- ordistep(dbrda.full, direction = "both")
    return(list(Full = dbrda.full, Step = dbrda.step))
  })
})

beta.epr.aovs.file <- file.path(dirs$save, "list_beta_div_epr_permanovas.rds")
beta.epr.aovs <- redo.if("beta.stats", beta.epr.aovs.file, {
  lapply(beta.methods, function(beta) {
    anova(beta.epr.res[[beta]]$Step, by = "term") %>%
      tidy() %>%
      as.data.table() %>%
      return()
  })
})
sig.epr.aov.id.dt <- lapply(beta.methods, function(beta) {
  aov <- beta.epr.aovs[[beta]]
  if (any(aov$p.value <= 0.05, na.rm = T)) {
    dt <- data.table(Method = beta)
  }
}) %>% rbindlist()
```


```{r beta-epr-tables, echo=FALSE, fig.height=4}
## Selected models

sig.epr.aov.mod.terms.dt <- lapply(1:nrow(sig.epr.aov.id.dt), function(i) {
  id <- sig.epr.aov.id.dt[i] %>% unlist()
  tbl <- beta.epr.aovs[[id]]
  cat(paste(id, collapse = " - "), sep = "\n")
  print(tbl)
  cat.sep()
  # flextable(tbl) %>%
  #   set_caption(caption = paste(ids, collapse = " - ")) %>%
  #   align(j = 1, align = "right") %>%
  #   colformat_double(j = 3:4, digits = 2) %>%
  #   colformat_double(j = 5, digits = 3) %>%
  #   autofit() %>%
  #   plot(zoom = 10)
  tbl[, names(id) := as.list(id)]
  return(tbl[term != "Residual", .(term, statistic, p.value, Method)])
}) %>% rbindlist(fill = T)

sig.epr.aov.mod.terms.dt[, Sig := p.value <= 0.05]
names(sig.epr.aov.mod.terms.dt)[1] <- "Variable"
```

```{r beta-epr-plot, echo=FALSE}
merge.cols <- c("Variable", "Method", "Sig")
biplot.data.lists.file <- file.path(dirs$save, "list_epr_beta_biplot_data.rds")
biplot.data.lists <- redo.if("beta.stats", biplot.data.lists.file, {
  lapply(1:nrow(sig.epr.aov.id.dt), function(i) {
    id <- sig.epr.aov.id.dt[i] %>% unlist()
    dbrda <- beta.epr.res[[id]]$Step
    biplot.data <- get.biplot.data(ps = ps.decontam, ord = dbrda)
    smpls.dt <- copy(biplot.data$sample.coords)
    smpls.dt[, names(id) := as.list(id)]
    names(smpls.dt)[2:3] <- c("Axis1", "Axis2")
    vctrs.dt0 <- biplot.data$vector.coords
    vctrs.dt0[, 2:3] <- vctrs.dt0[, 2:3] * biplot.data$coord.scale
    vctrs.dt0[, names(id) := as.list(id)]
    names(vctrs.dt0)[2:3] <- c("Axis1", "Axis2")
    vctrs.dt <- sig.epr.aov.mod.terms.dt[, ..merge.cols][
      vctrs.dt0,
      on = merge.cols[-3]
    ]
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
    axes.labs.dt[, names(id) := as.list(id)]
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
plot.list <- lapply(biplot.data.lists, function(x) {
  smpls.dt <- copy(x$Samples)
  smpls.dt <- smpls.dt[
    auc.dts$EPR[, .(Sample, AUC.trans)],
    on = "Origin==Sample",
    nomatch = 0
  ]
  vctrs.dt <- copy(x$Vectors)
  axes.labs.dt <- copy(x$Axes.labs)
  coord.scale <- x$Coord.scale
  plot <- ggplot(smpls.dt, aes(x = Axis1, y = Axis2))
  if (
    any(str_detect(vctrs.dt$Variable, ":")) |
    nrow(vctrs.dt) > 1 &
    all(vctrs.dt$Sig)
  ) {
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
        aes(
          xend = Axis1 / coord.scale,
          yend = Axis2 / coord.scale,
          color = Sig.lab
        ),
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
        geom_point(
          size = 2.5,
          shape = 21,
          color = "black",
          aes(fill = AUC.trans)
        ) +
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
      x = axes.labs.dt$Axis1,
      y = axes.labs.dt$Axis2,
      caption = wrap.caption(caption, w = 80)
    ) +
    theme(legend.position = "none")
  return(plot)
})

plot.list[[length(plot.list) + 1]] <- plot.env$full.legend
plot_grid(plotlist = plot.list, nrow = 2) %>%
  ggsave(
    filename = file.path(dirs$plots, "ordinations_epr.png"),
    width = img.width * 3,
    height = img.height * 2,
    dpi = img.dpi
  )
```

![](Plots/ordinations_epr.png)

### LPR

```{r beta-lpr-stats}


beta.lpr.res.file <- file.path(dirs$save, "list_beta_div_lpr_results.rds")
beta.lpr.res <- redo.if("beta.stats", beta.lpr.res.file, {
  lapply(beta.methods, function(beta) {
    lapply(lpr.cycles, function(cycle) {
      lapply(lpr.epochs, function(epoch) {
        if (beta == "Aitchison") {
          beta.lpr.df <- beta.lpr.dfs.list$clr.trans[[cycle]][[epoch]]
        } else {
          beta.lpr.df <- beta.lpr.dfs.list$rarefied[[cycle]][[epoch]]
        }
        dist.mat <- dist.list[[beta]] %>%
          dist_subset(rownames(beta.lpr.df))
        dbrda.full <- capscale(
          dist.mat ~ AUC.trans * BaP_uM,
          data = beta.lpr.df
        )
        dbrda.step <- ordistep(dbrda.full, direction = "both")
        return(list(Full = dbrda.full, Step = dbrda.step))
      })
    })
  })
})

beta.lpr.aovs.file <- file.path(dirs$save, "list_beta_div_lpr_permanovas.rds")
beta.lpr.aovs <- redo.if("beta.stats", beta.lpr.aovs.file, {
  lapply(beta.methods, function(beta) {
    lapply(lpr.cycles, function(cycle) {
      lapply(lpr.epochs, function(epoch) {
        anova(beta.lpr.res[[beta]][[cycle]][[epoch]]$Step, by = "term") %>%
          tidy() %>%
          as.data.table() %>%
          return()
      })
    })
  })
})
sig.lpr.aov.id.dt <- lapply(beta.methods, function(beta) {
  lapply(lpr.cycles, function(cycle) {
    lapply(lpr.epochs, function(epoch) {
      aov <- beta.lpr.aovs[[beta]][[cycle]][[epoch]]
      if (any(aov$p.value <= 0.05, na.rm = T)) {
        dt <- data.table(Method = beta, Cycle = cycle, Epoch = epoch)
      }
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()

```

```{r beta-lpr-stats-bapCat}
auc.dts$LPR1 <- copy(auc.dts$LPR)
auc.dts$LPR1[
  , BaP_cat := factor(
    paste0(BaP_uM, "uM"),
    levels = paste0(unique(auc.dts$LPR1$BaP_uM), "uM")
  )
]
beta.lpr.bapCat.dfs.list <- lapply(abundance.types, function(abund.type) {
  lapply(lpr.cycles, function(cycle) {
    lapply(lpr.epochs, function(epoch) {
      merge(
        sample.dts.list[[abund.type]][, .(Sequencing.ID, Sample)],
        auc.dts$LPR1[Cycle == cycle & Epoch == epoch],
        by = "Sample"
      ) %>%
        as.data.frame() %>%
        set_rownames(., .[["Sequencing.ID"]])
    })
  })
})

beta.lpr.bapCat.res.file <- file.path(
  dirs$save,
  "list_beta_div_lpr_bapCat_results.rds"
)
beta.lpr.bapCat.res <- redo.if("beta.stats", beta.lpr.bapCat.res.file, {
  lapply(beta.methods, function(beta) {
    lapply(lpr.cycles, function(cycle) {
      lapply(lpr.epochs, function(epoch) {
        beta.lpr.df <- beta.lpr.bapCat.dfs.list$rarefied
        if (beta == "Aitchison") {
          beta.lpr.df <- beta.lpr.bapCat.dfs.list$clr.trans
        }
        dist.mat <- dist.list[[beta]] %>%
          dist_subset(rownames(beta.lpr.df[[cycle]][[epoch]]))
        dbrda.full <- capscale(
          dist.mat ~ AUC.trans * BaP_cat,
          data = beta.lpr.df[[cycle]][[epoch]]
        )
        dbrda.step <- ordistep(dbrda.full, direction = "both")
        return(list(Full = dbrda.full, Step = dbrda.step))
      })
    })
  })
})

beta.lpr.bapCat.aovs.file <- file.path(
  dirs$save,
  "list_beta_div_lpr_bapCat_permanovas.rds"
)
beta.lpr.bapCat.aovs <- redo.if("beta.stats", beta.lpr.bapCat.aovs.file, {
  lapply(beta.methods, function(beta) {
    lapply(lpr.cycles, function(cycle) {
      lapply(lpr.epochs, function(epoch) {
        anova(beta.lpr.bapCat.res[[beta]][[cycle]][[epoch]]$Step, by = "term") %>%
          tidy() %>%
          as.data.table() %>%
          return()
      })
    })
  })
})
sig.lpr.aov.id.bapCat.dt <- lapply(beta.methods, function(beta) {
  lapply(lpr.cycles, function(cycle) {
    lapply(lpr.epochs, function(epoch) {
      aov <- beta.lpr.bapCat.aovs[[beta]][[cycle]][[epoch]]
      if (any(aov$p.value <= 0.05, na.rm = T)) {
        dt <- data.table(Method = beta, Cycle = cycle, Epoch = epoch)
      }
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()

```


```{r beta-lpr-tables, echo=FALSE, fig.height=4}
## Selected models

sig.lpr.aov.mod.terms.dt <- lapply(1:nrow(sig.lpr.aov.id.dt), function(i) {
  ids <- sig.lpr.aov.id.dt[i] %>% unlist()
  tbl <- beta.lpr.aovs[[ids[1]]][[ids[2]]][[ids[3]]]
  cat(paste(ids, collapse = " - "), sep = "\n")
  print(tbl)
  cat.sep()
  tbl[, names(ids) := as.list(ids)]
  return(
    tbl[term != "Residual", .(term, statistic, p.value, Method, Cycle, Epoch)]
  )
}) %>% rbindlist(fill = T)

sig.lpr.aov.mod.terms.dt[, Sig := p.value <= 0.05]
names(sig.lpr.aov.mod.terms.dt)[1] <- "Variable"
```

```{r beta-lpr-bapCat-tables, echo=FALSE, fig.height=4}
## Selected models

sig.lpr.aov.bapCat.mod.terms.dt <- lapply(1:nrow(sig.lpr.aov.id.bapCat.dt), function(i) {
  ids <- sig.lpr.aov.id.bapCat.dt[i] %>% unlist()
  tbl <- beta.lpr.bapCat.aovs[[ids[1]]][[ids[2]]][[ids[3]]]
  cat(paste(ids, collapse = " - "), sep = "\n")
  print(tbl)
  cat.sep()
  tbl[, names(ids) := as.list(ids)]
  return(
    tbl[term != "Residual", .(term, statistic, p.value, Method, Cycle, Epoch)]
  )
}) %>% rbindlist(fill = T)

sig.lpr.aov.bapCat.mod.terms.dt[, Sig := p.value <= 0.05]
names(sig.lpr.aov.bapCat.mod.terms.dt)[1] <- "Variable"
```

```{r beta-lpr-plot, echo=FALSE}
merge.cols <- c("Variable", "Method", "Cycle", "Epoch", "Sig")
biplot.data.lists.file <- file.path(dirs$save, "list_lpr_beta_biplot_data.rds")
biplot.data.lists <- redo.if("beta.stats", biplot.data.lists.file, {
  lapply(1:nrow(sig.lpr.aov.id.dt), function(i) {
    ids <- sig.lpr.aov.id.dt[i] %>% unlist()
    dbrda <- beta.lpr.res[[ids[1]]][[ids[2]]][[ids[3]]]$Step
    biplot.data <- get.biplot.data(ps = ps.decontam, ord = dbrda)
    smpls.dt <- copy(biplot.data$sample.coords)
    smpls.dt[, names(ids) := as.list(ids)]
    names(smpls.dt)[2:3] <- c("Axis1", "Axis2")
    vctrs.dt0 <- biplot.data$vector.coords
    vctrs.dt0[, 2:3] <- vctrs.dt0[, 2:3] * biplot.data$coord.scale
    vctrs.dt0[, names(ids) := as.list(ids)]
    names(vctrs.dt0)[2:3] <- c("Axis1", "Axis2")
    vctrs.dt <- sig.lpr.aov.mod.terms.dt[, ..merge.cols][
      vctrs.dt0,
      on = merge.cols[-5]
    ]
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
plot.list <- lapply(biplot.data.lists, function(x) {
  smpls.dt <- copy(x$Samples)
  # cat(smpls.dt$Method[1], sep = "\n")
  smpls.dt <- smpls.dt[
    auc.dts$LPR[, .(Sample, Cycle, Epoch, AUC.trans)],
    on = c("Origin==Sample", "Cycle", "Epoch"),
    nomatch = 0
  ]
  vctrs.dt <- copy(x$Vectors)
  axes.labs.dt <- copy(x$Axes.labs)
  coord.scale <- x$Coord.scale
  plot <- ggplot(smpls.dt, aes(x = Axis1, y = Axis2))
  if (
    any(str_detect(vctrs.dt$Variable, ":")) |
    nrow(vctrs.dt) > 1 &
    all(vctrs.dt$Sig)
  ) {
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
        aes(
          xend = Axis1 / coord.scale,
          yend = Axis2 / coord.scale,
          color = Sig.lab
        ),
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
        geom_point(
          size = 2.5,
          shape = 21,
          color = "black",
          aes(fill = AUC.trans)
        ) +
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
    filename = file.path(dirs$plots, "ordinations_lpr.png"),
    width = img.width * 5,
    height = img.height * 3,
    dpi = img.dpi
  )
```

![figure](Plots/ordinations_lpr.png)

## Random Forests


### BaP

```{r rfs-bap-rarVclr-asvVagg-tbls, echo=FALSE, fig.height=4}
# flextable(roc.aucs.dt[order(ROC.AUC, decreasing = T)]) %>%
#   plot(zoom = 10)
cat("Predicting BaP from taxon abundances (classification)", sep = "\n")
print(roc.aucs.dt[order(ROC.AUC, decreasing = T)])

# flextable(bap.rmses.dt[order(RMSE)])%>%
#   plot(zoom = 10)
cat("Predicting BaP from taxon abundance (regression)", sep = "\n")
print(bap.rmses.dt[order(RMSE)])
```

```{r rfs-bap-featureImportance}

```

### LPR AUC


```{r rfs-auc-rarVclr-asvVagg-tbl, echo=FALSE}
cat("Predicting LPR AUCs from taxon abundances (regression)", sep = "\n")
cat("All results:", sep = "\n")
lpr.rmses.dt[order(RMSE)]
cat("", sep = "\n")
cat("Best daate set per cycle and epoch:", sep = "\n")
lpr.rmses.dt[
  , .(
    Taxa.set = Taxa.set[which(RMSE == min(RMSE))],
    Abundance.type = Abundance.type[which(RMSE == min(RMSE))],
    RMSE = min(RMSE)
  ),
  by = c("LPR.cycle", "LPR.epoch")
]
```

```{r rfs-auc-featureImportance}

```

## LMs

### Foward-building models


#### Significant __interactions__

```{r lpr-lms-sigIntrxn-table, echo=FALSE}
get.taxa.labs <- function(plot.taxa) {
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
  return(taxa.labs.dt)
}

epoch.colors <- c("cadetblue1", "cadetblue3", "cadetblue4") %>%
  set_names(lpr.epochs)
print.dt0 <- copy(intrxn.reg.data)
taxa.labs.dt <- get.taxa.labs(unique(print.dt0$Taxon))
print.dt <- merge(
  print.dt0,
  taxa.labs.dt[, .(Orig.ID, Lab)],
  by.x = "Taxon",
  by.y = "Orig.ID"
)
setorderv(
  print.dt,
  cols = c("Epoch", "Cycle", "Importance", "Lab"),
  order = c(1, 1, -1, 1)
)
table.file <- file.path(dirs$plots, "table_LPR_LMs_sigInteractions.html")
flextable(print.dt[, c(6:5, 9, 2, 4, 3)]) %>%
  colformat_double(j = 5:6, digits = 3) %>%
  merge_v(j = 1:3) %>%
  bg(j = 1, bg = function(x) epoch.colors[x]) %>%
  bg(j = 2, bg = function(x) ifelse(x == "light", "white", "black")) %>%
  color(j = 2, color = function(x) ifelse(x == "light", "black", "white")) %>%
  bg(j = 4, bg = bap.colors) %>%
  color(j = 5, i = ~ Slope < 0, color = "red") %>%
  autofit() %>%
  theme_box() %>%
  save_as_html(path = table.file)
htmltools::includeHTML(table.file)
```


##### LPR dark cycle (epoch 1) only

```{r lpr-lms-sigIntrxn-darkCycle-E1-plots, include=FALSE}
dark.e1.reg.dt <- copy(intrxn.reg.data)[Cycle == "dark" & Epoch %in% "E1"]
impt.dt <- dark.e1.reg.dt[
  , .(Importance = Importance[1]),
  by = .(Taxon)
]
max.facets <- 16
for (i in seq(1, nrow(impt.dt), by = max.facets)) {
  cat(i, sep = "\n")
  end <- min(c(i + max.facets - 1, nrow(impt.dt)))
  plot.taxa <- impt.dt[order(Importance, decreasing = T)]$Taxon[i:end]
  taxa.labs.dt <- get.taxa.labs(plot.taxa)

  plot.taxa.dt <- merge(
    dark.e1.reg.dt[Taxon %in% plot.taxa],
    taxa.labs.dt[, .(Orig.ID, Lab)],
    by.x = "Taxon",
    by.y = "Orig.ID"
  )

  keep.cols <- c("BaP_uM", "AUC.trans", plot.taxa)
  lpr.taxa.set <- best.lpr.rf.datasets[LPR.cycle == "dark" & LPR.epoch == "E1"] %>%
    use_series(Taxa.set)
  lpr.abund.set <- best.lpr.rf.datasets[
    LPR.cycle == "dark" & LPR.epoch == "E1"
  ] %>%
    use_series(Abundance.type)
  plot.data <- rbind(
    copy(lpr.rf.data.list[[lpr.taxa.set]][[lpr.abund.set]]$dark$E1$Train),
    copy(lpr.rf.data.list[[lpr.taxa.set]][[lpr.abund.set]]$dark$E1$Test)
  )[, ..keep.cols] %>%
    melt(
      id.vars = c("BaP_uM", "AUC.trans"),
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
  p <- ggplot(plot.data, aes(x = AUC.trans, y = sqrt(Abund))) +
    geom_point(aes(color = as.factor(BaP_uM))) +
    geom_abline(
      data = plot.taxa.dt,
      aes(intercept = Intercept, slope = Slope, color = as.factor(BaP_uM))
    ) +
    facet_wrap(~ Lab, scales = "free_y", ncol = 4) +
    bap.color.scale +
    labs(
      caption = ,
      x = "AUC (Tukey-transformed)",
      y = "sqrt(Abundance)"
    )
  plot.file <- file.path(
    dirs$plots,
    paste0(
      "scatterPlots_darkCycle_E1_sigImptTaxa", i, "-", end,"_wRegressions.png"
    )
  )
  ggsave(
    p,
    filename = plot.file,
    width = img.width * 3,
    height = img.height * 2,
    dpi = img.dpi
  )
}
```

![](Plots/scatterPlots_darkCycle_E1_sigImptTaxa1-16_wRegressions.png)

![](Plots/scatterPlots_darkCycle_E1_sigImptTaxa17-32_wRegressions.png)

![](Plots/scatterPlots_darkCycle_E1_sigImptTaxa33-48_wRegressions.png)

![](Plots/scatterPlots_darkCycle_E1_sigImptTaxa49-64_wRegressions.png)

![](Plots/scatterPlots_darkCycle_E1_sigImptTaxa65-80_wRegressions.png)

![](Plots/scatterPlots_darkCycle_E1_sigImptTaxa81-96_wRegressions.png)

![](Plots/scatterPlots_darkCycle_E1_sigImptTaxa97-112_wRegressions.png)

![](Plots/scatterPlots_darkCycle_E1_sigImptTaxa113-127_wRegressions.png)

#### Significant __main effects__ only

The plots below show significant effects of BaP (changes in y-intercept) or AUC (slopes different from zero), or both. I have not yet indicated which effects are significant on the plots themselves, but some should be quite clear (for example differing y-intercepts but nearly flat slopes would indicate significant differences in BaP exposure only).

```{r lpr-lms-sigMainEffects-all-plots, include=FALSE}
impt.dt <- copy(main.effects.reg.data)[, .(Importance = Importance[1]), by = Taxon]
for (i in seq(1, nrow(impt.dt), by = 12)) {
  cat(i, sep = "\n")
  end <- min(c(i + 11, nrow(impt.dt)))
  plot.taxa <- impt.dt[order(Importance, decreasing = T)]$Taxon[i:end] %>%
    str_replace_all("env\\.OPS", "envOPS")
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

  plot.taxa.dt <- merge(
    main.effects.reg.data[Taxon %in% plot.taxa],
    taxa.labs.dt[, .(Orig.ID, Lab)],
    by.x = "Taxon",
    by.y = "Orig.ID"
  )

  data.info.dt <- unique(sig.main.effects.dt[, .(Taxon, Cycle, Epoch)])
  lpr.taxa.sets <- sapply(1:nrow(data.info.dt), function(j) {
    row <- unique(sig.main.effects.dt[, .(Taxon, Cycle, Epoch)])[j]
    best.lpr.rf.datasets[LPR.cycle == row$Cycle & LPR.epoch == row$Epoch]$Taxa.set
  })
  lpr.abund.sets <-sapply(1:nrow(data.info.dt), function(j) {
    row <- unique(sig.main.effects.dt[, .(Taxon, Cycle, Epoch)])[j]
    best.lpr.rf.datasets[
      LPR.cycle == row$Cycle & LPR.epoch == row$Epoch
    ]$Abundance.type
  })
  plot.data <- lapply(1:nrow(data.info.dt), function(j) {
    # cat(j, sep = "\n")
    keep.cols <- c("BaP_uM", "AUC.trans", data.info.dt[j]$Taxon)
    taxon.dt <- rbind(
      lpr.rf.data.list %>%
        extract2(lpr.taxa.sets[j]) %>%
        extract2(lpr.abund.sets[j]) %>%
        extract2(data.info.dt[j]$Cycle) %>%
        extract2(data.info.dt[j]$Epoch) %>%
        extract2("Train") %>%
        .[, ..keep.cols],
      lpr.rf.data.list %>%
        extract2(lpr.taxa.sets[j]) %>%
        extract2(lpr.abund.sets[j]) %>%
        extract2(data.info.dt[j]$Cycle) %>%
        extract2(data.info.dt[j]$Epoch) %>%
        extract2("Train") %>%
        .[, ..keep.cols]
    )
    names(taxon.dt)[3] <- "Abund"
    taxon.dt[, Taxon := data.info.dt[j]$Taxon]
    taxon.dt[data.info.dt[j], on = "Taxon"] %>%
      merge(taxa.labs.dt[, .(Orig.ID, Lab)],
            by.x = "Taxon",
            by.y = "Orig.ID",
            all.y = FALSE
      ) %>%
      return()
  }) %>% rbindlist()

  p <- ggplot(plot.data, aes(x = AUC.trans, y = sqrt(Abund))) +
    geom_point(aes(color = as.factor(BaP_uM))) +
    geom_abline(
      data = plot.taxa.dt,
      aes(intercept = Intercept, slope = Slope, color = as.factor(BaP_uM))
    ) +
    bap.color.scale +
    facet_wrap(~ Lab + Cycle + Epoch, scales = "free", nrow = 3) +
    panel_border() +
    labs(
      caption = paste(
        "All taxa that had significant relationships",
        "with an interaction between LPR cycles AUCs"
      ),
      x = "AUC (Tukey-transformed)",
      y = "sqrt(Abundance)"
    )
  plot.file <- file.path(
    dirs$plots,
    paste0("scatterPlots_mainEffects_sigImptTaxa", i, "-", end,"_wRegressions.png")
  )
  ggsave(
    p,
    filename = plot.file,
    width = img.width * 3,
    height = img.width * 2,
    dpi = img.dpi
  )
}
```

![](Plots/scatterPlots_mainEffects_sigImptTaxa1-12_wRegressions.png)

![](Plots/scatterPlots_mainEffects_sigImptTaxa13-24_wRegressions.png)

![](Plots/scatterPlots_mainEffects_sigImptTaxa25-36_wRegressions.png)

![](Plots/scatterPlots_mainEffects_sigImptTaxa37-42_wRegressions.png)

```{r lpr-lms-sigMainEffects-for-Tom, include=FALSE, eval=FALSE}
focal.taxa <- c("ASV00028")
focal.cycles <- c("light")
focal.epochs <- c("E1")
impt.dt <- copy(main.effects.reg.data)[
  Cycle %in% focal.cycles & Epoch %in% focal.epochs & Taxon %in% focal.taxa,
  .(Importance = Importance[1]),
  by = Taxon
]
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

plot.taxa.dt <- merge(
  main.effects.reg.data[
    Cycle %in% focal.cycles & Epoch %in% focal.epochs & Taxon %in% focal.taxa
  ],
  taxa.labs.dt[, .(Orig.ID, Lab)],
  by.x = "Taxon",
  by.y = "Orig.ID"
)

data.info.dt <- unique(sig.main.effects.dt[, .(Taxon, Cycle, Epoch)])
lpr.taxa.sets <- sapply(1:nrow(data.info.dt), function(j) {
  row <- unique(sig.main.effects.dt[, .(Taxon, Cycle, Epoch)])[j]
  best.lpr.rf.datasets[LPR.cycle == row$Cycle & LPR.epoch == row$Epoch]$Taxa.set
})
lpr.abund.sets <-sapply(1:nrow(data.info.dt), function(j) {
  row <- unique(sig.main.effects.dt[, .(Taxon, Cycle, Epoch)])[j]
  best.lpr.rf.datasets[LPR.cycle == row$Cycle & LPR.epoch == row$Epoch]$Abundance.type
})
plot.data <- lapply(1:nrow(data.info.dt), function(j) {
  # cat(j, sep = "\n")
  keep.cols <- c("BaP_uM", "AUC.trans", data.info.dt[j]$Taxon)
  taxon.dt <- rbind(
    lpr.rf.data.list %>%
      extract2(lpr.taxa.sets[j]) %>%
      extract2(lpr.abund.sets[j]) %>%
      extract2(data.info.dt[j]$Cycle) %>%
      extract2(data.info.dt[j]$Epoch) %>%
      extract2("Train") %>%
      .[, ..keep.cols],
    lpr.rf.data.list %>%
      extract2(lpr.taxa.sets[j]) %>%
      extract2(lpr.abund.sets[j]) %>%
      extract2(data.info.dt[j]$Cycle) %>%
      extract2(data.info.dt[j]$Epoch) %>%
      extract2("Train") %>%
      .[, ..keep.cols]
  )
  names(taxon.dt)[3] <- "Abund"
  taxon.dt[, Taxon := data.info.dt[j]$Taxon]
  taxon.dt[data.info.dt[j], on = "Taxon"] %>%
    merge(taxa.labs.dt[, .(Orig.ID, Lab)],
          by.x = "Taxon",
          by.y = "Orig.ID",
          all.y = FALSE
    ) %>%
    return()
}) %>%
  rbindlist() %>%
  .[Cycle %in% focal.cycles & Epoch %in% focal.epochs & Taxon %in% focal.taxa]

plot.file <- file.path(
  tomDir,
  paste0(
    "scatter_and_beeswarm_",
    paste(focal.taxa, collapse = "-"),
    "_lightCycle_E1.pdf"
  )
)
plot.list <- NULL
plot.list[[1]] <- ggplot(plot.data, aes(x = AUC, y = sqrt(Abund))) +
  geom_point(aes(color = as.factor(BaP_uM))) +
  geom_abline(
    data = plot.taxa.dt,
    aes(intercept = Intercept, slope = Slope, color = as.factor(BaP_uM))
  ) +
  bap.color.scale +
  facet_wrap(~ Lab, scales = "free", nrow = 3) +
  panel_border() +
  labs(
    title = paste(tools::toTitleCase(focal.cycles), "Cycle"),
    subtitle = focal.epochs,
    x = "AUC (µm × s)",
    y = "sqrt(Abundance)"
  )


plot.list[[2]] <- ggplot(plot.data, aes(x = BaP_uM, y = sqrt(Abund))) +
  geom_quasirandom(aes(color = as.factor(BaP_uM))) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", color = "black") +
  bap.color.scale +
  facet_wrap(~ Lab, scales = "free_y", ncol = 4) +
  labs(
    x = "BaP exposure (µM)",
    y = "sqrt(Abundance)"
  ) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank()
  )


plot_grid(
  plotlist = plot.list,
  nrow = 1,
  align = "hv",
  axis = "lt",
  rel_widths = c(2, 1)
) %>%
  ggsave(
    filename = plot.file,
    width = 16,
    height = 8
  )
```

## SpiecEasi

```{r speiceasi-genera}

```

```{r speiceasi-genera-plot-v1, include=FALSE}
set.redo.false("pick.random.colors")
v.colors.dt.file <- file.path(dirs$save, "dt_genus_vertex_phylum_colors.rds")
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

network.stats.dt <- data.table()
CairoPNG(
  file = file.path(dirs$plots, "network_speicEasi_genera_byBaP.png"),
  width = img.width * 7.5,
  height = img.width * 4,
  dpi = img.dpi,
  units = "in"
)
par(mfrow = c(2,2), mar = c(2,2,2,2))
for (conc in bap.concs) {
  ps <- prune_samples(
    sample.data.table(ps.genera)[BaP_uM == conc]$Sample,
    ps.genera
  ) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    return()
  genus.taxonomy <- taxa.data.table(ps)
  genus.taxonomy[, Taxon := Genus]
  genus.taxonomy[, Genus := NULL]
  setkey(genus.taxonomy, Taxon)
  ig.se <- adj2igraph(getRefit(se.bap.genera[[paste0("uM", conc)]]))
  vsize <- log10(colMeans(otu.matrix(ps)) + 1) + 2
  vcolor <- v.colors.dt[genus.taxonomy[taxa_names(ps)]$Phylum]$Color
  se.coord <- layout.kamada.kawai(ig.se)

  se.beta <- getOptBeta(se.bap.genera[[paste0("uM", conc)]]) %>%
    symBeta(mode='maxabs') %>%
    summary()
  beta.bins <- bin(sort(se.beta$x), nbins = n.bins)
  e.colors.dt <- data.table(
    Beta.bin = levels(beta.bins)
  )
  e.colors.dt$Start <- str_remove_all(e.colors.dt$Beta.bin, "\\(|\\]") %>%
    str_split(",") %>%
    sapply(`[`, 1) %>%
    as.numeric()
  e.colors.dt$End <- str_remove_all(e.colors.dt$Beta.bin, "\\(|\\]") %>%
    str_split(",") %>%
    sapply(tail, 1) %>%
    as.numeric()
  e.colors.dt$Color <- c(
    colorRampPalette(c(set2.colors[2], "grey80"))(sum(e.colors.dt$Start < 0)),
    colorRampPalette(c("grey80", set2.colors[1]))(sum(e.colors.dt$Start >= 0))
  )
  ecolor <- e.colors.dt[
    sapply(se.beta$x, function (beta) {
      which(e.colors.dt$Start < beta & e.colors.dt$End >= beta)
    })
  ]$Color

  plot(
    ig.se,
    layout = se.coord,
    vertex.size = vsize,
    vertex.label = NA,
    vertex.color = vcolor,
    edge.color = ecolor,
    edge.width = abs(se.beta$x) * 10
  )
  title(
    main = paste(conc, "µM B[a]P"),
    cex.main = 3
  )
  if (conc == 0) {
    legend(
      title = "Phylum",
      x = 1,
      y = 1,
      v.colors.dt$Phylum,
      pch = 21,
      col = "#777777",
      pt.bg = v.colors.dt$Color,
      pt.cex = 3,
      cex = 2.5,
      bty = "n",
      ncol = 2,
      x.intersp = 1,
      title.adj = 0
    )
  }
  node.degree.all <- degree(ig.se, mode = "all")
  network.stats.dt <- rbind(
    network.stats.dt,
    data.table(
      BaP_uM = conc,
      Edge.density = round(edge_density(ig.se, loops = F), 4),
      Mean.node.degree = round(mean(node.degree.all), 4),
      Median.node.degree = median(node.degree.all),
      Max.node.degree = max(node.degree.all),
      Min.node.degree = min(node.degree.all),
      Mean.distance = round(mean_distance(ig.se, directed = F), 4)
    )
  )
}
dev.off()
saveRDS(
  network.stats.dt,
  file = file.path(dirs$save, "dt_genus_spiecEasi_network_statistics.rds")
)
```

![](Plots/network_speicEasi_genera_byBaP.png)

```{r speiceasi-genera-table, echo=FALSE}
network.stats.dt
```


```{r speiceasi-genera-plot-v2, include=FALSE, eval=FALSE}
for (conc in bap.concs) {
  ps <- prune_samples(sample.data.table(ps.genera)[BaP_uM == conc]$Sample, ps.genera) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    return()
  ig.se <- adj2igraph(getRefit(se.bap.genera[[paste0("uM", conc)]]))
}
```

# Spearman correlations

```{r spearman-ASV-cors, include=FALSE, eval=FALSE}
taxon.pairs.file <- file.path(dirs$save, "list_all_taxon_pairs_by_abundType.rds")
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
spearman.cors.file <- file.path(
  dirs$save,
  "list_spearman_cors_by_abundType_and_bap.rds"
)
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
```

```{r spearman-genus-cors, include=FALSE, eval=FALSE}
genus.ps.file <- file.path(dirs$save, "physeq_genus_glom.rds")
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

genus.pairs.file <- file.path(dirs$save, "dt_all_genus_pairs_by_abundType.rds")
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


genus.cors.file <- file.path(dirs$save, "dt_genus_spearman_cors_by_bap.rds")
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
  ggsave(
    filename = file.path(
      tomDir, "histogram_genus_spearman_rhos_logCounts.pdf"
    )
  )
```

```{r spearman-genus-plots, include=FALSE, eval=FALSE}
genus.taxonomy <- taxa.data.table(genus.ps)
genus.taxonomy[, Taxon := Genus]
genus.taxonomy[, Genus := NULL]
set.redo.false("pick.random.colors")
v.colors.dt.file <- file.path(dirs$save, "dt_genus_vertex_phylum_colors.rds")
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
  file = file.path(dirs$save, "dt_genus_spearman_network_statistics.rds")
)
saveRDS(
  genus.cors.sig,
  file = file.path(dirs$save, "dt_genus_spearman_correlations_sig0_01.rds")
)
```

```{r spearman-top50pct-genus-plots, include=FALSE, eval=FALSE}
genus.taxonomy <- prune_taxa(top.genera, genus.ps) %>% taxa.data.table()
genus.taxonomy[, Taxon := Genus]
genus.taxonomy[, Genus := NULL]
set.redo.false("pick.random.colors")
v.colors.dt.file <- file.path(dirs$save, "dt_genus_vertex_phylum_colors.rds")
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
    edge.curved=.1,
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
  file = file.path(dirs$save, "dt_genus_spearman_network_statistics.rds")
)
saveRDS(
  top.genus.cors.sig,
  file = file.path(dirs$save, "dt_genus_spearman_correlations_sig0_01.rds")
)
```

