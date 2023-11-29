# functions.R
dbrda.permanova <- function(dbrda, beta, cycle = NULL) {
  aov <- anova(dbrda, by = "term") %>%
    helpers::anova2dt() %>%
    .[, Beta.metric := beta]
  if (!is.null(cycle)) {
    aov[, Cycle := cycle]
    neworder <- {
      ncol(aov) - 1
    }:ncol(aov)
  } else {
    neworder <- ncol(aov)
  }
  setcolorder(aov, neworder)
  return(aov)
}

reformat.data <- function(inFile) {
  dt0 <- read_csv(inFile) %>%
    as.data.table() %>%
    set_colnames(gsub("id", "ID", tools::toTitleCase(names(.)))) %>%
    subset(complete.cases(.))
  if ("Remove" %in% names(dt0)) {
    dt0[, Remove := NULL]
  }
  names(dt0)[which(names(dt0) == "Conc")] <- "BaP.uM"
  dt0[, `:=`(Chemical.ID = NULL, Bottle.ID = NULL)] # remove unnecessary columns
  idVars <- str_subset(names(dt0), "^T", negate = T)
  dt <- melt(
    dt0,
    id.vars = idVars,
    variable.name = "Timepoint",
    value.name = "Movement"
  ) # melting data into a usable format for analyses and plotting
  dt[, Timepoint := as.numeric(str_remove(Timepoint, "T"))] # remove the "T" from the timpoints and make numeric for plotting and analyses
  return(dt)
}

prep.behavior.data <- function(
    physeq,
    behavior.assays = c("EPR", "LPR"),
    behavior.files,
    metadata.file,
    output.dir
) {
  behavior.assays <- toupper(behavior.assays)
  if (!(any(assays %in% c("EPR", "LPR")))) {
    stop("Argument `assay` must be either EPR or LPR")
  }

  for (assay in behavior.assays) {
    sample.dt <- sample.data.table(physeq)
    base.dt <- reformat.data(behavior.files[assay])
    if (assay == "EPR") {
      saveRDS(base.dt, file = file.path(output.dir, paste0(assay, "_base_dt.rds")))
      auc.dt <- sample.dt[
        base.dt,
        on = c(SARL_plateID = "Plate.ID", Sample_well = "Well")
      ] %>%
        subset(Timepoint > 29 & Timepoint < 41 & Dissect_result != "no gut") %>%
        .[
          , .(EPR.auc = pracma::trapz(x = Timepoint, y = Movement)),
          by = names(sample.dt)
        ]
    } else {
      saveRDS(base.dt, file = file.path(output.dir, paste0(assay, "_base_dt.rds")))
      epoch.dt <- data.table(
        Epoch = c(1, 2, 3),
        Start = c(61, 120, 180),
        End = c(89, 149, 209)
      )
      base.dt[, Sec := floor((Timepoint - 61) / 5) * 30]
      cycle.epoch.dt <- sample.dt[
        base.dt,
        on = c(SARL_plateID = "Plate.ID", Sample_well = "Well")
      ] %>%
        subset(Timepoint >= epoch.dt$Start[1] & Dissect_result != "no gut") %>%
        .[
          , Epoch := ifelse(
            Timepoint < epoch.dt[[2]][1],
            "E0",
            ifelse(
              between(Timepoint, epoch.dt[[2]][1], epoch.dt[[2]][2]),
              "E1",
              ifelse(
                between(Timepoint, epoch.dt[[2]][2], epoch.dt[[2]][3]),
                "E2",
                "E3"
              )
            )
          )
        ] %>%
        subset(Epoch != "E0") %>%
        .[
          , Cycle := ifelse(
            between(Timepoint, epoch.dt[[2]][1], epoch.dt[[3]][1]) |
              between(Timepoint, epoch.dt[[2]][2], epoch.dt[[3]][2]) |
              between(Timepoint, epoch.dt[[2]][3], epoch.dt[[3]][3]),
            "light",
            "dark"
          )
        ]

      quant.dt <- lapply(unique(cycle.epoch.dt$SARL_plateID), function(plate) {
        lapply(c(0, 1, 5, 10), function(conc) {
          lapply(unique(cycle.epoch.dt$Timepoint), function(time) {
            sub.dt <- cycle.epoch.dt[
              SARL_plateID == plate & BaP_uM == conc & Timepoint == time
            ]
            quants <- quantile(sub.dt$Movement, c(0.05, 0.95))
            sub.dt[, Movement := ifelse(Movement %between% quants, Movement, NA)]
            return(sub.dt)
          }) %>% rbindlist()
        }) %>% rbindlist()
      }) %>% rbindlist()

      auc.dt <- NULL
      for (epoch in paste0("E", 1:3)) {
        for (cycle in c("light", "dark")) {
          auc.col.name <- paste("LPR", epoch, cycle, "auc", sep = ".")
          sub.auc.dt <- quant.dt[Epoch == epoch & Cycle == cycle] %>%
            .[
              complete.cases(.[, .(Timepoint, Movement)]),
              .(New.auc = pracma::trapz(Timepoint, Movement)),
              by = names(sample.dt)
            ]

          if (is.null(auc.dt)) {
            auc.dt <- sub.auc.dt
            names(auc.dt)[which(names(auc.dt) == "New.auc")] <- auc.col.name
          } else {
            auc.dt[, (auc.col.name) := sub.auc.dt$New.auc]
          }
        }
      }
      light.cols <- str_subset(names(auc.dt), "LPR.E[1-3].light.auc")
      dark.cols <- str_subset(names(auc.dt), "LPR.E[1-3].dark.auc")
      auc.dt[, LPR.light.auc := rowSums(.SD), .SDcols = light.cols]
      auc.dt[, LPR.dark.auc := rowSums(.SD), .SDcols = dark.cols]
    }
    missing.smpls <- sample.dt$Sample[!(sample.dt$Sample %in% auc.dt$Sample)]
    auc.dt <- rbind(auc.dt, sample.dt[Sample %in% missing.smpls], fill = TRUE)
    auc.df <- as.data.frame(auc.dt) %>%
      magrittr::set_rownames(auc.dt$Sample)
    auc.df$Sample <- NULL
    sample_data(physeq) <- sample_data(auc.df)
  }
  return(physeq)
}

prep.behavior.data.noPS <- function(
    behavior.assays = c("EPR", "LPR"),
    behavior.files,
    metadata.file,
    output.dir
) {
  behavior.assays <- toupper(behavior.assays)
  if (!(any(assays %in% c("EPR", "LPR")))) {
    rlang::abort("Argument `assay` must be either EPR or LPR")
  }
  sample.dt <- read.csv(metadata.file) %>% as.data.table()
  by.cols <- c("Well", "Plate_ID", "BaP.uM")
  epr.dt <- NULL
  lpr.dt <- NULL
  for (assay in behavior.assays) {
    base.dt <- reformat.data(behavior.files[assay])

    if (assay == "EPR") {
      saveRDS(base.dt, file = file.path(
        output.dir, paste0("cv_cvz_gf_", assay, "_base_dt.rds")
      ))
      epr.dt <- merge(
        sample.dt,
        base.dt,
        all.y = TRUE,
        by.x = "Plate_ID",
        by.y = "Plate.ID"
      ) %>%
        subset(Timepoint > 29 & Timepoint < 41) %>%
        .[
          , .(EPR.auc = pracma::trapz(x = Timepoint, y = Movement)),
          by = by.cols
        ]
      final.by.cols <- names(epr.dt)
    } else {
      saveRDS(
        base.dt,
        file = file.path(output.dir, paste0("cv_cvz_gf_", assay, "_base_dt.rds"))
      )
      epoch.dt <- data.table(
        Epoch = c(1, 2, 3),
        Start = c(61, 120, 180),
        End = c(89, 149, 209)
      )
      base.dt[, Sec := floor((Timepoint - 61) / 5) * 30]
      cycle.epoch.dt <- merge(
        sample.dt,
        base.dt,
        all.y = TRUE,
        by.x = "Plate_ID",
        by.y = "Plate.ID"
      ) %>%
        subset(Timepoint >= epoch.dt$Start[1]) %>%
        .[
          , Epoch := ifelse(
            Timepoint < epoch.dt[[2]][1],
            "E0",
            ifelse(
              between(Timepoint, epoch.dt[[2]][1], epoch.dt[[2]][2]),
              "E1",
              ifelse(
                between(Timepoint, epoch.dt[[2]][2], epoch.dt[[2]][3]),
                "E2",
                "E3"
              )
            )
          )
        ] %>%
        subset(Epoch != "E0") %>%
        .[
          , Cycle := ifelse(
            between(Timepoint, epoch.dt[[2]][1], epoch.dt[[3]][1]) |
              between(Timepoint, epoch.dt[[2]][2], epoch.dt[[3]][2]) |
              between(Timepoint, epoch.dt[[2]][3], epoch.dt[[3]][3]),
            "light",
            "dark"
          )
        ]

      quant.dt <- lapply(unique(cycle.epoch.dt$Plate_ID), function(plate) {
        lapply(c(0, 1, 5, 10), function(conc) {
          lapply(unique(cycle.epoch.dt$Timepoint), function(time) {
            sub.dt <- cycle.epoch.dt[
              Plate_ID == plate & BaP.uM == conc & Timepoint == time
            ]
            quants <- quantile(sub.dt$Movement, c(0.05, 0.95))
            sub.dt[, Movement := ifelse(Movement %between% quants, Movement, NA)]
            return(sub.dt)
          }) %>% rbindlist()
        }) %>% rbindlist()
      }) %>% rbindlist()

      lpr.dt <- NULL
      for (epoch in paste0("E", 1:3)) {
        for (cycle in c("light", "dark")) {
          auc.col.name <- paste("LPR", epoch, cycle, "auc", sep = ".")
          sub.lpr.dt <- quant.dt[Epoch == epoch & Cycle == cycle] %>%
            .[
              complete.cases(.[, .(Timepoint, Movement)]),
              .(New.auc = pracma::trapz(Timepoint, Movement)),
              by = by.cols
            ]

          if (is.null(lpr.dt)) {
            lpr.dt <- sub.lpr.dt
            names(lpr.dt)[which(names(lpr.dt) == "New.auc")] <- auc.col.name
          } else {
            lpr.dt <- merge(
              lpr.dt,
              sub.lpr.dt[, .(Well, Plate_ID, New.auc)],
              by = c("Well", "Plate_ID"),
              all = TRUE
            )
            names(lpr.dt)[which(names(lpr.dt) == "New.auc")] <- auc.col.name
          }
        }
      }
      light.cols <- str_subset(names(lpr.dt), "LPR.E[1-3].light.auc")
      dark.cols <- str_subset(names(lpr.dt), "LPR.E[1-3].dark.auc")
      lpr.dt[, LPR.light.auc := rowSums(.SD), .SDcols = light.cols]
      lpr.dt[, LPR.dark.auc := rowSums(.SD), .SDcols = dark.cols]
    }
  }
  auc.dt <- merge(epr.dt, lpr.dt, all = T, by = by.cols)
  auc.dt[, Sample := paste(Well, Plate_ID, sep = "_")]
  names(auc.dt)[which(names(auc.dt) == "BaP.uM")] <- "BaP_uM"
  setkey(auc.dt, Sample)
  setcolorder(auc.dt, ncol(auc.dt))
  return(auc.dt)
}

prep.phyloseq <- function(
    input.ps.file,
    qubit.data.file,
    output.path
) {
  input.ps <- readRDS(input.ps.file)
  input.stem <- str_split(input.ps.file, "/")[[1]] %>%
    tail(1) %>%
    str_remove(".rds")
  sample0.dt <- sample.data.table(input.ps)
  qubit.dt <- readxl::read_excel(qubit.data.file) %>%
    as.data.table()

  # include qubit results for `decontam` processing
  sample1.df <- sample0.dt[qubit.dt, on = .(Barcode, Run), nomatch = 0] %>%
    as.data.frame()
  row.names(sample1.df) <- sample1.df$Sample
  sample1.df[, c("Sample", "Sample.1")] <- NULL
  sample_data(input.ps) <- sample1.df

  # De-dupe samples that got sequenced more than once (due to some PCR issues)
  keep.samples <- sapply(unique(sample0.dt$Origin), function(origin) {
    sample_sums(input.ps)[sample0.dt[Origin == origin]$Sample] %>%
      which.max() %>%
      names()
  })
  ps1 <- prune_samples(keep.samples, input.ps) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    rename.NA.taxa() %>%
    subset_taxa(
      !is.na(Kingdom) &
        Kingdom != "Eukaryota" &
        Order != "Mitochondria" &
        Order != "Chloroplast" &
        Family != "Mitochondria" &
        Family != "Chloroplast"
    )
  ps1@sam_data$Is_ctrl <- ifelse(is.na(ps1@sam_data$BaP_uM), TRUE, FALSE)
  saveRDS(ps1, file = file.path(output.path, paste0(input.stem, "_prepped.rds")))
  return(ps1)
}

decontam.phyloseq <- function(
    physeq,
    conc.col,
    neg.col,
    decontam.method = "combined",
    output.path
) {
  require(decontam)
  require(phyloseqCompanion)
  require(stringr)
  physeq.stem <- str_replace(deparse(substitute(physeq)), "\\.", "_")

  decontam.res.df <- isContaminant(
    seqtab = physeq,
    method = decontam.method,
    conc = conc.col,
    neg = neg.col,
    batch = physeq@sam_data$Run
  )
  decontam.res.df <- decontam.res.df[sort(row.names(decontam.res.df)), ]
  contam.df <- decontam.res.df[decontam.res.df$contaminant, ]
  saveRDS(
    contam.df,
    file = file.path(
      dirs$save,
      paste0(physeq.stem, "_decontam_results_df.rds")
    )
  )

  freq.plot <- plot_frequency(
    physeq,
    row.names(contam.df),
    conc = conc.col
  ) +
    xlab("DNA Concentration (PicoGreen fluorescent intensity)")
  ggsave(
    freq.plot,
    file = file.path(dirs$plots, "decontam_frequencies_scatterPlot.pdf"),
    width = 14,
    height = 12
  )

  ps.pa <- transform_sample_counts(physeq, function(abund) 1 * (abund > 0))
  ps.pa.neg <- prune_samples(sample_data(ps.pa)[[neg.col]], ps.pa)
  ps.pa.pos <- prune_samples(!sample_data(ps.pa)[[neg.col]], ps.pa)
  # Make data.frame of prevalence in positive and negative samples
  df.pa0 <- data.frame(
    pa.pos = taxa_sums(ps.pa.pos),
    pa.neg = taxa_sums(ps.pa.neg)
  )
  df.pa <- merge(df.pa0, decontam.res.df, by = 0)
  prev.plot <- ggplot(
    data = df.pa,
    aes(x = pa.neg, y = pa.pos, color = contaminant)
  ) +
    geom_quasirandom() +
    xlab("Prevalence (Negative Controls)") +
    ylab("Prevalence (True Samples)")
  ggsave(prev.plot, file = file.path(dirs$plots, "decontam_prevalences_beeswarm.pdf"))

  result.ps <- prune_taxa(!(taxa_names(physeq) %in% row.names(contam.df)), physeq) %>%
    prune_samples(sample_names(physeq)[!physeq@sam_data[[neg.col]]], .) %>%
    prune_taxa(taxa_sums(.) > 0, .)
  result.ps.file <- file.path(
    output.path,
    paste0(physeq.stem, "_decontaminated.rds")
  )
  saveRDS(result.ps, file = result.ps.file)
  return(result.ps)
}

transform.aucs <- function(dt.list, files.list, show.plot = FALSE, force = FALSE) {
  if (!identical(sort(names(dt.list)), sort(names(files.list)))) {
    stop("The names() in the objects passed to `dt.list` and `files.list` must be identical.")
  }
  lapply(setNames(names(dt.list), names(dt.list)), function(dt.name) {
    auc.dt <- dt.list[[dt.name]]
    auc.dt.file <- files.list[dt.name]
    no.zeroes.file <- str_replace(auc.dt.file, "AUCs_", "no0AUCs_")
    only.zeroes.file <- str_replace(auc.dt.file, "AUCs_", "only0AUCs_")
    if (!all(file.exists(c(no.zeroes.file, only.zeroes.file))) | force) {
      if ("Cycle" %in% names(auc.dt)) {
        auc.no0.dt <- lapply(lpr.cycles, function(cycle) {
          subset <- auc.dt[Cycle == cycle & AUC > 0]
          if (show.plot) {
            fitdistrplus::descdist(subset$AUC)
            rcompanion::transformTukey(subset$AUC) %>%
              fitdistrplus::descdist()
          }
          subset[
            , AUC.trans := rcompanion::transformTukey(AUC, plotit = F, quiet = F)
          ]
        }) %>% rbindlist()
        auc.only0.dt <- lapply(lpr.cycles, function(cycle) {
          auc.dt[Cycle == cycle & AUC == 0]
        }) %>% rbindlist()
      } else {
        auc.no0.dt <- auc.dt[AUC > 0]
        if (show.plot) {
          fitdistrplus::descdist(auc.dt$AUC)
          fitdistrplus::descdist(rcompanion::transformTukey(auc.dt$AUC))
        }
        auc.no0.dt[
          , AUC.trans := rcompanion::transformTukey(AUC, plotit = F, quiet = F)
        ]
        auc.only0.dt <- auc.dt[AUC == 0]
      }
      saveRDS(auc.no0.dt, no.zeroes.file)
      cat.n(paste("File written:", no.zeroes.file))
      saveRDS(auc.only0.dt, only.zeroes.file)
      cat.n(paste("File written:", only.zeroes.file))
      results <- list(No0s = auc.no0.dt, Only0s = auc.only0.dt)
    } else {
      cat.n(
        paste(
          "Files:\n\t", no.zeroes.file, "\nt", only.zeroes.file,
          "\nexist and `force = FALSE`; nothing done."
        )
      )
      results <- list(No0s = readRDS(no.zeroes.file), Only0s = readRDS(only.zeroes.file))
    }
    return(results)
  }) %>% return()
}

transform.alpha.scores <- function(base.dt, methods, show.plot = FALSE) {
  alpha.trans.dt <- data.table(
    Method = methods,
    Transformed = FALSE
  ) %>% setkeyv("Method")
  for (alpha in methods) {
    if (show.plot) {
      cat(alpha, sep = "\n")
    }
    alpha.vec <- base.dt[[alpha]]
    if (show.plot) {
      fitdistrplus::descdist(alpha.vec)
    }
    if (shapiro.test(alpha.vec)$p.value <= 0.05) {
      alpha.trans <- rcompanion::transformTukey(
        x = alpha.vec,
        plotit = show.plot,
        quiet = !show.plot
      )
      base.dt[[alpha]] <- alpha.trans
      alpha.trans.dt[alpha]$Transformed <- TRUE
    }
    if (show.plot) {
      fitdistrplus::descdist(alpha.trans)
    }
  }
  # Simpson after transformation is looking like a uniform distribution...
  # not sure how to implement that...
  # Untransformed it's best fit by the beta distribution, and is already between
  # 0 and 1. So let's try that.
  if ("Simpson" %in% names(base.dt)) {
    base.dt$Simpson <- base.dt$Simpson
    alpha.trans.dt["Simpson"]$Transformed <- FALSE
  }
  return(list(Alpha.dt = base.dt, Trans.dt = alpha.trans.dt))
}

# prettify.stat.table <- function(stat.table) {
#   helpers::anova2dt(stat.table) %>%
#     suppressWarnings() %>%
#     as.data.table() %>%
#     set_names(., toTitleCase(names(.))) %>%
#     return()
# }

run.random.forest <- function(
    mod.data,
    response.var,
    rf.type = c("classification", "regression"),
    n.cores,
    rnd.seed = 42
) {
  if (!(rf.type %in% c("classification", "regression"))) {
    stop("`rf.type` must be one of 'classification' or 'regression'")
  }
  mod.recipe <- recipe(mod.data) %>%
    update_role(!!response.var, new_role = "outcome") %>%
    update_role(-!!response.var) %>%
    step_nzv(
      all_predictors(),
      freq_cut = nrow(mod.data) - 1,
      unique_cut = 1 / nrow(mod.data) * 100
    ) %>%
    step_impute_median(all_predictors()) %>%
    step_naomit(!!response.var)
  n.remove <- prep(mod.recipe) %>%
    use_series(steps) %>%
    extract2(1) %>%
    use_series(removals) %>%
    length()
  mtries <- floor((ncol(mod.data) - n.remove - 1) * 0.5)
  while (mtries[1] > 50) {
    mtries <- c(floor(mtries[1] / 2), mtries)
  }
  splitrules <- c("gini", "extratrees")
  minNodeSize <- 1
  if (rf.type == "classification" & length(unique(mod.data[[response.var]])) > 2) {
    train.ctrl <- trainControl(
      method = "cv",
      number = 5,
      summaryFunction = multiClassSummary,
      classProbs = TRUE,
      savePredictions = "all",
      returnResamp = "all",
      verboseIter = TRUE
    )
  } else if (rf.type == "classification") {
    train.ctrl <- trainControl(
      method = "cv",
      number = 5,
      summaryFunction = twoClassSummary,
      classProbs = TRUE,
      savePredictions = "all",
      returnResamp = "all",
      verboseIter = TRUE
    )
  } else {
    splitrules <- c("variance", "extratrees", "maxstat")
    minNodeSize <- 5
    train.ctrl <- trainControl(
      method = "cv",
      number = 5,
      savePredictions = "all",
      returnResamp = "all",
      verboseIter = TRUE
    )
  }

  set.seed(rnd.seed)
  cat(paste("Random seed:", rnd.seed), sep = "\n")
  cl <- makeCluster(n.cores, type = "FORK")
  registerDoParallel(cl, n.cores)
  mod <- caret::train(
    mod.recipe,
    data = mod.data,
    method = "ranger",
    importance = "permutation",
    tuneGrid = expand.grid(
      mtry = mtries,
      splitrule = splitrules,
      min.node.size = minNodeSize
    ),
    trControl = train.ctrl
  ) %>% try(silent = TRUE)
  stopCluster(cl)
  if ("try-error" %in% class(mod)) {
    cat(mod)
  } else {
    return(mod)
  }
}

id.sig.important.features <- function(
    rf.model,
    response.var,
    n.cores,
    n.perm = 100,
    rnd.seed = 42
) {
  require(rlist)

  x <- rf.model$finalModel
  data <- as.data.frame(rf.model$trainingData)
  data <- prep(rf.model$recipe, training = rf.model$trainingData) %>%
    bake(rf.model$trainingData) %>%
    as.data.frame()
  dat_x <- data[, -c(which(names(data) == response.var))]
  dat_x <- dat_x[, names(x$variable.importance)] %>% as.matrix()
  set.seed(rnd.seed)
  cl <- makeCluster(n.cores, type = "FORK")
  registerDoParallel(cl, n.cores)
  vimp <- foreach(
    i = 1:n.perm,
    .final = rlist::list.cbind,
    .verbose = TRUE
  ) %dopar%
    {
      dat_y <- data[sample(nrow(data)), response.var]
      if (class(dat_y) == "character") {
        dat_y <- factor(dat_y)
      }
      ranger::ranger(
        x = dat_x,
        y = dat_y,
        num.trees = x$num.trees,
        mtry = x$mtry,
        min.node.size = x$min.node.size,
        importance = x$importance.mode,
        replace = x$replace,
        num.threads = 1
      )$variable.importance
    } %>% try(silent = TRUE)
  stopCluster(cl)
  if ("try-error" %in% class(vimp)) {
    cat(vimp)
  } else {
    pval <- sapply(1:nrow(vimp), function(i) {
      (sum(vimp[i, ] >= x$variable.importance[i]) + 1) / (ncol(vimp) + 1)
    })

    impt.mat <- cbind(x$variable.importance, pval)
    colnames(impt.mat) <- c("Importance", "P.value")

    sig.impt.dt <- impt.mat[impt.mat[, "P.value"] <= 0.05, ] %>%
      as.data.table(keep.rownames = "Taxon")
    sig.impt.dt <- sig.impt.dt[order(Importance, decreasing = T)]
    return(sig.impt.dt)
  }
}

parse.chunks <- function(RMD) {
  require(stringr)
  require(readr)
  if (!dir.exists("Chunks")) {
    dir.create("Chunks")
  }
  if (length(list.files("Chunks")) > 0) {
    file.remove(list.files("Chunks", full.names = T))
  }
  lines <- readLines(RMD)
  while (!str_detect(lines[1], "^```\\{r")) {
    lines <- lines[-1]
  }
  curr.file <- NULL
  curr.chunk <- 0
  write.lines <- TRUE
  for (line in lines) {
    if (str_detect(line, "^```\\{r") & !str_detect(line, "eval=FALSE")) {
      curr.chunk <- curr.chunk + 1
      curr.file <- str_extract(line, "r [\\w\\d\\-]+") %>%
        str_remove("^r ") %>%
        str_replace_all("-", "_") %>%
        paste0(
          "Chunks/", ifelse(curr.chunk < 10, paste0("0", curr.chunk), curr.chunk), "_", ., ".R"
        )
      file.create(curr.file)
      write.lines <- TRUE
    } else if (str_detect(line, "^```$")) {
      write.lines <- FALSE
    } else if (write.lines) {
      write_lines(x = line, file = curr.file, append = T)
    }
  }
}

gen.clr.matrix <- function(
    asv.mat,
    min_reads,
    min_prop = 0.001,
    min_occur = 0,
    smpls_by_row = TRUE,
    method = "CZM",
    lab = 0
) {
  require(CoDaSeq)
  require(zCompositions)
  asv.mat.f <- codaSeq.filter(
    asv.mat,
    min.reads = min_reads,
    min.prop = min_prop,
    min.occurrence = min_occur,
    samples.by.row = smpls_by_row
  )
  asv.mat.f <- asv.mat.f[, colSums(asv.mat.f) > 0]
  # replace 0 values with an estimate
  asv.mat.fn0 <- cmultRepl(t(asv.mat.f), method = method, label = lab)
  return(codaSeq.clr(asv.mat.fn0))
}

rename.taxa.dt <- function(dt, lvls) {
  renamed.dt <- data.table(Phylum = paste0("p__", dt[["Phylum"]]))
  for (l in 2:length(lvls)) {
    level <- lvls[l]
    renamed.dt[[level]] <- paste(
      renamed.dt[[l - 1]],
      paste0(names(lvls[lvls == level]), dt[[level]]),
      sep = "."
    )
  }
  renamed.dt[["Taxon"]] <- dt$Taxon

  return(renamed.dt)
}

generate.auc.rf.data <- function(
    physeq,
    split.strata,
    lpr.cycles = c("light", "dark")
) {
  require(stringr)
  levels <- setNames(
    c("Phylum", "Class", "Order", "Family", "Genus"),
    c("p__", "c__", "o__", "f__", "g__")
  )
  return.list <- NULL
  sample.dt <- sample.data.table(physeq)
  asv.dt <- otu.data.table(physeq) %>% setkeyv(smpl.col)
  taxa.dt <- taxa.data.table(physeq) %>% setkeyv("Taxon")
  renamed.taxa.dt <- rename.taxa.dt(taxa.dt, levels)
  asv.cols <- names(asv.dt)[grepl("ASV", names(asv.dt))]
  noSample.dt <- asv.dt[, ..asv.cols]
  agg.dt <- copy(asv.dt)
  for (level in levels) {
    lvl.dt <- sapply(sort(unique(renamed.taxa.dt[[level]])), function(assignment) {
      asvs <- renamed.taxa.dt[renamed.taxa.dt[[level]] == assignment]$Taxon
      return(rowSums(noSample.dt[, ..asvs]))
    }) %>% as.data.table()
    agg.dt <- cbind(agg.dt, lvl.dt)
  }

  for (cycle in lpr.cycles) {
    keep.cols <- c("Sample", split.strata, str_subset(names(sample.dt), paste0("LPR.", cycle)))
    auc.asv.dt <- sample.dt[, ..keep.cols] %>%
      merge(asv.dt, by = "Sample")
    auc.asv.dt[, Sample := NULL]
    set.seed(user.seed)
    asv.split <- initial_split(auc.asv.dt, prop = 0.7, strata = all_of(split.strata))
    return.list$ASV_only[[cycle]]$Train <- training(asv.split)
    return.list$ASV_only[[cycle]]$Test <- testing(asv.split)

    auc.agg.dt <- sample.dt[, ..keep.cols] %>%
      merge(agg.dt, by = "Sample")
    auc.agg.dt[, Sample := NULL]
    set.seed(user.seed)
    agg.split <- initial_split(auc.agg.dt, prop = 0.7, strata = all_of(split.strata))
    return.list$aggregated[[cycle]]$Train <- training(agg.split)
    return.list$aggregated[[cycle]]$Test <- testing(agg.split)
  }
  return(return.list)
}

aggregate.taxa <- function(physeq, verbose = FALSE) {
  levels <- setNames(
    c("Phylum", "Class", "Order", "Family", "Genus"),
    c("p__", "c__", "o__", "f__", "g__")
  )
  asv.dt <- otu.data.table(physeq) %>% setkeyv(smpl.col)
  taxa.dt <- taxa.data.table(physeq) %>% setkeyv("Taxon")
  renamed.taxa.dt <- rename.taxa.dt(taxa.dt, levels)
  asv.cols <- names(asv.dt)[grepl("ASV", names(asv.dt))]
  noSample.dt <- asv.dt[, ..asv.cols]
  agg.dt <- copy(asv.dt)
  for (level in levels) {
    if (verbose) {
      cat.n(level)
    }
    lvl.dt <- sapply(sort(unique(renamed.taxa.dt[[level]])), function(assignment) {
      asvs <- renamed.taxa.dt[renamed.taxa.dt[[level]] == assignment]$Taxon
      return(rowSums(noSample.dt[, ..asvs]))
    }) %>% as.data.table()
    agg.dt <- cbind(agg.dt, lvl.dt)
  }
  return(agg.dt)
}

generate.bap.rf.data <- function(
    physeq,
    split.strata
) {
  require(stringr)
  return.list <- NULL
  asv.dt <- otu.data.table(physeq) %>% setkeyv(smpl.col)
  agg.dt <- aggregate.taxa(physeq, verbose = T) %>% setkeyv(smpl.col)
  smpl.dt <- sample.data.table(physeq)
  bap.asv.dt <- merge(smpl.dt[, .(Sample, BaP_uM)], asv.dt, by = "Sample")
  bap.asv.dt[, Sample := NULL]
  set.seed(user.seed)
  asv.split <- initial_split(bap.asv.dt, prop = 0.7, strata = all_of(split.strata))
  return.list$ASV_only$Train <- training(asv.split)
  return.list$ASV_only$Test <- testing(asv.split)

  bap.agg.dt <- merge(smpl.dt[, .(Sample, BaP_uM)], agg.dt, by = "Sample")
  bap.agg.dt[, Sample := NULL]
  set.seed(user.seed)
  agg.split <- initial_split(bap.agg.dt, prop = 0.7, strata = all_of(split.strata))
  return.list$aggregated$Train <- training(agg.split)
  return.list$aggregated$Test <- testing(agg.split)
  return(return.list)
}

find.bad.batches <- function(ps, metric = "Chao1", threshhold = 90) {
  sample.dt <- sample.data.table(
    ps,
    sample.column.name = "Sequencing.ID"
  )
  tax.dt <- estimate_richness(
    physeq = ps,
    measures = alpha.methods[1:3]
  ) %>%
    as.data.table(keep.rownames = "Sequencing.ID") %>%
    setkeyv("Sequencing.ID")
  tax.dt[, se.chao1 := NULL]

  phy.dt <- pd(
    samp = otu.matrix(ps),
    tree = phy_tree(ps)
  ) %>%
    as.data.table(keep.rownames = "Sequencing.ID") %>%
    setkeyv("Sequencing.ID")
  names(phy.dt)[2:3] <- alpha.methods[4:5]
  alpha.base <- tax.dt[phy.dt]

  alpha.prep.dt <- merge(
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
    ggplot(alpha.prep.dt, aes(x = as.numeric(Extract_batch), y = Score)) +
      geom_quasirandom() +
      stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", color = "red") +
      facet_wrap(~Alpha, scales = "free_y")
    } %>%
    ggsave(
      filename = file.path(
        dirs$plots,
        "QC_alphaDiv_by_extractionBatch_psPrepped_beeswarms.pdf"
      ),
      width = 14,
      height = 12
    )
  batch.chao1s <- alpha.prep.dt[
    Alpha == metric,
    .(mean.chao1 = mean(Score)),
    keyby = Extract_batch
  ]
  return(batch.chao1s[mean.chao1 > threshhold]$Extract_batch)
}

save.stat.tables <- function(flex, filename, exts = c(".png", ".html"), ...) {
  exts <- match.arg(
    exts,
    choices = c(".png", ".pdf", ".jpeg", ".docx", ".html", ".pptx"),
    several.ok = TRUE
  )
  to.save <- lapply(exts, function(ext) {
    if (ext %in% c(".png", ".pdf", ".jpeg")) {
      save_as_image(flex, path = paste0(filename, ext), webshot = "webshot2", ...)
    } else if (ext == ".docx") {
      save_as_docx(flex, path = paste0(filename, ext))
    } else if (ext == ".html") {
      save_as_html(flex, path = paste0(filename, ext))
    } else if (ext == ".pptx") {
      save_as_pptx(flex, path = paste0(filename, ext))
    } else {
      rlang::abort(message = "Exention choice(s) not recgonized")
    }
  })
}

print.html.table <- function(file) {
  XML::htmlParse(file = file) %>%
    print()
}

run.lme <- function(fixed.frm, random.frm, mod.data) {
  base.args <- list(
    na.action = na.omit,
    method = "ML",
    control = lmeControl(opt = "optim")
  )
  all.args <- c(
    list(
      fixed = fixed.frm,
      random = random.frm,
      data = as.name(deparse(substitute(mod.data)))
    ),
    base.args
  )
  do.call(nlme::lme, args = all.args, quote = TRUE) %>%
    return()
}

# Adam Burns - 2/10/2015
# aburns2@uoregon.edu
# From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
# Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE. For use in R.
# spp: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
# pool: A community table for defining source community (optional; Default=NULL).
# taxon: A table listing the taxonomic calls for each otu, with OTU ids as row names and taxonomic classifications as columns.
# If stats=TRUE the function will return fitting statistics.
# If stats=FALSE the function will return a table of observed and predicted values for each otu.

sncm.fit <- function(spp, pool = NULL, stats = TRUE, taxon = NULL) {
  require(minpack.lm)
  require(Hmisc)
  require(stats4)

  options(warn = -1)

  # Calculate the number of individuals per community
  N <- mean(apply(spp, 1, sum))

  # Calculate the average relative abundance of each taxa across communities
  if (is.null(pool)) {
    p.m <- apply(spp, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m / N
  } else {
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m / N
  }

  # Calculate the occurrence frequency of each taxa across communities
  spp.bi <- 1 * (spp > 0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]

  # Combine
  C <- merge(p, freq, by = 0)
  C <- C[order(C[, 2]), ]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ] # Removes rows with any zero (absent in either source pool or local communities)
  p <- C.0[, 2]
  freq <- C.0[, 3]
  names(p) <- C.0[, 1]
  names(freq) <- C.0[, 1]

  # Calculate the limit of detection
  d <- 1 / N

  ## Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(
    freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE),
    start = list(m = 0.1)
  )
  m.ci <- confint(m.fit, "m", level = 0.95)

  ## Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  sncm.LL <- function(m, sigma) {
    R <- freq - pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE)
    R <- dnorm(R, 0, sigma)
    - sum(log(R))
  }
  m.mle <- mle(
    sncm.LL,
    start = list(m = 0.1, sigma = 0.1),
    nobs = length(p)
  ) %>% try(silent = T)
  if ("try-error" %in% class(m.mle)) {
    m.mle <- NULL
  }

  ## Calculate Akaike's Information Criterion (AIC)
  if (is.null(m.mle)) {
    aic.fit <- NULL
    bic.fit <- NULL
  } else {
    aic.fit <- AIC(m.mle, k = 2)
    bic.fit <- BIC(m.mle)
  }

  ## Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2)) / (sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq - freq.pred)^2) / (length(freq) - 1))

  pred.ci <- binconf(
    freq.pred * nrow(spp),
    nrow(spp),
    alpha = 0.05,
    method = "wilson",
    return.df = TRUE
  )

  ## Calculate AIC for binomial model
  bino.LL <- function(mu, sigma) {
    R <- freq - pbinom(d, N, p, lower.tail = FALSE)
    R <- dnorm(R, mu, sigma)
    - sum(log(R))
  }
  bino.mle <- mle(bino.LL, start = list(mu = 0, sigma = 0.1), nobs = length(p))

  aic.bino <- AIC(bino.mle, k = 2)
  bic.bino <- BIC(bino.mle)

  ## Goodness of fit for binomial model
  bino.pred <- pbinom(d, N, p, lower.tail = FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2)) / (sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2) / (length(freq) - 1))

  bino.pred.ci <- binconf(
    bino.pred * nrow(spp),
    nrow(spp),
    alpha = 0.05,
    method = "wilson",
    return.df = TRUE
  )

  ## Calculate AIC for Poisson model
  pois.LL <- function(mu, sigma) {
    R <- freq - ppois(d, N * p, lower.tail = FALSE)
    R <- dnorm(R, mu, sigma)
    - prep.behavior.data.noPSsum(log(R))
  }
  pois.mle <- mle(pois.LL, start = list(mu = 0, sigma = 0.1), nobs = length(p))

  aic.pois <- AIC(pois.mle, k = 2)
  bic.pois <- BIC(pois.mle)

  ## Goodness of fit for Poisson model
  pois.pred <- ppois(d, N * p, lower.tail = FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2)) / (sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2) / (length(freq) - 1))

  pois.pred.ci <- binconf(
    pois.pred * nrow(spp),
    nrow(spp),
    alpha = 0.05,
    method = "wilson",
    return.df = TRUE
  )

  ## Results
  if (stats == TRUE) {
    fitstats <- data.frame(
      m = numeric(),
      m.ci = numeric(),
      m.mle = numeric(),
      maxLL = numeric(),
      binoLL = numeric(),
      poisLL = numeric(),
      Rsqr = numeric(),
      Rsqr.bino = numeric(),
      Rsqr.pois = numeric(),
      RMSE = numeric(),
      RMSE.bino = numeric(),
      RMSE.pois = numeric(),
      AIC = numeric(),
      BIC = numeric(),
      AIC.bino = numeric(),
      BIC.bino = numeric(),
      AIC.pois = numeric(),
      BIC.pois = numeric(),
      N = numeric(),
      Samples = numeric(),
      Richness = numeric(),
      Detect = numeric()
    )
    fitstats[1, ] <- c(
      coef(m.fit),
      coef(m.fit) - m.ci[1],
      ifelse(is.null(m.mle), NA, m.mle@coef["m"]),
      ifelse(is.null(m.mle), NA, m.mle@details$value),
      bino.mle@details$value,
      pois.mle@details$value,
      Rsqr,
      Rsqr.bino,
      Rsqr.pois,
      RMSE,
      RMSE.bino,
      RMSE.pois,
      ifelse(is.null(aic.fit), NA, aic.fit),
      ifelse(is.null(bic.fit), NA, bic.fit),
      aic.bino,
      bic.bino,
      aic.pois,
      bic.pois,
      N,
      nrow(spp),
      length(p),
      d
    )
    return(fitstats)
  } else {
    A <- cbind(p, freq, freq.pred, pred.ci[, 2:3], bino.pred, bino.pred.ci[, 2:3])
    A <- as.data.frame(A)
    colnames(A) <- c(
      "p", "freq", "freq.pred", "pred.lwr", "pred.upr", "bino.pred", "bino.lwr", "bino.upr"
    )
    if (is.null(taxon)) {
      B <- A[order(A[, 1]), ]
    } else {
      B <- merge(A, taxon, by = 0, all = TRUE)
      row.names(B) <- B[, 1]
      B <- B[, -1]
      B <- B[order(B[, 1]), ]
    }
    return(B)
  }
}
