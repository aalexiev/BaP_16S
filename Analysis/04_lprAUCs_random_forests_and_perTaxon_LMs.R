# 04_lprAUCs_random_forests.R

if (!("dirs" %in% ls())) {
  source("Analysis/00_setup.R")
}

script.var <- "lpr.auc.rfs"
script.var.state <- get.redo.state(script.var)
if (is.null(script.var.state)) {
  assign.redo(script.var, state = default.redo.state)
}
sample.dt <- copy(sample.data.table(ps.processed))

lpr.rf.data.list.file <- file.path(
  dirs$save,
  "lprAUCs_rf_data_splits_list.rds"
)
lpr.rf.data.list <- redo.if(script.var, lpr.rf.data.list.file, {
  generate.auc.rf.data(
    physeq = ps.processed,
    split.strata = "BaP_factor"
  )
})

lpr.rf.list <- lapply(taxa.sets, function(set) {
  lapply(lpr.cycles, function(cycle) {
    # cat.n(paste(set, "-", cycle))
    filename.parts <- c("rf_lprAUCs_", cycle, set, "_train.rds")
    rf.data <- copy(
      lpr.rf.data.list[[set]][[cycle]]$Train
    )
    rf.data[, BaP_factor := NULL]
    rf.file <- file.path(dirs$save, paste(filename.parts, collapse = "_"))
    rf <- redo.if(script.var, rf.file, {
      run.random.forest(
        mod.data = rf.data[complete.cases(rf.data)],
        response.var = str_subset(names(sample.dt), paste0("LPR.", cycle)),
        rf.type = "regression",
        n.cores = nCores
      )
    })
    return(rf)
  })
})

lpr.rmses.dt.file <- file.path(
  dirs$save, "lprAUCs_rf_importantFeatures_dt.rds"
)
lpr.rmses.dt <- redo.if(script.var, lpr.rmses.dt.file, {
  lapply(taxa.sets, function(set) {
    lapply(lpr.cycles, function(cycle) {
      rf.train <- lpr.rf.list[[set]][[cycle]]
      test.data <- copy(
        lpr.rf.data.list[[set]][[cycle]]$Test
      )
      test.data[, BaP_factor := NULL]
      test.data <- test.data[complete.cases(test.data)]
      preds <- predict(rf.train, test.data)
      response.col <- rf.train$recipe$var_info %>%
        filter(role == "outcome") %>%
        select(variable) %>%
        as.character()
      return(
        data.table(
          Taxa.set = set,
          Cycle = cycle,
          RMSE = sqrt(mean((preds - test.data[[response.col]])^2))
        )
      )
    }) %>% rbindlist()
  }) %>% rbindlist()
})

best.lpr.rf.datasets <- lpr.rmses.dt[
  , .(RMSE = min(RMSE)),
  by = "Cycle"
] %>%
  merge(lpr.rmses.dt[, .(RMSE, Taxa.set)], by = "RMSE")

lpr.impt.dts.list <- lapply(lpr.cycles, function(cycle) {
  taxa.set <- best.lpr.rf.datasets[Cycle == cycle]$Taxa.set
  impt.dt.file <- file.path(dirs$save, paste0("lprAUCs_rf_", cycle,"_importantFeatures.rds"))
  impt.dt <- redo.if(script.var, impt.dt.file, {
    id.sig.important.features(
      rf.model = lpr.rf.list[[taxa.set]][[cycle]],
      response.var = str_subset(names(sample.dt), paste0("LPR.", cycle)),
      n.cores = nCores,
      rnd.seed = user.seed
    )
  })
  return(impt.dt)
})

lpr.impt.gt0.dts.list.file <- file.path(
  dirs$save,
  "lprAUCs_rf_importantFeatures_importanceGt0_dts_list.rds"
)
lpr.impt.gt0.dts.list <- redo.if(script.var, lpr.impt.gt0.dts.list.file, {
  lapply(lpr.impt.dts.list, function(dt) {
    return(dt[Importance > 0])
  })
})

####

bap.lpr.rf.list <- lapply(taxa.sets, function(set) {
  lapply(lpr.cycles, function(cycle) {
    lapply(rf.types, function(type) {
      # cat.n(paste(set, "-", cycle))
      filename.parts <- c("rf_bap_lprAUCs", type, cycle, set, "train.rds")
      rf.data <- copy(
        lpr.rf.data.list[[set]][[cycle]]$Train
      )
      if (type == "classification") {
        rf.data[, BaP_factor := factor(paste0("uM", BaP_factor), levels = paste0("uM", bap.concs))]
      } else {
        rf.data[, BaP_uM := as.numeric(as.character(BaP_factor))]
        rf.data[, BaP_factor := NULL]
      }
      rf.file <- file.path(dirs$save, paste(filename.parts, collapse = "_"))
      rf <- redo.if(script.var, rf.file, {
        run.random.forest(
          mod.data = rf.data[complete.cases(rf.data)],
          response.var = ifelse(type == "classification", "BaP_factor", "BaP_uM"),
          rf.type = type,
          n.cores = nCores
        )
      })
      return(rf)
    })
  })
})

bap.lpr.rocAUCs.dt.file <- file.path(dirs$save, "bap_lpr_rf_roc_aucs_dt.rds")
bap.lpr.rocAUCs.dt <- redo.if(script.var, bap.lpr.rocAUCs.dt.file, {
  lapply(taxa.sets, function(set) {
    lapply(lpr.cycles, function(cycle) {
      rf.train <- bap.lpr.rf.list[[set]][[cycle]]$classification
      test.data <- copy(lpr.rf.data.list[[set]][[cycle]]$Test)
      test.data[, BaP_factor := factor(paste0("uM", BaP_factor), levels = paste0("uM", bap.concs))]
      preds <- predict(rf.train, test.data, type = "prob")
      roc <- multiclass.roc(test.data$BaP_factor, preds)
      return(
        data.table(
          Taxa.set = set,
          Cycle = cycle,
          ROC.AUC = as.numeric(auc(roc))
        )
      )
    }) %>% rbindlist()
  }) %>% rbindlist()
})

bap.lpr.rmses.dt.file <- file.path(dirs$save, "bap_lpr_rf_rmses_dt.rds")
bap.lpr.rmses.dt <- redo.if(script.var, bap.lpr.rmses.dt.file, {
  lapply(taxa.sets, function(set) {
    lapply(lpr.cycles, function(cycle) {
      rf.train <- bap.lpr.rf.list[[set]][[cycle]]$regression
      test.data <- copy(lpr.rf.data.list[[set]][[cycle]]$Test)
      test.data[, BaP_uM := as.numeric(as.character(BaP_factor))]
      test.data[, BaP_factor := NULL]
      preds <- predict(rf.train, test.data)
      return(
        data.table(
          Taxa.set = set,
          Cycle = cycle,
          RMSE = sqrt(mean((preds - test.data$BaP_uM)^2))
        )
      )
    }) %>% rbindlist()
  }) %>% rbindlist()
})

lapply(lpr.cycles, function(cycle) {
  train.data <- lpr.rf.data.list$ASV_only[[cycle]]$Train[, 1:2]
  train.data <- train.data[complete.cases(train.data)]
  frm <- paste("BaP_factor ~", names(train.data)[2])
  mod <- nnet::multinom(as.formula(frm), data = train.data )
  test.data <- lpr.rf.data.list$ASV_only[[cycle]]$Test[, 1:2]
  test.data <- test.data[complete.cases(test.data)]
  preds <- predict(mod, test.data, type = "prob")
  roc <- multiclass.roc(test.data$BaP_factor, preds)
  data.table(
    Cycle = cycle,
    ROC.AUC = as.numeric(auc(roc))
  ) %>% return()
}) %>% rbindlist()

####

remove.pattern <- "[-/\\(\\)\\.]"
lpr.sig.lms.file <- file.path(
  dirs$save,
  "lprAUCs_sig_lme_and_mediation_dts_list.rds"
)
lpr.sig.lms <- redo.if(script.var, lpr.sig.lms.file, {
  res.list <- NULL
  bap.var <- "BaP_uM"
  random.effect <- "SARL_plateID"
  random.eff.frm <- as.formula(paste0("~ ", bap.var, "|", random.effect))
  for (cycle in lpr.cycles) {
    cat.n(paste("###", cycle, "###"))
    auc.var <- str_subset(names(sample.dt), paste0("LPR.", cycle))
    sig.taxa <- lpr.impt.gt0.dts.list[[cycle]]$Taxon %>% unique()
    lpr.taxa.set <- best.lpr.rf.datasets[Cycle == cycle]$Taxa.set

    smpl.dt.cols <- c("Sample", auc.var, "BaP_factor", "BaP_uM", "SARL_plateID")
    lm.data <- rbind(
      copy(lpr.rf.data.list[[lpr.taxa.set]][[cycle]]$Train),
      copy(lpr.rf.data.list[[lpr.taxa.set]][[cycle]]$Test)
    ) %>% merge(
      sample.dt[, ..smpl.dt.cols],
      by = c(auc.var, "BaP_factor")
    )

    sig.taxa <- sig.taxa[sig.taxa %in% names(lm.data)]

    lm.res <- data.table()
    lm.pred <- data.table()
    med.res <- data.table()
    for (taxon in sig.taxa) {
      cat.n(taxon)
      old.id <- taxon
      if (str_detect(taxon, remove.pattern)) {
        new.id <- str_remove_all(taxon, remove.pattern)
        names(lm.data)[which(names(lm.data) == taxon)] <- new.id
        taxon <- new.id
      }

      ## Mixed Effects Linear Regressions
      frms <- NULL
      mods <- NULL
      select.mod.dt <- NULL
      frms[[1]] <- paste0(auc.var, " ~ ", bap.var, " + I(", bap.var, "^2)")
      frms[[2]] <- paste(frms[[1]], "+", taxon)
      frms[[3]] <- paste(frms[[1]], "*", taxon)
      for (m in 1:3) {
        mods[[m]] <- run.lme(
          fixed.frm = as.formula(frms[[m]]),
          random.frm = random.eff.frm,
          mod.data = lm.data
        ) %>% try(silent = TRUE)
      }
      # sapply(mods, class)
      mods <- mods[sapply(mods, function(m) class(m) != "try-error")]
      if (length(mods) > 1) {
        if (length(mods) == 3) {
          mod.anova <- anova.lme(mods[[1]], mods[[2]], mods[[3]]) %>%
            as.data.table() %>%
            try(silent = TRUE)
        } else if (length(mods) == 2) {
          mod.anova <- anova(mods[[1]], mods[[2]]) %>%
            as.data.table() %>%
            try(silent = TRUE)
        }
        # print(mod.anova)
        mod.anova$`p-value`[1] <- 0
        mod.choice <- max(which(mod.anova$`p-value` <= 0.05))
        mod <- mods[[mod.choice]]
      } else {
        mod.choice <- 1
      }
      cat.n(paste("model choice:", mod.choice))
      if (mod.choice > 1) {
        mod.pval <- mod.anova$`p-value`[mod.choice]
      } else {
        mod0 <- run.lme(
          fixed.frm = as.formula(paste(auc.var, "~ 1")),
          random.frm = random.eff.frm,
          mod.data = lm.data
        ) %>% try(silent = TRUE)
        mod.pval <- anova(mod0, mods[[1]]) %>%
          as.data.table() %>%
          use_series(`p-value`) %>%
          magrittr::extract(2) %>%
          try(silent = TRUE)
      }

      select.mod.dt <- tidy(mod) %>% as.data.table()
      select.mod.dt[
        , `:=`(Taxon = old.id, Mod.pval = mod.pval)
      ]
      names(select.mod.dt) <- toTitleCase(names(select.mod.dt))
      lm.res <- rbind(lm.res, select.mod.dt)

      ## Model predictions
      predictions.dt <- NULL
      if (mod.choice > 1) {
        predictions.dt <- lapply(bap.concs, function(conc) {
          pred.dt <- data.table(
            Taxon = seq(min(lm.data[[taxon]]), max(lm.data[[taxon]]), length = 1000),
            BaP_uM = conc
          )
          names(pred.dt)[1] <- taxon
          pred.dt[, (auc.var) := predict(mod, newdata = pred.dt, level = 0)]
          names(pred.dt)[1] <- "Abundance"
          pred.dt[, Taxon := taxon]
        }) %>% rbindlist()
      }
      lm.pred <- rbind(lm.pred, predictions.dt)

      # Mediation
      keep.cols <- c(taxon, bap.var, auc.var)
      mod.data <- lm.data[complete.cases(lm.data[, ..keep.cols])]
      frm.med <- paste0(taxon, " ~ ", bap.var, " + I(", bap.var, "^2)")
      res.sum <- mediate(
        glm(as.formula(frm.med), data = mod.data, family = poisson(link = "sqrt")),
        glm(as.formula(frms[[2]]), data = mod.data, family = "gaussian"),
        treat = bap.var,
        mediator = taxon,
        boot = TRUE,
        sims = 500
      ) %>%
        summary() %>%
        try(silent = TRUE)
      if ("try-error" %in% class(res.sum)) {
        med.dt <- NULL
      } else {
        prop <- c(2, 3, 1, 4)
        p.vals <- unlist(res.sum[str_subset(names(res.sum), "[gu]\\.p$")])[prop]
        if (any(p.vals <= 0.05)) {
          med.dt <- data.table(
            Taxon = old.id,
            Term = c("ACME", "ADE", "Total Effect", "Prop. Mediated"),
            Estimate = unlist(res.sum[str_subset(names(res.sum), "avg$|coef$")][prop]),
            CI.95.Lo = res.sum[str_subset(names(res.sum), "[gu]\\.ci$")][prop] %>%
              sapply(`[`, 1),
            CI.95.Up = res.sum[str_subset(names(res.sum), "[gu]\\.ci$")][prop] %>%
              sapply(`[`, 2),
            P.value = p.vals
          )
          med.res <- rbind(med.res, med.dt)
        }
      }
    }
    res.list[[cycle]]$LME <- lm.res
    res.list[[cycle]]$PRD <- lm.pred
    res.list[[cycle]]$MED <- med.res
  }
  res.list
})

lm.sig.dt.file <- file.path(
  dirs$save,
  "lprAUCs_taxonLMs_significant_associations_all_dt.rds"
)
lm.sig.dt <- lapply(lpr.cycles, function(cycle) {
  dt <- lpr.sig.lms[[cycle]]$LME
  p.dt <- dt[, .(Mod.pval = Mod.pval[1]), by = Taxon]
  dt[, `:=`(Cycle = cycle)]
  dt[Taxon %in% p.dt[Mod.pval <= 0.05]$Taxon] %>%
    return()
}) %>% rbindlist()

interactions.taxa <- lapply(lpr.cycles, function(cycle) {
  unique(lm.sig.dt[str_detect(Term, ":") & Cycle == cycle]$Taxon)
})
main.effects.taxa <- lapply(lpr.cycles, function(cycle) {
  taxa <- unique(
    lm.sig.dt[Effect == "fixed" & !str_detect(Term, "BaP") & Cycle == cycle]$Taxon
  )
  taxa[!(taxa %in% interactions.taxa[[cycle]])] %>% return()
})
sig.intrxn.dt <- lm.sig.dt[Taxon %in% interactions.taxa]
sig.mainEf.dt <- lm.sig.dt[Taxon %in% main.effects.taxa]

lpr.sig.mediation <- lapply(lpr.cycles, function(cycle) {
  dt <- lpr.sig.lms[[cycle]]$MED
}) %>% rbindlist()
