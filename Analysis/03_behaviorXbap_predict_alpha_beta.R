# 03_bapXbehav_predict_alpha_beta.R

if (!("dirs" %in% ls())) {
  source("Analysis/00_setup.R")
}

script.var <- "behavior.by.bap"
script.var.state <- get.redo.state(script.var)
if (is.null(script.var.state)) {
  assign.redo(script.var, state = default.redo.state)
}
sample.dt <- copy(sample.data.table(ps.processed))
### Alpha-diversity by BaP exposure and movement AUCs
### EPR
alpha.epr.lmes.dt.file <- file.path(dirs$save, "bapXbehav_alpha_epr_lmes_dt.rds")
alpha.epr.lmes.dt <- redo.if(script.var, alpha.epr.lmes.dt.file, {
  lapply(alpha.methods, function(alpha) {
    lms <- NULL
    predictors <- c("BaP_uM", alpha)
    model.data <- sample.dt[complete.cases(sample.dt[, ..predictors])]
    random.effect <- str_subset(names(sample.dt), "SARL")
    .GlobalEnv$fixed1 <- paste("EPR.auc ~", predictors[1])
    .GlobalEnv$randEff <- paste("~1 |", random.effect)
    lms[[1]] <- nlme::lme(
      data = model.data,
      fixed = as.formula(fixed1),
      random = as.formula(randEff),
      na.action = na.omit,
      method = "ML",
      control = lmeControl(opt = 'optim')
    )
    .GlobalEnv$fixed2 <- paste(fixed1, "+", predictors[2])
    lms[[2]] <- update(lms[[1]], fixed = as.formula(fixed2))
    .GlobalEnv$fixed3 <- paste(fixed1, "*", predictors[2])
    lms[[3]] <- update(lms[[1]], fixed3)
    res.dt <- anova(lms[[1]], lms[[2]], lms[[3]]) %>%
      as.data.table() %>%
      set_names(toTitleCase(str_replace_all(names(.), "-", ".")))
    if (all(res.dt$P.value[2:3] > 0.05)) {
      mod <- lms[[1]]
    } else if (res.dt$P.value[3] <= 0.05) {
      mod <- lms[[3]]
    } else {
      mod <- lms[[2]]
    }

    prettify.stat.table(mod) %>%
      .[
        , `:=`(
          Sig = ifelse(P.value <= 0.05, "*", ""),
          Alpha.metric = alpha
        )
      ] %>%
      setcolorder(neworder = ncol(.)) %>%
      return()
  }) %>% rbindlist()
})

### LPR
alpha.lpr.lmes.dt.file <- file.path(dirs$save, "bapXbehav_alpha_lpr_lmes_dt.rds")
alpha.lpr.lmes.dt <- redo.if(script.var, alpha.lpr.lmes.dt.file, {
  lapply(alpha.methods, function(alpha) {
    lapply(lpr.cycles, function(cycle) {
      lms <- NULL
      response <- str_subset(names(sample.dt), paste0("LPR.", cycle))
      predictors <- c("BaP_uM", alpha)
      model.data <- sample.dt[complete.cases(sample.dt[, ..predictors])]
      random.effect <- str_subset(names(sample.dt), "SARL")
      .GlobalEnv$fixed1 <- paste(response, "~", predictors[1])
      .GlobalEnv$randEff <- paste("~1 |", random.effect)
      lms[[1]] <- nlme::lme(
        data = model.data,
        fixed = as.formula(fixed1),
        random = as.formula(randEff),
        na.action = na.omit,
        method = "ML",
        control = lmeControl(opt = 'optim')
      )
      .GlobalEnv$fixed2 <- paste0(fixed1, " + I(", predictors[1], "^2)")
      lms[[2]] <- update(lms[[1]], fixed = as.formula(fixed2))
      aov1 <- anova(lms[[1]], lms[[2]]) %>% as.data.table()
      if (aov1$`p-value`[2] <= 0.05) {
        base.lm <- lms[[2]]
        .GlobalEnv$fixed3 <- paste(fixed2, "+", predictors[2])
        lms[[3]] <- update(lms[[1]], fixed = as.formula(fixed3))
        .GlobalEnv$fixed4 <- paste(fixed2, "*", predictors[2])
        lms[[4]] <- update(lms[[1]], fixed = as.formula(fixed4))
      } else {
        base.lm <- lms[[1]]
        .GlobalEnv$fixed3 <- paste(fixed1, "+", predictors[2])
        lms[[3]] <- update(lms[[1]], fixed = as.formula(fixed3))
        .GlobalEnv$fixed4 <- paste(fixed1, "*", predictors[2])
        lms[[4]] <- update(lms[[1]], fixed = as.formula(fixed4))
      }
      res.dt <- anova(base.lm, lms[[3]], lms[[4]]) %>%
        as.data.table() %>%
        set_names(toTitleCase(str_replace_all(names(.), "-", ".")))
      if (all(res.dt$P.value[2:3] > 0.05)) {
        mod <- base.lm
      } else if (res.dt$P.value[3] <= 0.05) {
        mod <- lms[[4]]
      } else {
        mod <- lms[[3]]
      }
      predictions.dt <- lapply(bap.concs, function(conc) {
        pred.dt <- data.table(
          Alpha = seq(min(sample.dt[[alpha]]), max(sample.dt[[alpha]]), length = 1000),
          BaP_uM = conc
        )
        names(pred.dt)[1] <- alpha
        pred.dt[, (response) := predict(mod, newdata = pred.dt, level = 0)]
      }) %>% rbindlist()
      saveRDS(
        predictions.dt,
        file = str_replace(alpha.lpr.lmes.dt.file, "dt\\.rds", paste0(alpha, "_preds_dt.rds"))
      )

      prettify.stat.table(mod) %>%
        .[
          , `:=`(
            Sig = ifelse(P.value <= 0.05, "*", ""),
            Alpha.metric = alpha,
            Cycle = cycle
          )
        ] %>%
        setcolorder(neworder = {ncol(.) - 1}:ncol(.)) %>%
        return()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

alpha.lpr.lmes.alt.dt.file <- file.path(dirs$save, "bapXbehav_alpha_lpr_lmes_alt_dt.rds")
alpha.lpr.lmes.alt.dt <- redo.if(script.var, alpha.lpr.lmes.alt.dt.file, {
  lapply(alpha.methods, function(alpha) {
    lapply(lpr.cycles, function(cycle) {
      cat(paste(alpha, cycle), sep = "\n")
      lms <- NULL
      response <- alpha
      predictors <- c("BaP_uM", str_subset(names(sample.dt), paste0("LPR.", cycle)))
      model.data <- sample.dt[complete.cases(sample.dt[, ..predictors])]
      random.effect <- str_subset(names(sample.dt), "SARL")
      .GlobalEnv$fixed1 <- paste(response, "~", predictors[1])
      .GlobalEnv$randEff <- paste("~1 |", random.effect)
      lms[[1]] <- nlme::lme(
        data = model.data,
        fixed = as.formula(fixed1),
        random = as.formula(randEff),
        na.action = na.omit,
        method = "ML",
        control = lmeControl(opt = 'optim')
      )
      .GlobalEnv$fixed2 <- paste0(fixed1, " + I(", predictors[1], "^2)")
      lms[[2]] <- update(lms[[1]], fixed = as.formula(fixed2))
      aov1 <- anova(lms[[1]], lms[[2]]) %>% as.data.table()
      if (aov1$`p-value`[2] <= 0.05) {
        base.lm <- lms[[2]]
        .GlobalEnv$fixed3 <- paste(fixed2, "+", predictors[2])
        lms[[3]] <- update(lms[[1]], fixed = as.formula(fixed3))
        .GlobalEnv$fixed4 <- paste(fixed2, "*", predictors[2])
        lms[[4]] <- update(lms[[1]], fixed = as.formula(fixed4))
      } else {
        base.lm <- lms[[1]]
        .GlobalEnv$fixed3 <- paste(fixed1, "+", predictors[2])
        lms[[3]] <- update(lms[[1]], fixed = as.formula(fixed3))
        .GlobalEnv$fixed4 <- paste(fixed1, "*", predictors[2])
        lms[[4]] <- update(lms[[1]], fixed = as.formula(fixed4))
      }
      res.dt <- anova(base.lm, lms[[3]], lms[[4]]) %>%
        as.data.table() %>%
        set_names(toTitleCase(str_replace_all(names(.), "-", ".")))
      if (all(res.dt$P.value[2:3] > 0.05)) {
        mod <- base.lm
      } else if (res.dt$P.value[3] <= 0.05) {
        mod <- lms[[4]]
      } else {
        mod <- lms[[3]]
      }
      predictions.dt <- lapply(bap.concs, function(conc) {
        pred.dt <- data.table(
          V1 = seq(
            min(sample.dt[[predictors[2]]], na.rm = T),
            max(sample.dt[[predictors[2]]], na.rm = T),
            length = 1000
          ),
          BaP_uM = conc
        )
        names(pred.dt)[1] <- predictors[2]
        pred.dt[, (response) := predict(mod, newdata = pred.dt, level = 0)]
      }) %>% rbindlist()
      saveRDS(
        predictions.dt,
        file = str_replace(alpha.lpr.lmes.alt.dt.file, "dt\\.rds", paste0(alpha, "_preds_dt.rds"))
      )

      tidy(mod) %>%
        as.data.table() %>%
        set_names(., toTitleCase(names(.))) %>%
        .[
          , `:=`(
            Sig = ifelse(P.value <= 0.05, "*", ""),
            Alpha.metric = alpha,
            Cycle = cycle
          )
        ] %>%
        setcolorder(neworder = {ncol(.) - 1}:ncol(.)) %>%
        return()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

### Beta-diversity by BaP exposure and movement AUCs
### EPR
beta.epr.res.file <- file.path(
  dirs$save,
  "bapXbehav_beta_epr_dbrdaData_and_anovaTables_list.rds"
)
beta.epr.res <- redo.if(script.var, beta.epr.res.file, {
  lapply(beta.methods, function(beta) {
    cat.n(beta)
    results.list <- NULL

    random.effect <- str_subset(names(sample.dt), "SARL")
    predictor <- str_subset(names(sample.dt), "EPR.")
    model.data <- sample.dt[complete.cases(sample.dt[, ..predictor])]
    dist.mat <- dist.list[[beta]] %>% dist_subset(model.data$Sample)
    frm0 <- paste("dist.mat ~ BaP_factor *", predictor)
    dbrda.base <- capscale(
      as.formula(frm0),
      data = model.data,
      na.action = na.exclude
    )
    frm.cond <- paste0(". ~ . + Condition(", random.effect, ")")
    dbrda.full <- update(dbrda.base, frm.cond)
    results.list$Full.ord.data <- get.biplot.data(ps = ps.processed, ord = dbrda.full)
    results.list$Full.permanova <- dbrda.permanova(dbrda.full, beta = beta)

    dbrda.step <- ordistep(dbrda.base, direction = "both")
    if (any(as.character(dbrda.step$call) == "dist.mat ~ 1")) {
      results.list$Step.ord.data <- NA
      results.list$Step.permanova <- NA
    } else {
      dbrda.step.cond <- update(dbrda.step, frm.cond)
      results.list$Step.ord.data <- get.biplot.data(ps = ps.processed, ord = dbrda.step.cond)
      results.list$Step.permanova <- dbrda.permanova(dbrda.step.cond, beta = beta)
    }
    return(results.list)
  })
})

beta.epr.full.aov.dt.file <- file.path(
  dirs$save,
  "bapXbehav_beta_epr_full_permanovas_dt.rds"
)
beta.epr.full.aov.dt <- redo.if(script.var, beta.epr.full.aov.dt.file, {
  lapply(beta.methods, function(beta) {
    beta.epr.res[[beta]]$Full.permanova %>% return()
  }) %>% rbindlist()
})

beta.epr.step.aov.dt.file <- file.path(
  dirs$save,
  "bapXbehav_beta_epr_step_permanovas_dt.rds"
)
beta.epr.step.aov.dt <- redo.if(script.var, beta.epr.step.aov.dt.file, {
  lapply(beta.methods, function(beta) {
    aov.dt <- beta.epr.res[[beta]]$Step.permanova
    { if (!is.na(aov.dt)) { return(aov.dt) } } %>% suppressWarnings()
  }) %>% rbindlist()
})

beta.epr.full.ord.dts.file <- file.path(
  dirs$save,
  "bapXbehav_beta_epr_full_dbrda_dt_list.rds"
)
beta.epr.full.ord.dts <- redo.if(script.var, beta.epr.full.ord.dts.file, {
  sections <- names(beta.epr.res[[1]]$Full.ord.data) %>% set_names()
  lapply(sections, function(section) {
    lapply(beta.methods, function(beta) {
      sec.data <- beta.epr.res[[beta]]$Full.ord.data[[section]]
      if ("data.table" %in% class(sec.data)) {
        sec.dt <- copy(sec.data)
        axes.cols <- which(grepl("CAP|MDS", names(sec.dt)))
        if (length(axes.cols) != 2) {
          stop("Number of columns that match CAP or MDS does not equal 2")
        }
        names(sec.dt)[axes.cols] <- c("Axis1", "Axis2")
        sec.dt[, Beta.method := beta]
      } else if (section == "axes.labs") {
        sec.dt <- data.table(
          Beta.method = beta,
          X.lab = sec.data[1],
          Y.lab = sec.data[2]
        )
      } else {
        sec.dt <- data.table(
          Beta.method = beta,
          Scale = sec.data
        )
      }
      return(sec.dt)
    }) %>% rbindlist()
  })
})


### LPR

beta.lpr.mantel.res.file <- file.path(
  dirs$save,
  "bapXbehav_beta_lpr_mantel_list.rds"
)

beta.lpr.mantel.res <- redo.if(script.var, beta.lpr.mantel.res.file, {
  lapply(beta.methods, function(beta) {
    cat.n(beta)
    lapply(lpr.cycles, function(cycle) {
      cat.n(cycle)
      lapply(setNames(bap.concs, paste0("uM", bap.concs)), function(conc) {
        response <- str_subset(names(sample.dt), paste0("LPR.", cycle))
        model.data <- sample.dt[complete.cases(sample.dt[, ..response])][BaP_uM == conc]
        beta.dist <- dist.list[[beta]] %>% dist_subset(model.data$Sample)
        auc.dist <- model.data[, ..response] %>%
          as.matrix() %>%
          set_rownames(model.data$Sample) %>%
          dist()
        res <- mantel(xdis = auc.dist, ydis = beta.dist, method = "spearman")
        data.table(
          Beta = beta,
          Cycle = cycle,
          BaP_uM = conc,
          Mantel.rho = res$statistic,
          Sig = res$signif
        ) %>% return()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
})



beta.lpr.dblm.res.file <- file.path(
  dirs$save,
  "bapXbehav_beta_lpr_dblm_list.rds"
)

beta.lpr.dblm.res <- redo.if(script.var, beta.lpr.dblm.res.file, {
  train.split <- 0.7
  lapply(beta.methods, function(beta) {
    cat.n(beta)
    lapply(lpr.cycles, function(cycle) {
      cat.n(cycle)
      lapply(setNames(bap.concs, paste0("uM", bap.concs)), function(conc) {
        response <- str_subset(names(sample.dt), paste0("LPR.", cycle))
        model.data <- sample.dt[complete.cases(sample.dt[, ..response])] %>%
          .[BaP_uM == conc]
        train.idx <- sample(1:nrow(model.data), size = train.split * nrow(model.data), replace = FALSE)
        train.data <- model.data[sort(train.idx)]
        test.data <- model.data[-sort(train.idx)]
        beta.D2 <- dist.list[[beta]] %>% dist_subset(model.data$Sample) %>% disttoD2()
        train.D2 <- beta.D2[sort(train.idx), sort(train.idx)]
        class(train.D2) <- "D2"
        test.D2 <- beta.D2[-sort(train.idx), sort(train.idx)]
        train.mod <- dbstats::dblm(D2 = train.D2, y = train.data[[response]])
        sum.train.mod <- summary(train.mod)
        test.pred <- predict(train.mod, newdata = test.D2, type.var = "D2")[, 1]
        rmse <- sqrt(mean({test.data[[response]] - test.pred}^2))
        data.table(
          Beta = beta,
          Cycle = cycle,
          BaP_uM = conc,
          RMSE = rmse,
          Gvar = sum.train.mod$gvar[, 1],
          Eff.rank = sum.train.mod$eff.rank,
          Rel.gvar = sum.train.mod$rel.gvar,
          Crit.val = sum.train.mod$crit.value
        ) %>% return()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
})
# beta.lpr.dblm.res[
#   , `:=`(
#     Cycle = factor(Cycle, levels = lpr.cycles),
#     Beta = factor(Beta, levels = beta.methods)
#   )
# ]
# ggplot(beta.lpr.dblm.res, aes(x = BaP_uM, y = RMSE)) +
#   geom_col(aes(fill = as.factor(BaP_uM))) +
#   bap.fill.scale +
#   scale_x_continuous(breaks = bap.concs) +
#   facet_grid(Beta ~ Cycle, scales = "free")

beta.lpr.res.file <- file.path(
  dirs$save,
  "bapXbehav_beta_lpr_dbrdaData_and_anovaTables_list.rds"
)

beta.lpr.res <- redo.if(script.var, beta.lpr.res.file, {
  lapply(beta.methods, function(beta) {
    cat.n(beta)
    lapply(lpr.cycles, function(cycle) {
      cat.n(cycle)
      results.list <- NULL

      random.effect <- str_subset(names(sample.dt), "SARL")
      predictor <- str_subset(names(sample.dt), paste0("LPR.", cycle))
      model.data <- sample.dt[complete.cases(sample.dt[, ..predictor])]
      dist.mat <- dist.list[[beta]] %>% dist_subset(model.data$Sample)
      frm0 <- paste("dist.mat ~ BaP_uM *", predictor)
      dbrda.base <- capscale(
        as.formula(frm0),
        data = model.data,
        na.action = na.exclude
      )
      adonis2(formula = dist.mat ~ BaP_uM * LPR.dark.auc, data = model.data, na.action = na.exclude)
      frm.cond <- paste0(". ~ . + Condition(", random.effect, ")")
      dbrda.full <- update(dbrda.base, frm.cond)
      results.list$Full.ord.data <- get.biplot.data(smpls = model.data, ord = dbrda.full)
      results.list$Full.permanova <- dbrda.permanova(dbrda.full, beta = beta, cycle = cycle)

      dbrda.step <- ordistep(dbrda.base, direction = "both")
      results.list$Step.uncond.ord.data <- get.biplot.data(smpls = model.data, ord = dbrda.step)
      if (any(as.character(dbrda.step$call) == "dist.mat ~ 1")) {
        results.list$Step.ord.data <- NA
        results.list$Step.permanova <- NA
      } else {
        dbrda.step.cond <- update(dbrda.step, frm.cond)
        results.list$Step.ord.data <- get.biplot.data(
          smpls = model.data,
          ord = dbrda.step.cond
        )
        results.list$Step.permanova <- dbrda.permanova(
          dbrda.step.cond,
          beta = beta,
          cycle = cycle
        )
      }
      results.list$PCoA <- pcoa(dist.mat)
      return(results.list)
    })
  })
})

beta.lpr.full.aov.dt.file <- file.path(
  dirs$save,
  "bapXbehav_beta_lpr_full_permanovas_dt.rds"
)
beta.lpr.full.aov.dt <- redo.if(script.var, beta.lpr.full.aov.dt.file, {
  lapply(beta.methods, function(beta) {
    lapply(lpr.cycles, function(cycle) {
      beta.lpr.res[[beta]][[cycle]]$Full.permanova %>% return()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

beta.lpr.step.aov.dt.file <- file.path(
  dirs$save,
  "bapXbehav_beta_lpr_step_permanovas_dt.rds"
)
beta.lpr.step.aov.dt <- redo.if(script.var, beta.lpr.step.aov.dt.file, {
  lapply(beta.methods, function(beta) {
    lapply(lpr.cycles, function(cycle) {
      aov.dt <- beta.lpr.res[[beta]][[cycle]]$Step.permanova
      { if (!is.na(aov.dt)) { return(aov.dt) } } %>% suppressWarnings()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

beta.lpr.full.ord.dts.file <- file.path(
  dirs$save,
  "bapXbehav_beta_lpr_full_dbrda_dt_list.rds"
)
beta.lpr.full.ord.dts <- redo.if(script.var, beta.lpr.full.ord.dts.file, {
  sections <- names(beta.lpr.res[[1]][[1]]$Full.ord.data) %>% set_names()
  lapply(sections, function(section) {
    lapply(beta.methods, function(beta) {
      lapply(lpr.cycles, function(cycle) {
        cat.n(paste(section, beta, cycle, sep = " - "))
        sec.data <- beta.lpr.res[[beta]][[cycle]]$Full.ord.data[[section]]
        if ("data.table" %in% class(sec.data)) {
          sec.dt <- copy(sec.data)
          axes.cols <- which(grepl("CAP|MDS", names(sec.dt)))
          if (length(axes.cols) != 2) {
            stop("Number of columns that match CAP or MDS does not equal 2")
          }
          names(sec.dt)[axes.cols] <- c("Axis1", "Axis2")
          sec.dt[, `:=`(Beta.method = beta, Cycle = cycle)]
          auc.col <- str_subset(names(sec.dt), paste0("LPR.", cycle))
          if (length(auc.col) == 1) { sec.dt$AUC <- sec.dt[[auc.col]] }
        } else if (section == "axes.labs") {
          sec.dt <- data.table(
            Beta.method = beta,
            Cycle = cycle,
            X.lab = sec.data[1],
            Y.lab = sec.data[2]
          )
        } else {
          sec.dt <- data.table(
            Beta.method = beta,
            Cycle = cycle,
            Scale = sec.data
          )
        }
        sec.dt[, Cycle := factor(Cycle, levels = c("light", "dark"))]
        return(sec.dt)
      }) %>% rbindlist()
    }) %>% rbindlist()
  })
})
beta.lpr.step.ord.dts.file <- file.path(
  dirs$save,
  "bapXbehav_beta_lpr_step_dbrda_dt_list.rds"
)
beta.lpr.step.ord.dts <- redo.if(script.var, beta.lpr.step.ord.dts.file, {
  good.beta <- min(which(!is.na(sapply(beta.lpr.res, function(x) x$Step.ord.data))))
  sections <- names(beta.lpr.res[[good.beta]][[1]]$Step.ord.data) %>% set_names()
  lapply(sections, function(section) {
    lapply(beta.methods, function(beta) {
      lapply(lpr.cycles, function(cycle) {
        beta.data <- beta.lpr.res[[beta]][[cycle]]$Step.ord.data
        if (!all(is.na(beta.data))) {
          sec.data <- beta.data[[section]]
          if ("data.table" %in% class(sec.data)) {
            sec.dt <- copy(sec.data)
            axes.cols <- which(grepl("CAP|MDS", names(sec.dt)))
            if (length(axes.cols) != 2) {
              stop("Number of columns that match CAP or MDS does not equal 2")
            }
            names(sec.dt)[axes.cols] <- c("Axis1", "Axis2")
            sec.dt[, `:=`(Beta.method = beta, Cycle = cycle)]
            auc.col <- str_subset(names(sec.dt), paste0("LPR.", cycle))
            if (length(auc.col) == 1) { sec.dt$AUC <- sec.dt[[auc.col]] }
          } else if (section == "axes.labs") {
            sec.dt <- data.table(
              Beta.method = beta,
              Cycle = cycle,
              X.lab = sec.data[1],
              Y.lab = sec.data[2]
            )
          } else {
            sec.dt <- data.table(
              Beta.method = beta,
              Cycle = cycle,
              Scale = sec.data
            )
          }
          sec.dt[, Cycle := factor(Cycle, levels = c("light", "dark"))]
          return(sec.dt)
        }
      }) %>% rbindlist()
    }) %>% rbindlist()
  })
})

