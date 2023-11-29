# 01_bapOnly_predict_behavior_alpha_beta.R

if (!("dirs" %in% ls())) {
  source("Analysis/00_setup.R")
}

script.var <- "bap.stats"
script.var.state <- get.redo.state(script.var)
if (is.null(script.var.state)) {
  assign.redo(script.var, state = default.redo.state)
}
sample.dt <- copy(sample.data.table(ps.processed))

### Movement by BaP exposure
### EPR
epr.bap.stat.tbl.file <- file.path(dirs$save, "bap_epr_lms_dt.rds")
epr.bap.stat.tbl <- redo.if(script.var, epr.bap.stat.tbl.file, {
  rec <- recipe(EPR.auc ~ BaP_factor, data = sample.dt)
  mod <- linear_reg()
  wflow <- workflow() %>%
    add_recipe(rec) %>%
    add_model(mod)
  stat.tbl <- fit(wflow, data = sample.dt) %>%
    prettify.stat.table() %>%
    .[, Sig := ifelse(P.value <= 0.05, "*", "")]
  stat.tbl[
    Term != "(Intercept)",
    Term := str_replace(Term, "_factor", " ") %>%
      paste0(., "µM")
  ]
})

### LPR
lpr.bap.kw.tbl.file <- file.path(dirs$save, "bap_lpr_kws_dt.rds")
lpr.bap.kw.tbl <- redo.if(script.var, lpr.bap.kw.tbl.file, {
  lapply(lpr.cycles, function(cycle) {
    response.var <- str_subset(names(sample.dt), paste0("LPR.", cycle))
    frm <- as.formula(paste(response.var, "~ BaP_factor"))
    kw.dt <- kruskal.test(
      frm,
      data = sample.dt
    ) %>% tidy()
    stats.dt <- data.table(
      Cycle = cycle,
      ChiSq = round(kw.dt$statistic, 2),
      DF = kw.dt$parameter,
      P.value = round(kw.dt$p.value, 3),
      Sig = ifelse(kw.dt$p.value <= 0.05, "*", "")
    )
  }) %>% rbindlist()
})

lpr.bap.lm.tbl.file <- file.path(dirs$save, "bap_lpr_lms_dt.rds")
lpr.bap.lm.tbl <- redo.if(script.var, lpr.bap.lm.tbl.file, {
  lapply(lpr.cycles, function(cycle) {
    response <- str_subset(names(sample.dt), paste0("LPR.", cycle))
    predictor <- "BaP_uM"
    random.effect <- "SARL_plateID"
    frm0 <- paste(response, "~", predictor)
    mod0 <- nlme::lme(
      data = sample.dt,
      fixed = as.formula(frm0),
      random = as.formula(paste0("~ ", predictor, "|", random.effect)),
      na.action = na.omit
    )
    # anova(mod0)
    stat.tbl0 <- prettify.stat.table(mod0) %>%
      .[
        , `:=`(
          Sig = ifelse(P.value <= 0.05, "*", ""),
          Model = "1st order",
          Cycle = cycle
        )
      ]
    stat.tbl0[
      Term != "(Intercept)" & Effect == "fixed",
      Term := str_remove_all(Term, "I\\(|\\)") %>%
        str_replace("_", " ") %>%
        str_replace("u", "µ")
    ]
    setcolorder(stat.tbl0, (ncol(stat.tbl0) - 1):ncol(stat.tbl0))
    frm1 <- paste0(frm0, " + ", "I(", predictor, "^2)")
    mod1 <- nlme::lme(
      data = sample.dt,
      fixed = as.formula(frm1),
      random = as.formula(paste0("~ ", predictor, "|", random.effect)),
      na.action = na.omit
    )
    stat.tbl1 <- prettify.stat.table(mod1) %>%
      .[
        , `:=`(
          Sig = ifelse(P.value <= 0.05, "*", ""),
          Model = "2nd order",
          Cycle = cycle
        )
      ]
    stat.tbl1[
      Term != "(Intercept)" & Effect == "fixed",
      Term := str_remove_all(Term, "I\\(|\\)") %>%
        str_replace("_", " ") %>%
        str_replace("u", "µ")
    ]
    setcolorder(stat.tbl1, (ncol(stat.tbl1) - 1):ncol(stat.tbl1))
    return(rbind(stat.tbl0, stat.tbl1))
  }) %>% rbindlist()
})


### Alpha-diversity by BaP Exposure
alpha.bap.res.file <- file.path(dirs$save, "bap_alpha_glms_dt.rds")
alpha.bap.res <- redo.if(script.var, alpha.bap.res.file, {
  lapply(alpha.methods, function(alpha) {
    predictor <- "BaP_uM"
    random.effect <- "SARL_plateID"
    .GlobalEnv$frm0 <- paste(alpha, "~", predictor)
    mod0 <- nlme::lme(
      data = sample.dt,
      fixed = as.formula(frm0),
      random = as.formula(paste0("~ ", predictor, "|", random.effect)),
      na.action = na.omit,
      method = "ML",
      control = lmeControl(opt = 'optim')
    )
    .GlobalEnv$frm1 <- paste0(frm0, " + I(", predictor, "^2)")
    mod1 <- nlme::lme(
      data = sample.dt,
      fixed = as.formula(frm1),
      random = as.formula(paste0("~ ", predictor, "|", random.effect)),
      na.action = na.omit,
      method = "ML",
      control = lmeControl(opt = 'optim')
    )
    mod.comp <- anova(mod0, mod1) %>% as.data.table()
    if (mod.comp[2, 9] <= 0.05) {
      best.model <- mod1
    } else {
      best.model <- mod0
    }
    stat.tbl <- prettify.stat.table(best.model) %>%
      .[
        , `:=`(
          Sig = ifelse(P.value <= 0.05, "*", ""),
          Alpha.metric = alpha
        )
      ]
    stat.tbl[
      Term != "(Intercept)" & Effect == "fixed",
      Term := str_remove_all(Term, "I\\(|\\)") %>%
        str_replace("_", " ") %>%
        str_replace("u", "µ")
    ]
    setcolorder(stat.tbl, ncol(stat.tbl))
    return(stat.tbl)
  }) %>% rbindlist()
})

### Beta-diversity by BaP exposure
beta.bap.res.file <- file.path(
  dirs$save,
  "bap_beta_dbrdaData_and_anovaTables_list.rds"
)
beta.bap.res <- redo.if(script.var, beta.bap.res.file, {
  lapply(beta.methods, function(beta) {
    cat(beta, sep = "\n")
    dist.mat <- dist.list[[beta]]
    predictor <- "BaP_uM"
    random.effect <- "SARL_plateID"
    frm0 <- paste0("dist.mat ~ ", predictor)
    dbrda0 <- capscale(as.formula(frm0), data = sample.dt)
    frm1 <- paste0(". ~ . + I(", predictor, "^2)")
    dbrda1 <- update(dbrda0, frm1)
    mod.comp <- anova(dbrda0, dbrda1)
    frm <- paste0(". ~ . + Condition(", random.effect, ")")
    if (mod.comp[2, 6] <= 0.05) {
      dbrda <- update(dbrda1, as.formula(frm))
    } else {
      dbrda <- update(dbrda0, as.formula(frm))
    }
    dbrda.data <- get.biplot.data(ps = ps.processed, ord = dbrda)
    aov <- anova(dbrda, by = "term") %>%
      prettify.stat.table() %>%
      .[, Beta.metric := beta] %>%
      setcolorder(ncol(.))
    return(list(Ord.data = dbrda.data, Permanova = aov))
  })
})

beta.bap.aov.dt.file <- file.path(dirs$save, "bap_beta_permanovas_dt.rds")
beta.bap.aov.dt <- redo.if(script.var, beta.bap.aov.dt.file, {
  lapply(beta.methods, function(beta) {
    beta.bap.res[[beta]]$Permanova %>% return()
  }) %>% rbindlist()
})
beta.bap.aov.dt[, Sig := ifelse(P.value <= 0.05, "*", "")]

beta.bap.ord.dts.file <- file.path(dirs$save, "bap_beta_dbrda_dt_list.rds")
beta.bap.ord.dts <- redo.if(script.var, beta.bap.ord.dts.file, {
  sections <- names(beta.bap.res[[1]]$Ord.data) %>% set_names()
  lapply(sections, function(section) {
    lapply(beta.methods, function(beta) {
      sec.data <- beta.bap.res[[beta]]$Ord.data[[section]]
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

### BETADISPER ANALYSIS

beta.disp.dt.file <- file.path(dirs$save, "bapOnly_beta_dispersion.rds")
beta.disp.dt <- redo.if(script.var, beta.disp.dt.file, {
  lapply(beta.methods, function(beta) {
    bd.res <- vegan::betadisper(d = dist.list[[beta]], group = sample.dt$BaP_uM)
    bd.aov <- anova(bd.res) %>% tidy()
    bd.dt <- data.table(
      Sample = names(bd.res$distances),
      BaP_uM = bd.res$group,
      Dist.to.centr = bd.res$distances,
      Beta.metric = beta,
      Sig = bd.aov$p.value[1] <= 0.05
    )
  }) %>% rbindlist()
})

## NO SIGNIFICANT DIFFERENCES IN BETA-DISPERSION FOR ANY BETA METRIC
