# 06_cvcvzgf_behavior.R

if (!("dirs" %in% ls())) {
  source("Analysis/00_setup.R")
}

script.var <- "cvcvzgf"
script.var.state <- get.redo.state(script.var)
if (is.null(script.var.state)) {
  assign.redo(script.var, state = default.redo.state)
}
microbiome.plate.treat <- data.table(
  Plate.ID = c(
    25484,
    25485,
    25486,
    25487,
    25488,
    25489,
    25490,
    25492,
    25493
  ),
  Microbiome = c(
    "CV",
    "CV",
    "CV",
    "CVZ",
    "CVZ",
    "CVZ",
    "GF",
    "GF",
    "GF"
  )
)
cvcvzgf.behaviorDir <- "Input/Behavior_results/CV_CVZ_GF_redo"
cvcvzgf.behavior.files <- list(
  EPR = file.path(cvcvzgf.behaviorDir, "EPR.csv"),
  LPR = file.path(cvcvzgf.behaviorDir, "LPR.csv")
)
treatments <- c(CV = "conventional", CVZ = "conventionalized", GF = "germ-free")

cvcvzgf.dt.file <- file.path(dirs$save, "cvcvzgf_behavior_dt.rds")
cvcvzgf.dt <- redo.if(script.var, cvcvzgf.dt.file, {
  prep.behavior.data.noPS(
    behavior.assays = assays,
    behavior.files = cvcvzgf.behavior.files,
    metadata.file = file.path(dirs$input, "cv_cvz_gf_metadata.csv"),
    output.dir = dirs$input
  )
}) %>%
  merge(microbiome.plate.treat, by.x = "Plate_ID", by.y = "Plate.ID")
# View(cvcvzgf.dt)
# table(cvcvzgf.dt$Plate_ID)
cvcvzgf.dt$Treatment <- sapply(
  cvcvzgf.dt$Plate_description,
  function(x) names(treatments)[which(treatments == x)]
) %>% factor(levels = names(treatments))
cvcvzgf.dt[
  , `:=`(
    BaP_factor = factor(BaP_uM),
    Plate_ID = factor(Plate_ID),
    Microbiome = factor(Microbiome, levels = c("CV", "CVZ", "GF"))
  )
]

### Movement by BaP exposure
### EPR
# epr.ccg.bap.stat.tbl.file <- file.path(dirs$save, "cvcvzgf_bap_epr_lms_dt.rds")
# epr.ccg.bap.stat.tbl <- redo.if(script.var, epr.ccg.bap.stat.tbl.file, {
#   response <- "EPR.auc"
#   predictors <- c("BaP_factor", "Treatment")
#   random.effect <- "Plate_ID"
#   frm.full <- paste(response, "~", paste(predictors, collapse = " * "))
#   random.eff.frm <- paste0("~ 1|", random.effect)
#   mod.full <- nlme::lme(
#     data = cvcvzgf.dt,
#     fixed = as.formula(frm.full),
#     random = as.formula(random.eff.frm),
#     na.action = na.omit,
#     method = "ML",
#     control = lmeControl(opt = 'optim')
#   )
#   stat.tbl <- prettify.stat.table(mod.full) %>%
#     .[, Sig := ifelse(P.value <= 0.05, "*", "")]
#   stat.tbl[
#     Term != "(Intercept)",
#     Term := str_remove_all(Term, "_factor|Treatment")
#   ]
# })

### LPR
# lpr.ccg.bap.kw.tbl.file <- file.path(dirs$save, "bap_lpr_kws_dt.rds")
# lpr.ccg.bap.kw.tbl <- redo.if(script.var, lpr.ccg.bap.kw.tbl.file, {
#   lapply(lpr.cycles, function(cycle) {
#     response.var <- str_subset(names(sample.dt), paste0("LPR.", cycle))
#     frm <- as.formula(paste(response.var, "~ BaP_factor"))
#     kw.dt <- kruskal.test(
#       frm,
#       data = sample.dt
#     ) %>% tidy()
#     stats.dt <- data.table(
#       Cycle = cycle,
#       ChiSq = round(kw.dt$statistic, 2),
#       DF = kw.dt$parameter,
#       P.value = round(kw.dt$p.value, 3),
#       Sig = ifelse(kw.dt$p.value <= 0.05, "*", "")
#     )
#   }) %>% rbindlist()
# })

# set.redo.true(script.var)
lpr.ccg.bap.lm.res.file <- file.path(dirs$save, "bap_ccg_lpr_lms_list.rds")
lpr.ccg.bap.lm.res <- redo.if(script.var, lpr.ccg.bap.lm.res.file, {
  lapply(lpr.cycles, function(cycle) {
    response <- str_subset(names(cvcvzgf.dt), paste0("LPR.", cycle))
    predictors <- c("BaP_uM", "Microbiome")

    ### 1st order
    frm0 <- paste(response, "~", paste(predictors, collapse = " * "))
    mod0 <- glm( # microbiome treatment and plate ID are conflated in this design, so we can't use an lme here
      data = cvcvzgf.dt,
      formula = as.formula(frm0),
      na.action = na.omit
    )
    # anova(mod0)

    stat.tbl0 <- tidy(mod0) %>%
      as.data.table() %>%
      set_names(., toTitleCase(names(.))) %>%
      .[
        , `:=`(
          Sig = ifelse(P.value <= 0.05, "*", ""),
          Model = "1st order",
          Cycle = cycle
        )
      ]
    stat.tbl0[
      Term != "(Intercept)",
      Term := str_remove_all(Term, "I\\(|\\)|Microbiome") %>%
        str_replace("_", " ") %>%
        str_replace("u", "µ")
    ]
    setcolorder(stat.tbl0, (ncol(stat.tbl0) - 1):ncol(stat.tbl0))

    pred.dt0 <- lapply(c("CV", "CVZ", "GF"), function(microb) {
      pred.dt <- data.table(
        BaP_uM = seq(0, 10, length = 1000),
        Microbiome = microb
      )
      pred.dt[, (response) := predict(mod0, newdata = pred.dt, level = 0)]
    }) %>% rbindlist()
    names(pred.dt0)[3] <- "AUC"
    pred.dt0[, `:=`(Model = "1st order", Cycle = cycle)]

    ### 2nd order
    frm1 <- str_replace(
      frm0,
      predictors[1],
      paste0(predictors[1], " + I(", predictors[1], "^2)")
    )
    mod1 <- glm(
      data = cvcvzgf.dt,
      formula = as.formula(frm1),
      na.action = na.omit
    )
    stat.tbl1 <- tidy(mod1) %>%
      as.data.table() %>%
      set_names(., toTitleCase(names(.))) %>%
      .[
        , `:=`(
          Sig = ifelse(P.value <= 0.05, "*", ""),
          Model = "2nd order",
          Cycle = cycle
        )
      ]
    stat.tbl1[
      Term != "(Intercept)",
      Term := str_remove_all(Term, "I\\(|\\)|Microbiome") %>%
        str_replace("_", " ") %>%
        str_replace("u", "µ")
    ]
    setcolorder(stat.tbl1, (ncol(stat.tbl1) - 1):ncol(stat.tbl1))

    pred.dt1 <- lapply(c("CV", "CVZ", "GF"), function(microb) {
      pred.dt <- data.table(
        BaP_uM = seq(0, 10, length = 1000),
        Microbiome = microb
      )
      pred.dt[, (response) := predict(mod1, newdata = pred.dt, level = 0)]
    }) %>% rbindlist()
    names(pred.dt1)[3] <- "AUC"
    pred.dt1[, `:=`(Model = "2nd order", Cycle = cycle)]

    list(
      GLM = rbind(stat.tbl0, stat.tbl1),
      PRD = rbind(pred.dt0, pred.dt1)
    ) %>%
      return()
  })
})

lpr.ccg.bap.lm.tbl <- lapply(lpr.cycles, function(cycle) {
  lpr.ccg.bap.lm.res[[cycle]]$GLM
}) %>% rbindlist()
lpr.ccg.bap.lm.pred <- lapply(lpr.cycles, function(cycle) {
  lpr.ccg.bap.lm.res[[cycle]]$PRD
}) %>% rbindlist()
