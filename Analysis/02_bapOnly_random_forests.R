# 02_bapOnly_random_forests.R

if (!("dirs" %in% ls())) {
  source("Analysis/00_setup.R")
}

script.var <- "bap.stats"
script.var.state <- get.redo.state(script.var)
if (is.null(script.var.state)) {
  assign.redo(script.var, state = default.redo.state)
}

# Generate data for predicting BaP exposure with random forests
bap.rf.data.list.file <- file.path(dirs$save, "bap_rf_data_splits_list.rds")
bap.rf.data.list <- redo.if(script.var, bap.rf.data.list.file, {
  generate.bap.rf.data(
    physeq = ps.processed,
    split.strata = "BaP_uM"
  )
})

bap.rf.list <- lapply(taxa.sets, function(set) {
  lapply(rf.types, function(type) {
    filename.parts <- c("rf_BaP", type, set, "train.rds")
    rf.data <- copy(bap.rf.data.list[[set]]$Train)
    if (type == "classification") {
      rf.data[, BaP_uM := factor(paste0("uM", BaP_uM))]
    }
    rf.file <- file.path(dirs$save, paste(filename.parts, collapse = "_"))
    redo.if("bap.stats", rf.file, {
      run.random.forest(
        mod.data = rf.data,
        response.var = "BaP_uM",
        rf.type = type,
        n.cores = nCores
      )
    })
  })
})

bap.rmses.dt.file <- file.path(dirs$save, "bap_rf_rmses_dt.rds")
bap.rmses.dt <- redo.if(script.var, bap.rmses.dt.file, {
  lapply(taxa.sets, function(set) {
    rf.train <- bap.rf.list[[set]]$regression
    test.data <- copy(bap.rf.data.list[[set]]$Test)
    preds <- predict(rf.train, test.data)
    return(
      data.table(
        Taxa.set = set,
        RMSE = sqrt(mean((preds - test.data$BaP_uM)^2))
      )
    )
  }) %>% rbindlist()
})

bap.rocAUCs.dt.file <- file.path(dirs$save, "bap_rf_roc_aucs_dt.rds")
bap.rocAUCs.dt <- redo.if(script.var, bap.rocAUCs.dt.file, {
  lapply(taxa.sets, function(set) {
    rf.train <- bap.rf.list[[set]]$classification
    test.data <- copy(bap.rf.data.list[[set]]$Test)
    test.data[, BaP_uM := factor(paste0("uM", BaP_uM))]
    preds <- predict(rf.train, test.data, type = "prob")
    roc <- multiclass.roc(
      test.data$BaP_uM,
      preds
    )
    return(
      data.table(
        Taxa.set = set,
        ROC.AUC = as.numeric(auc(roc))
      )
    )
  }) %>% rbindlist()
})

bap.reg.impt.dt.file <- file.path(
  dirs$save, "bap_rf_regress_importantFeatures_dt.rds"
)
bap.reg.impt.dt <- redo.if(script.var, bap.reg.impt.dt.file, {
  taxa.set <- bap.rmses.dt[order(RMSE)]$Taxa.set[1]
  id.sig.important.features(
    rf.model = bap.rf.list[[taxa.set]]$regression,
    response.var = "BaP_uM",
    n.cores = nCores,
    rnd.seed = user.seed
  )
})

bap.class.impt.dt.file <- file.path(
  dirs$save, "bap_rf_classify_importantFeatures_dt.rds"
)
bap.class.impt.dt <- redo.if(script.var, bap.class.impt.dt.file, {
  taxa.set <- bap.rocAUCs.dt[order(ROC.AUC, decreasing = T)]$Taxa.set[1]
  id.sig.important.features(
    rf.model = bap.rf.list[[taxa.set]]$classification,
    response.var = "BaP_uM",
    n.cores = nCores,
    rnd.seed = user.seed
  )
})


abund.prev.res.file <- file.path(
  dirs$save,
  "bapOnly_taxa_abundance_prevalence_list.rds"
)
sample.dt <- copy(sample.data.table(ps.processed))
abund.prev.res <- redo.if(script.var, abund.prev.res.file, {
  sncm.dt <- lapply(bap.concs, function(conc) {
    conc.ps <- prune_samples(sample.dt[BaP_uM == conc]$Sample, ps.processed) %>%
      prune_taxa(taxa_sums(.) > 0, .)
    conc.asv.mat <-  otu.matrix(conc.ps)
    conc.tax.mat <- as(tax_table(conc.ps), "matrix")
    conc.res <- sncm.fit(spp = conc.asv.mat, taxon = conc.tax.mat, stats = F) %>%
      as.data.table()
    conc.res[
      , `:=`(
        BaP_uM = conc,
        Status = ifelse(freq > pred.upr, "Pro", ifelse(freq < pred.lwr, "Con", "Neut"))
      )
    ]
    return(conc.res)
  }) %>%
    rbindlist() %>%
    setkey(Taxon)

  changes.dt <- lapply(unique(sncm.dt$Taxon), function(taxon) {
    taxon.dt <- sncm.dt[taxon]
    data.table(
      Taxon = taxon,
      Changes = sapply(bap.concs, function(conc) {
        stat <- taxon.dt[BaP_uM == conc]$Status
        ifelse(length(stat) == 0, "Abs", stat) %>% return()
      }) %>% paste(collapse = "-")
    ) %>% return()
  }) %>% rbindlist()
  list(Modeling = sncm.dt, Changes = changes.dt)
})

no.changes <- sapply(c("Neut", "Pro", "Con", "Abs"), function(stat) {
  paste(rep(stat, 4), collapse = "-") %>% return()
})
unboring.changes.dt <- copy(abund.prev.res$Changes)[!(Changes %in% no.changes)] %>%
  .[
    !{
      str_detect(Changes, "Abs") &
        str_detect(Changes, "Neut") &
        !str_detect(Changes, "Pro") &
        !str_detect(Changes, "Con")
    }
  ]

overdisp.samples.dt.file <- file.path(
  dirs$save,
  "bapOnly_overdispersion_byBaPbySample_dt.rds"
)
overdisp.samples.dt <- redo.if(script.var, overdisp.samples.dt.file, {
  setkey(alpha.base, Sample)
  ref.tree <- prune_samples(sample.dt[BaP_uM == 0]$Sample, ps.processed) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    phy_tree()
  # ref.tree <- phy_tree(ps.processed)
  exp.pd <- expected.pd(ref.tree)
  var.pd <- variance.pd(ref.tree)

  lapply(bap.concs, function(conc) {
    conc.smpls <- sample.dt[BaP_uM == conc]$Sample
    alpha.dt <- copy(alpha.base)[conc.smpls]
    data.table(
      BaP_uM = conc,
      Obs.taxa = alpha.dt$Observed,
      Obs.PD = alpha.dt$Phylogenetic,
      Exp.PD = exp.pd[alpha.dt$Observed, 2],
      Var.PD = var.pd[alpha.dt$Observed, 2]
    ) %>% return()
  }) %>% rbindlist()
})

partitions <- c("Con", "Neut", "Pro") %>% set_names()
partition.data <- copy(abund.prev.res$Modeling)
overdisp.partition.dt.file <- file.path(
  dirs$save,
  "bapOnly_overdispersion_byBaPbyPartition_dt.rds"
)
overdisp.partition.dt <- redo.if(script.var, overdisp.partition.dt.file, {
  setkey(alpha.base, Sample)
  ref.tree <- prune_samples(sample.dt[BaP_uM == 0]$Sample, ps.processed) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    phy_tree()
  # ref.tree <- phy_tree(ps.processed)
  exp.pd <- rbind(data.frame(n = 0, expected.pd = 0), expected.pd(ref.tree))

  lapply(bap.concs, function(conc) {
    conc.smpls <- sample.dt[BaP_uM == conc]$Sample
    lapply(partitions, function(part) {
      part.taxa <- partition.data[BaP_uM == conc & Status == part]$Taxon
      sub.ps <- prune_samples(conc.smpls, ps.processed) %>%
        prune_taxa(part.taxa, .)
      sub.pd <- picante::pd(samp = otu.matrix(sub.ps), tree = phy_tree(sub.ps))
      data.table(
        BaP_uM = conc,
        Partition = part,
        Obs.taxa = sub.pd$SR,
        Obs.PD = sub.pd$PD,
        Exp.PD = exp.pd[sub.pd$SR + 1, 2]
      ) %>% return()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

# library(microbiomeMarker)
#
# marker.res.file <- file.path(dirs$save, "bapOnly_marker_taxa.rds")
# marker.res <- redo.if(script.var, marker.res.file, {
#   ps.mm <- copy(ps.processed)
#   tax_table(ps.mm) <- tax_table(ps.mm)[, -str_which(colnames(tax_table(ps.mm)), "Taxon")]
#   mm.res <- microbiomeMarker::compare_DA(
#     ps.mm,
#     group = "BaP_uM",
#     methods = c("ancombc", "deseq2", "lefse")
#   )
# })
# smp.tbl <- sample_data(ps.mm)
# smp.tbl <- smp.tbl[, c("BaP_uM", "Origin")]
# asv.tbl <- otu_table(ps.mm)
# tax.tbl <- tax_table(ps.mm)
# phy.tree <- phy_tree(ps.mm)
# ps.example <- phyloseq(
#   smp.tbl,
#   otu_table(ps.mm),
#   tax_table(ps.mm),
#   phy_tree(ps.mm)
# )
# saveRDS(ps.example, file = "~/Dropbox/Share/example_ps.rds")


# ps.url <- "https://www.dropbox.com/s/2ckxg7djky2nbfv/example_ps.rds?dl=1"
# ps.file <- "example_ps.rds"
# download.file(ps.url, destfile = ps.file)
# ps <- readRDS(ps.file)
# res <- microbiomeMarker::compare_DA(
#   ps,
#   group = "BaP_uM",
#   methods = c("ancombc", "deseq2", "lefse")
# )
# identical(taxa_names(otu_table(ps)), taxa_names(tax_table(ps)))
# identical(taxa_names(otu_table(ps)), taxa_names(phy_tree(ps)))
