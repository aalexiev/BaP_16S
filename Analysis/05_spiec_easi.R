# 05_speicEasi.R

if (!("dirs" %in% ls())) {
  source("Analysis/00_setup.R")
}

script.var <- "spiec.easi"
script.var.state <- get.redo.state(script.var)
if (is.null(script.var.state)) {
  assign.redo(script.var, state = default.redo.state)
}

taxa.levels <- c("phylum", "class", "order", "family", "genus", "ASV")
levels.dt <- data.table(
  Level = factor(taxa.levels, levels = taxa.levels),
  Plural = c("phyla", "classes", "orders", "families", "genera", "ASVs")
)
levels.dt[, Name := toTitleCase(as.character(Level))]
setkey(levels.dt, Level)

ps.good.batches <- readRDS(file.path(dirs$input, "good_batches_phyloseq.rds"))

for (lvl in as.character(levels.dt$Level)) {
  cat.n(paste("###", lvl, "###"))
  if (lvl != "ASV") {
    ps.glom.file <- file.path(dirs$save, paste0("phyloseq_goodBatches_", lvl, "Glom.rds"))
    physeq <- redo.if(script.var, ps.glom.file, {
      ps.glom <- tax_glom(ps.good.batches, taxrank = levels.dt[lvl]$Name)
      assignment.dt <- taxa.data.table(ps.good.batches) %>% setkeyv("Taxon")
      taxa_names(ps.glom) <- assignment.dt[taxa_names(ps.glom)][[levels.dt[lvl]$Name]]
      ps.glom
    })
  } else {
    physeq <- ps.good.batches
  }
  se.bap.file <- file.path(
    dirs$save,
    paste0("spiecEasi_", levels.dt[lvl]$Plural, "_byBaP_MB_lmr1e2_nlambda20_results_list.rds")
  )
  se.bap <- redo.if(script.var, se.bap.file, lf.dir = file.path(dirs$save, "Large_files"), {
    res.list <- lapply(bap.concs, function(conc) {
      ps <- prune_samples(sample.data.table(physeq)[BaP_uM == conc]$Sample, physeq) %>%
        prune_taxa(taxa_sums(.) > 0, .)
      spiec.easi(
        ps,
        method = 'mb',
        lambda.min.ratio = 1e-2,
        nlambda = 20,
        pulsar.params = list(rep.num = 50, ncores = nCores))
    }) %>% return()
    names(res.list) <- paste0("uM", bap.concs)
    res.list
  })
}
