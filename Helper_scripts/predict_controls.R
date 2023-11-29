require(stringr)

inDir <- "Input"
saveDir <- "Saved_objects"
ps.stem <- "phyloseq_prepped"

ps.prepped <- readRDS(file.path(inDir, paste0(ps.stem, ".rds")))

sample.dt <- sample.data.table(ps.prepped)
asv.dt0 <- otu.data.table(ps.prepped)
asv.cols <- str_subset(names(asv.dt0), "^ASV")
asv.dt <- asv.dt0[, c(.(Sample = Sample), lapply(.SD, function(x) x / max(x))), .SDcols = asv.cols]
# asv.mat <- otu.matrix(ps.prepped)
# asv.mat[asv.mat > 0] <- 1
# rowSums(asv.mat)
taxa.dt <- taxa.data.table(ps.prepped)

rf.dt <- merge(sample.dt[, .(Sample, Is_ctrl)], asv.dt, by = "Sample")
rf.dt[, Sample := NULL]
rf.dt[, Is_ctrl := factor(ifelse(Is_ctrl, "Control", "Sample"))]
rf.file <- file.path(saveDir, "randForest_predict_controls.rds")
rf <- redo.if("random.forests", rf.file, {
  # from Helper_scripts/random_forest_functions.R
  run.random.forest(
    mod.data = rf.dt,
    response.var = "Is_ctrl",
    rf.type = "classification",
    n.cores = nCores
  )
})

rf.impt.dt.file <- file.path(saveDir, "dt_predict_controls_importantFeatures.rds")
rf.impt.dt <- redo.if("random.forests", rf.impt.dt.file, {
  id.sig.important.features(
    rf.model = rf,
    response.var = "Is_ctrl",
    n.cores = nCores,
    rnd.seed = user.seed
  )
})

taxa.dt[Taxon %in% rf.impt.dt$Taxon]

plot.cols <- c("Is_ctrl", rf.impt.dt$Taxon)
plot.dt <- rf.dt[, ..plot.cols] %>%
  melt(id.vars = "Is_ctrl", variable.name = "Taxon", value.name = "Abund") %>%
  merge(taxa.dt, by = "Taxon")
plot.dt[, Label := paste(Genus, Taxon, sep = "|")]

saveRDS(plot.dt, file = file.path(saveDir, "dt_predict_controls_taxa.rds"))

plot.dt <- readRDS(file.path(saveDir, "dt_predict_controls_taxa.rds"))

ggplot(plot.dt, aes(x = Is_ctrl, y = Abund)) +
  geom_quasirandom() +
  # stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", color ="red") +
  facet_wrap(~ Label) +
  labs(x = "Control vs Gut Sample", y = "Taxon Fractional Abundance")

control.taxa <- rf.dt[
  Is_ctrl == "Control",
  c(lapply(.SD, sum)),
  .SDcols = str_subset(names(rf.dt), "^ASV")
  ] %>%
  t() %>%
  as.data.table(keep.rownames = "Taxon") %>%
  .[V1 != 0] %>%
  extract2("Taxon")
control.taxa.cols <- c("Is_ctrl", control.taxa)
control.taxa.dt <- rf.dt[, ..control.taxa.cols] %>%
  melt(id.vars = "Is_ctrl", variable.name = "Taxon", value.name = "Abund") %>%
  merge(taxa.dt, by = "Taxon")
control.taxa.dt[, Label := paste(Genus, Taxon, sep = "|")]
ggplot(control.taxa.dt, aes(x = Is_ctrl, y = Abund)) +
  geom_quasirandom() +
  # stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", color ="red") +
  facet_wrap(~ Label) +
  labs(x = "Control vs Gut Sample", y = "Taxon Fractional Abundance")
