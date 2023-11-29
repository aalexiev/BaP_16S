# 00_setup.R

maxCores <- 120
redos <- c(
  "data.prep",
  "bap.stats",
  "behavior.bap.bap",
  "lpr.auc.rfs",
  "spiec.easi",
  "cvcvzgf"
)
default.redo.state <- FALSE
set.to.false <- c()
set.to.true <- c()

### Directory and file names
dirs <- list(
  scripts = "Helper_scripts",
  save    = "Saved_objects",
  input   = "Input",
  plots   = "Plots",
  logs    = "Logs"
)
for (dir in dirs) {
  if (!dir.exists(dir)) { dir.create(dir) }
}
behaviorDir <- "Input/Behavior_results/GUI files/RemovedAffected"
behavior.files <- list(
  EPR = file.path(
    behaviorDir,
    "behavior_24h_FOR_GUI_HMAT_With_Removal_10-23-2020 112541.csv"
  ),
  LPR = file.path(
    behaviorDir,
    "behavior_5d_FOR_GUI_LPR_10-23-2020 112541.csv"
  )
)

source(file.path(dirs$scripts, "packages_sources.R"))
source(file.path(dirs$scripts, "functions.R"))

nCores <- ifelse(Sys.getenv("HOSTNAME") == "", 3, maxCores)
setDTthreads(nCores)
user.seed <- 33
rarefaction.min <- 10000

assign.redo(redos, state = default.redo.state)
if (!is.null(set.to.true)) { set.redo.true(set.to.true) }
if (!is.null(set.to.false)) { set.redo.false(set.to.false) }

smpl.col <- "Sample"
assays <- c("EPR", "LPR") %>% set_names()
lpr.cycles <- c("light", "dark") %>% set_names()
lpr.epochs <- paste0("E", 1:3) %>% set_names()
taxa.sets <- c("ASV_only", "aggregated") %>% set_names()
bap.concs <- c(0, 1, 5, 10)
conc.sets <- c("Concs.all", "Concs.0v10") %>% set_names()
alpha.methods <- c("Chao1", "Shannon", "Simpson", "Phylogenetic", "Observed") %>%
  set_names()
rf.types <- c("classification", "regression") %>% set_names()

ps.processed.file <- file.path(dirs$save, "processed_phyloseq.rds")
ps.processed <- redo.if("data.prep", ps.processed.file, {
  ps.prepped <- prep.phyloseq(
    input.ps.file = file.path(dirs$input, "phyloseq.rds"),
    qubit.data.file = file.path(dirs$input, "Qubit_data_edited.xlsx"),
    output.path = dirs$input
  )

  ps.behavior <-  prep.behavior.data(
    physeq = ps.prepped,
    behavior.assays = assays,
    behavior.files = behavior.files,
    metadata.file = file.path(dirs$input, "metadata.csv"),
    output.dir = dirs$input
  ) # AUC is calculated by TIMEPOINT here

  ps.decontam <- decontam.phyloseq(
    physeq = ps.behavior,
    conc.col = "ng.uL",
    neg.col = "Is_ctrl",
    decontam.method = "either",
    output.path = dirs$input
  )
  find.bad.batches(ps.decontam, metric = "Chao1", threshhold = 85)
  bad.batches <- c(1, 21, 22, 23) # from visual inspection of plots
  # produced by above function
  ps.good.batches <- subset_samples(
    ps.decontam,
    !(Extract_batch %in% bad.batches) & Is_ctrl == FALSE
  ) %>%
    prune_taxa(taxa_sums(.) > 0, .) %>%
    prune_samples(sample_sums(.) > 0, .)
  saveRDS(ps.good.batches, file = file.path(dirs$input, "good_batches_phyloseq.rds"))
  min.smpl.size <- min(
    sample_sums(ps.good.batches)[sample_sums(ps.good.batches) >= rarefaction.min]
  )
  saveRDS(min.smpl.size, file = file.path(dirs$save, "min_smpl_size.rds"))
  # Create rarefied and CLR-transformed phyloseqs
  rarefy_even_depth(
    physeq = ps.good.batches,
    sample.size = min.smpl.size,
    trimOTUs = TRUE,
    rngseed = user.seed
  )
})

### Generate beta-diversity distance matrices
dist.list.file <- file.path(dirs$save, "betaDiv_matrices_list.rds")
dist.list <- redo.if("data.prep", dist.list.file, {
  gen.dist.matrices(
    ps = ps.processed,
    methods = c("taxonomic", "phylogenetic"),
    cores = nCores
  )
})
beta.methods <- names(dist.list)[c(1:2, 4, 6, 3, 5)] %>%
  set_names()
beta.setNames <- function(x) { setNames(x, beta.methods) } # for use with `foreach`

### Estimate alpha-diversity metric scores
alpha.base.file <- file.path(dirs$save, "alphaDiv_metrics_dt.rds")
alpha.base <- redo.if("data.prep", alpha.base.file, {
  tax.dt <- estimate_richness(
    physeq = ps.processed,
    measures = alpha.methods[1:3]
  ) %>% as.data.table(keep.rownames = "Sample") %>%
    setkeyv("Sample")
  tax.dt[, se.chao1 := NULL]

  phy.dt <- pd(
    samp = otu.matrix(ps.processed),
    tree = phy_tree(ps.processed)
  ) %>%
    as.data.table(keep.rownames = "Sample") %>%
    setkeyv("Sample")
  names(phy.dt)[2:3] <- alpha.methods[4:5]
  tax.dt[phy.dt]
})
sample_data(ps.processed) <- merge(
  x = sample.data.table(ps.processed),
  y = alpha.base,
  by = "Sample"
) %>%
  .[, BaP_factor := as.factor(BaP_uM)] %>%
  as.data.frame() %>%
  set_rownames(., .[["Sample"]]) %>%
  select(-Sample) %>%
  sample_data()

list.redos() %>% print()
cat.n(paste("\nNumber of cores to be used:", nCores))
