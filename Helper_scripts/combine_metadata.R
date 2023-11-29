# combine_metadata.R
require(phyloseq)
require(phyloseqCompanion)
require(magrittr)

colNames <- c("Sample", "FileID", "Origin", "BaP_uM", "RunID")

bap16s.dt0 <- readRDS("Input/metadata_withFileIDs.rds")
pilot.dt0 <- readRDS("Input/Pilot_data/phyloseq_obj_dada1.12.1_2019-11-12.rds") %>%
  sample.data.table()
pilot.epr.dt <- read.csv("Input/Pilot_data/EPR_data.csv") %>%
  as.data.table() %>%
  .[, c(2, 4)]
pilot.epr.dt[, conc := conc / 10]
pilot.dt1 <- pilot.dt0[pilot.epr.dt, on = "Sample==well"]
methods.dt0 <- readRDS("Input/Methods_data/dt_sample_data.rds")


bap16s.dt1 <- bap16s.dt0[, .(Sample, Barcode, Origin, BaP_uM, Run)]
names(bap16s.dt1) <- colNames
bap16s.dt1[, RunID := paste0("BaP16S.", RunID)]

pilot.dt2 <- pilot.dt1[, .(Sample, BaP_uM = conc)]
pilot.dt2[, `:=`(FileID = Sample, Origin = Sample, RunID = "BaP16SPilot")]
setcolorder(pilot.dt2, colNames)

methods.dt1 <- methods.dt0[, .(Sample, Barcode, OriginID, BaP.Treat, RunID)]
names(methods.dt1) <- colNames
methods.dt1[
  , `:=`(
    BaP_uM = ifelse(BaP_uM == "control", 0, 10),
    RunID = paste0("MicrobMethods.", str_replace(RunID, "R", "Run"))
    )
  ]

combined.metadata.dt <- rbind(bap16s.dt1, pilot.dt2, methods.dt1)
nrow(combined.metadata.dt)
length(unique(combined.metadata.dt$Sample))
unique(combined.metadata.dt$RunID)
saveRDS(combined.metadata.dt, file = "Input/dt_combined_metadata.rds")

fastq.locations <- c(
  BaP16S.Run1 = "/dfs/Sharpton_Lab/keaton/CV_BaP_16S_FASTQs/Run1",
  BaP16S.Run2 = "/dfs/Sharpton_Lab/keaton/CV_BaP_16S_FASTQs/Run2",
  BaP16SPilot = "/dfs/Sharpton_Lab/keaton/BaP_16S_Pilot/FASTQs",
  MicrobMethods.Run1 = "/dfs/Sharpton_Lab/keaton/Microbiome_methods_FASTQs/Run1",
  MicrobMethods.Run2 = "/dfs/Sharpton_Lab/keaton/Microbiome_methods_FASTQs/Run2"
)
saveRDS(fastq.locations, file = "Input/vec_combined_runs_fastq_locations.rds")
