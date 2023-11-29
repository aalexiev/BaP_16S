# taxa_in_kit_controls.R

kit.ctrl.taxa <- colSums(asv.dts.list[[1]][sample.dt[is.na(SARL_plateID)]$Sample, -1]) %>%
  subset(. > 1) %>%
  names() %>%
  unique()
keep.cols <- c("Sample", kit.ctrl.taxa)
kit.ctrl.taxa.abunds <- asv.dts.list[["rarefied"]][sample.dt[is.na(SARL_plateID)]$Sample, ..keep.cols] %>%
  melt(id.vars = "Sample", variable.name = "Taxon", value.name = "Abundance")
write.table(
  dcast(... ~ Sample, data = taxa.dt[kit.ctrl.taxa.abunds, on = "Taxon"], value.var = "Abundance"),
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  file = file.path(saveDir, "bap_16S_kit_controls_taxa.csv")
)

plot.data <-  asv.dts.list[["rarefied"]][, ..keep.cols][sample.dt[, .(Sample, Dissect_result)], on = "Sample"] %>%
  melt(id.vars = c("Sample", "Dissect_result"), variable.name = "Taxon", value.name = "Abundance") %>%
  merge(taxa.dt[, .(Taxon, Family, Genus)], by = "Taxon")
plot.data[, Sample.type := ifelse(Dissect_result == "kit_control", "kit_control", "gut_sample")]

ggplot(plot.data, aes(x = Taxon, y = Abundance)) +
  geom_boxplot(aes(color = Sample.type)) +
  facet_wrap(~ Family, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave(
  filename = file.path(plotDir, "boxplots_kit_control_taxa_by_sampleType_and_Family.pdf"),
  width = 14,
  height = 40
  )

ggplot(plot.data, aes(x = Genus, y = Abundance)) +
  geom_boxplot(aes(color = Sample.type)) +
  facet_wrap(~ Family, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
ggsave(
  filename = file.path(plotDir, "boxplots_kit_control_Genera_by_sampleType_and_Family.pdf"),
  width = 14,
  height = 40
)
