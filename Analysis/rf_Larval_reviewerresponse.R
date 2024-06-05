## get confusion matrices for models
library(caret)
library(randomForest)
library(dplyr)
library(ranger)
library(janitor)
library(pROC)

setwd("/Users/alexieva/Documents/Projects/Analysis/KeatonLarvalBaP/BaP_16S")

## Random forest models determined how well we could predict BaP exposure from taxon abundances
# load test and training data for each type of model below
bap_rf_data_splits_list <- readRDS("Saved_objects/bap_rf_data_splits_list.rds")

## ASV only
rf_BaP_classification_ASV_only_train <- readRDS("Saved_objects/rf_BaP_classification_ASV_only_train.rds")
# mtry = 53, splitrule = extratrees was best
preds <- predict(rf_BaP_classification_ASV_only_train, 
                 bap_rf_data_splits_list$ASV_only$Test)

exp_dat <- bap_rf_data_splits_list$ASV_only$Test %>%
  mutate(BaP_uM = case_when(BaP_uM == "0" ~ "uM0",
                            BaP_uM == "1" ~ "uM1",
                            BaP_uM == "5" ~ "uM5",
                            TRUE ~ "uM10"))

confusionMatrix(preds, as.factor(exp_dat$BaP_uM))

# check roc matches what's in text
pred_prob <- predict(rf_BaP_classification_ASV_only_train, 
                     bap_rf_data_splits_list$ASV_only$Test,
                     type = "prob")
roc_obj <- multiclass.roc(as.factor(exp_dat$BaP_uM), pred_prob)
auc <- auc(roc_obj)
# 0.5416

## aggregated (we not only included abundances of individual ASVs, but also their higher-level taxonomic assignments (e.g., genus, family, class, etc.))
rf_BaP_classification_aggregated_train <- readRDS("Saved_objects/rf_BaP_classification_aggregated_train.rds")
preds_agg <- predict(rf_BaP_classification_aggregated_train, 
                 bap_rf_data_splits_list$aggregated$Test)

exp_dat_agg <- bap_rf_data_splits_list$aggregated$Test %>%
  mutate(BaP_uM = case_when(BaP_uM == "0" ~ "uM0",
                            BaP_uM == "1" ~ "uM1",
                            BaP_uM == "5" ~ "uM5",
                            TRUE ~ "uM10"))

confusionMatrix(preds_agg, as.factor(exp_dat_agg$BaP_uM))
# check roc matches what's in text
pred_prob_agg <- predict(rf_BaP_classification_aggregated_train, 
                     bap_rf_data_splits_list$aggregated$Test,
                     type = "prob")
roc_obj <- multiclass.roc(as.factor(exp_dat_agg$BaP_uM), pred_prob_agg)
auc <- auc(roc_obj)
# 0.5424


## Random forest regression models to predict LPR AUC values from taxon abundances
## RMSEs for models
feats <- readRDS("Saved_objects/lprAUCs_rf_importantFeatures_dt.rds")
feats

# this is the test and train data itself in case helpful
data_split <- readRDS("Saved_objects/lprAUCs_rf_data_splits_list.rds")


## make a panel for Fig 3
fig3 <- readRDS("Saved_objects/bap_beta_dbrda_dt_list.rds") 
data_panC <- as.data.frame(fig3$sample.coords) %>%
  dplyr::filter(Beta.method == "Unweighted UniFrac") %>%
  dplyr::select(Sample, Axis1, BaP_factor)
pal <- c("#006837", "#fdae61", "#f46d43", "#a50026")

fig3C <- ggplot(data = data_panC, aes(x = Axis1,
                                      y = BaP_factor,
                                      fill = BaP_factor)) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 20),
        legend.position = "none") +
  labs(x = "CAP1", y = "BaP Dosage") +
  scale_fill_manual(values = pal)

# ggsave("/Users/alexieva/Desktop/fig3barplot.png", 
#        height = 5, 
#        width = 6, 
#        units = "in", 
#        dpi = 300)


## p-value adjustment
tablebeta1 <- readRDS("Saved_objects/bapXbehav_beta_lpr_full_permanovas_dt.rds")
tablebeta1$p.adj <- p.adjust(tablebeta1$P.value, method = "BH")

tablebeta2 <- readRDS("Saved_objects/bap_beta_permanovas_dt.rds")
tablebeta2$p.adj <- p.adjust(tablebeta2$P.value, method = "BH")

