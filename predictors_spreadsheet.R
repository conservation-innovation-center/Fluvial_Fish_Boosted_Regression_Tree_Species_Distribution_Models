# making mega-spreadsheet with predictor variables & current predicted percentages

agap_pred_safo = read.csv("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/USGS_data/fluvial_fish_brt_predictions_v2_0_SAFO_ONLY.csv")

predictors_landscape<-read_parquet("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/nhdplusv2_agap_landscape_characterstics.parquet.gzip")

dam_metrics <- read_csv_arrow("Dam_metrics/dam_fragmentation_metrics_nhdplusv21.csv")

predictors_fluvial <- predictors_landscape %>%
  merge(dam_metrics, by.x = "comid", by.y = "COMID")

predictors_fluvial <- mutate(predictors_fluvial,
                             NB_nlcd11b_41_43 = NB_nlcd11b_41+ NB_nlcd11b_42+ NB_nlcd11b_43,
                             N_nlcd11_90_95 = N_nlcd11_90+N_nlcd11_95, #CW: changed this so it's using N_ instead of NB_ columns
                             N_nlcd11_21_24 = N_nlcd11_21+ N_nlcd11_22+ N_nlcd11_23+ N_nlcd11_24, #CW: fixed this too
                             N_nlcd11_11c = N_nlcd11_11,
                             totww_mgalc = N_totww,
                             DM2D_Fishtail = DM2D)

predictors_predictions_safo = predictors_fluvial %>%
  merge(agap_pred_safo, by="comid")

head(predictors_predictions_safo)

nrow(predictors_fluvial)
nrow(agap_pred_safo)
nrow(predictors_predictions_safo)

include <- c(
  "comid",
  "N_areasqkm", # this won't change
  "N_bfi",
  "N_precip",
  "L_temp",
  "L_fl_slope",
  "L_maxelev",
  "NB_nlcd11b_41_43", #*
  "N_nlcd11_11c",
  "N_nlcd11_90_95",
  "N_nlcd11_21_24", #*
  "N_nlcd11_81",
  "N_nlcd11_82",
  "N_pop11den",
  "N_allepa_den",
  "N_allmine_den",
  "N_total_p_yield",
  "totww_mgalc", # is this the same as totww_mgalc? N_totww
  "UDOR", #TODO: look for these.
  "UNDR",
  "DMD",
  "DM2D_Fishtail",
  "N_rx_stlen_den",
  "predict_prob",
  "predict_pa"
  )

predictors_predictions_safo_sub = predictors_predictions_safo %>% 
  dplyr::select(all_of(include))

write.csv(predictors_predictions_safo_sub, "K:/GIS/AFWA_BrookTrout/Data/Analysis/Sensitivity_analysis/covars_probabilities_safo_sub.csv")
