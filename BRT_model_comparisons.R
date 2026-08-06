# comparing our Brook Trout Model against AGAP's model

library(dplyr)

our_predictions = read.csv("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/results/BRT_Salvelinus_fontinalis_brook_trout_native_prediction.csv")
our_var_contributions = read.csv("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/results/BRT_VarContributions.csv")

AGAP_pred = read.csv("C:/Users/cweinstein/Documents/Projects/AFS_workshop/fluvial_fish_brt_predictions_v2_0/fluvial_fish_brt_predictions_v2_0.csv")
#AGAP_pred_safo = subset(AGAP_pred, scientific_name == "Salvelinus fontinalis")
#AGAP_pred_safo2 = subset(AGAP_pred, itis_tsn == "162003")

# Let's save a copy of just the brook trout predictions so I don't have to keep re-generating this
#write.csv(AGAP_pred_safo, "K:/GIS/AFWA_BrookTrout/Data/Raw_Data/USGS_data/fluvial_fish_brt_predictions_v2_0_SAFO_ONLY.csv")
AGAP_pred_safo = read.csv("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/USGS_data/fluvial_fish_brt_predictions_v2_0_SAFO_ONLY.csv")

# cursory comparison between our outputs & AGAP outputs
nrow(AGAP_pred_safo)
nrow(our_predictions)

# do set difference on comids between two sets of outputs
length(unique(AGAP_pred_safo$comid)) # 277,105
length(unique(our_predictions$comid)) # 2,615,694
setdiff(AGAP_pred_safo$comid, our_predictions$comid) # comid 4352486 is in AGAP predictions and not ours
length(setdiff(our_predictions$comid, AGAP_pred_safo$comid)) # 2,338,590 comids are included in our predictions and not AGAP's

# our predictions have wayyyy more rows than subsetted AGAP predictions. That doesn't seem right? Subsetting by itis produced the same number of observations as filtering by species name though

# maybe let's try using filter() function to see if that's any better?
#AGAP_pred_safo3 = AGAP_pred %>%
#  filter(itis_tsn == "162003")

# no, that didn't make a difference...what the heck?
# let's do a quick summary of each to see if that helps...
summary(AGAP_pred_safo$predict_prob)
summary(our_predictions$predict_prob)

# that didn't really help...there are a lot more tiny probabilities in our outputs though
# ok whatever let's just do the join and go from there

# Join our outputs & AGAP's outputs by comid to facilitate comparison
# let's try an outer join
probs_joined = our_predictions %>%
  full_join(AGAP_pred_safo, by = "comid") %>%
  #inner_join(AGAP_pred_safo, by = "comid") %>%
  select(comid, predict_prob.x, predict_prob.y) #%>%
  #dplyr::rename(predict_prob.x = AFWA_predict_prob)

# change column names to something that actually makes sense
colnames(probs_joined) = c("comid", "AFWA_predict_prob", "AGAP_predict_prob")

#probs_joined = merge(x = our_predictions, y = AGAP_pred_safo, by = "comid")
head(probs_joined)

# calculate RMSE
library("Metrics")
rmse(probs_joined$AFWA_predict_prob, probs_joined$AGAP_predict_prob) # 0.09571297...doesn't seem terrible?
# after adding dam fragmentation, rmse lowered to 0.08320491...better I guess but not that much better?
# NOTE: this doesn't work when doing a full join. Only when doing an inner join

# calculate deltas
probs_joined$delta = probs_joined$AGAP_predict_prob-probs_joined$AFWA_predict_prob
hist(probs_joined$delta) # normally distributed, with most close to zero...so that seems good?

# save joined tables for further examination in Arc
#write.csv(probs_joined, "K:/GIS/AFWA_BrookTrout/Data/Analysis/BRT_analysis/full_joined_AGAP_AFWA_outputs.csv")


hist(AGAP_pred_safo$predict_prob)
hist(our_predictions$predict_prob)
  
# join outputs to stream reaches, visualize, calculate RMSE
head(our_predictions)
head(AGAP_pred_safo)

# read in joined predictions after joining to flowlines to see if there's any correlation between delta and stream length
flowline_pred = read.csv("K:/GIS/AFWA_BrookTrout/Data/Analysis/BRT_analysis/SA_Flowlines_join_table.csv")


library(ggplot2)
ggplot(flowline_pred, aes(x=Shape_Leng, y=delta)) +
  geom_point()


##### comparing our input stream reaches against AGAP's output stream reaches
#### our input stream reaches: HUC_predictorsNATIVE$comid (from AGAP_NAS_BRTs_mje.R)
#### AGAP's output stream reaches: AGAP_pred_safo$comid (from earlier in this script)

## make a new data frame combining the two, with indications of which dataset has which comid

CIC_input_comids_df = data.frame(comid = unique(HUC_predictorsNATIVE$comid), CIC_in = 1)
agap_output_comids_df = data.frame(comid = unique(AGAP_pred_safo$comid), AGAP_out = 1)
comid_in_out_compare = merge(CIC_input_comids_df, agap_output_comids_df, by = "comid", all=TRUE)
#write.csv(comid_in_out_compare, "comid_in_out_compare.csv")


# read in NHD flowlines that are included in our inputs and excluded from AGAP's outputs
excluded_reaches = read.csv("K:/GIS/AFWA_BrookTrout/Data/Analysis/BRT_analysis/NHDFlowline_Network_comid_in_out_compare_AGAP_excluded_reaches.csv")
head(excluded_reaches)
nrow(excluded_reaches)
