# comparing our Brook Trout Model against AGAP's model

library(dplyr)



our_predictions = read.csv("K:/GIS/AFWA_BrookTrout/Data/Analysis/Sensitivity_analysis/model_validation/CIC_validation_prediction.csv")
#our_var_contributions = read.csv("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/results/BRT_VarContributions.csv")

AGAP_pred = read.csv("C:/Users/cweinstein/Documents/Projects/AFS_workshop/fluvial_fish_brt_predictions_v2_0/fluvial_fish_brt_predictions_v2_0.csv")
#AGAP_pred_safo = subset(AGAP_pred, scientific_name == "Salvelinus fontinalis")
#AGAP_pred_safo2 = subset(AGAP_pred, itis_tsn == "162003")

# Let's save a copy of just the brook trout predictions so I don't have to keep re-generating this
#write.csv(AGAP_pred_safo, "K:/GIS/AFWA_BrookTrout/Data/Raw_Data/USGS_data/fluvial_fish_brt_predictions_v2_0_SAFO_ONLY.csv")
AGAP_pred_safo = read.csv("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/USGS_data/fluvial_fish_brt_predictions_v2_0_SAFO_ONLY.csv")


# cursory comparison between our outputs & AGAP outputs
nrow(AGAP_pred_safo) # 277105
nrow(our_predictions) # 277104

# do set difference on comids between two sets of outputs
length(unique(AGAP_pred_safo$comid))  # 277105
length(unique(our_predictions$comid)) # 277104
setdiff(AGAP_pred_safo$comid, our_predictions$comid) # comid 4352486 is in AGAP predictions and not ours
length(setdiff(our_predictions$comid, AGAP_pred_safo$comid)) # 0 comids are included in our predictions and not AGAP's

# our predictions have wayyyy more rows than subsetted AGAP predictions. That doesn't seem right? Subsetting by itis produced the same number of observations as filtering by species name though

# maybe let's try using filter() function to see if that's any better?
#AGAP_pred_safo3 = AGAP_pred %>%
#  filter(itis_tsn == "162003")

# no, that didn't make a difference...what the heck?
# let's do a quick summary of each to see if that helps...
summary(AGAP_pred_safo$predict_prob)
summary(our_predictions$predict_prob_CIC)

# some changes in our predicted probability distributions...hopefully not a big deal?

# that didn't really help...there are a lot more tiny probabilities in our outputs though
# ok whatever let's just do the join and go from there

# Join our outputs & AGAP's outputs by comid to facilitate comparison
# let's try an outer join
probs_joined = our_predictions %>%
  full_join(AGAP_pred_safo, by = "comid") %>%
  #inner_join(AGAP_pred_safo, by = "comid") %>%
  select(comid, predict_prob_CIC, predict_prob) #%>%
  #dplyr::rename(predict_prob.x = AFWA_predict_prob)

# change column names to something that actually makes sense
colnames(probs_joined) = c("comid", "predict_prob_CIC", "predict_prob_AGAP")

#probs_joined = merge(x = our_predictions, y = AGAP_pred_safo, by = "comid")
head(probs_joined)

# calculate RMSE
library("Metrics")
# remove NA - should just remove one row
probs_joined_clean <- na.omit(probs_joined)
rmse(probs_joined_clean$predict_prob_CIC, probs_joined_clean$predict_prob_AGAP) # 0.09571297...doesn't seem terrible?
# after adding dam fragmentation, rmse lowered to 0.08320491...better I guess but not that much better?
# NOTE: this doesn't work when doing a full join. Only when doing an inner join
# UPDATE after using Hao's model directly and whatever data cleaning we did, rmse is now 0.06844332. Better than before, but still seems high maybe?

# calculate standard deviation to compare against rmse
sd(probs_joined_clean$predict_prob_CIC) # 0.2831714
sd(probs_joined_clean$predict_prob_AGAP) # 0.2903501

# calculate deltas
probs_joined_clean$delta = probs_joined_clean$predict_prob_AGAP-probs_joined_clean$predict_prob_CIC
hist(probs_joined_clean$delta) # normally distributed, with most close to zero...so that seems good?
range(probs_joined_clean$delta)
hist(abs(probs_joined_clean$delta), breaks=40)

# check how many rows have large positive or negative delta values
nrow(probs_joined_clean) # 277104 to start with
nrow(subset(probs_joined_clean, abs(delta)>.5))

# 6715, or 2.4% have absolute delta > 0.2. That isn't terrible at least
# 72037, or 26% have absolute delta > 0.05. That is worse, but also .05 is pretty small?
# only 135, or .05%, have absolute delta > 0.5. So that's nice. Only a small proportion of them are terrible.
# would still be nice to know why those are off, though.

# save joined tables for further examination in Arc
write.csv(probs_joined_clean, "K:/GIS/AFWA_BrookTrout/Data/Analysis/Sensitivity_analysis/model_validation/joined_validation_predictions.csv")

hist(AGAP_pred_safo$predict_prob)
hist(our_predictions$predict_prob)
  
# join outputs to stream reaches, visualize, calculate RMSE
head(our_predictions)
head(AGAP_pred_safo)

# read in joined predictions after joining to flowlines to see if there's any correlation between delta and stream length
#flowline_pred = read.csv("K:/GIS/AFWA_BrookTrout/Data/Analysis/BRT_analysis/SA_Flowlines_join_table.csv")
flowline_pred = read.csv("C:/Users/cweinstein/Documents/Projects/AFWA_2026/Sensitivity_analysis/model_validation/NHDFlowline_Network_validation_join_v0.csv")

large_delta = subset(flowline_pred, delta > 0.5)
large_delta_predictors <- large_delta %>%
  merge(predictors_predictions_safo, by.x = "COMID", by.y = "comid")

library(ggplot2)
ggplot(flowline_pred, aes(x=DivDASqKM, y=abs(delta))) +
  geom_point()

library(tidyverse)

# maybe let's try making some more plots to see if anything obviously correlates with delta or abs(delta)
flowline_pred_long = flowline_pred %>%
  pivot_longer(
    cols = c(LENGTHKM, WBAREACOMI, StreamLeve, StreamOrde, Pathlength, ArbolateSu, AreaSqKM, TotDASqKM, DivDASqKM, TOTMA, MAXELEVRAW, MINELEVRAW, SLOPE, predict_prob_CIC, predict_prob_AGAP), # The other variables to compare against
    names_to = "variable", 
    values_to = "value"
  )

ggplot(flowline_pred_long, aes(x = value, y = abs(delta))) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_wrap(~ variable, scales = "free_x") +
  theme_minimal() +
  labs(
    y = "abs(delta)", 
    x = "Other Variables",
    title = "Relationship of delta absolute value against other features"
  )

# let's try joining our deltas to our covariates table. That might reveal if any of our interpretations of the inputs are problematic
HUC_predictorsNATIVE <- read.csv("K:/GIS/AFWA_BrookTrout/Data/Analysis/Sensitivity_analysis/covars_probabilities_safo_sub.csv")
covars_pred = probs_joined_clean %>%
  full_join(HUC_predictorsNATIVE, by = "comid") #%>%
  #inner_join(AGAP_pred_safo, by = "comid") %>%
  #select(comid, predict_prob_CIC, predict_prob)

# now let's try making scatter plots, focusing on the covariates I flagged when switching to Hao's model
ggplot(covars_pred, aes(x=UDOR, y=delta)) +
  geom_point()


# let's try making several plots at once actually, this is getting tedious
covars_pred_longer = covars_pred %>%
  pivot_longer(cols = -delta, names_to = "variable", values_to = "value")

ggplot(covars_pred_longer, aes(x = value, y = abs(delta))) +
  geom_point(color = "steelblue", alpha = 0.7) +
  
  # 3. Create a separate panel for each variable. 
  # 'scales = "free_x"' is crucial so each plot has its own X-axis range!
  facet_wrap(~ variable, scales = "free_x") +
  theme_minimal() +
  labs(title = "Scatter Plots of delta absolute value Against covariates")

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
