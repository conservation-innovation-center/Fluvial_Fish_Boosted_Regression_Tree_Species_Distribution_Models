## scripting up workflow for sensitivity analysis

#Load and attach necessary libraries
library(dismo) #'gbm.step' function to generate BRT models
library(labdsv) #'matrify' function to flip species data table orientation
library(dplyr)
library(stringr)
library(gbm)
library(arrow)
library(nhdplusTools)
library(ggplot2)

setwd("K:/GIS/AFWA_BrookTrout/Data/Raw_Data")

# 1. import model
BR <- readRDS("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/AGAP_BT_BRT.rds")
BR_CIC <- readRDS("results/162003.rds")

# check which inputs are important
summary(BR, plotit=FALSE) # out of the ones that we can easily control, most important is NB_nlcd11_41_43. That's comforting.
summary(BR_CIC, plotit=FALSE)

gbm.plot(BR)
gbm.plot(BR_CIC)
gbm.plot(BR_CIC, 6)

# check for interactions (idk, Claude said this was important)
int <- gbm.interactions(BR)
int$interactions # ...but what counts as a strong interaction?
int$rank.list

# generate predictor variables

include <- c(
  "N_areasqkm", # this won't change
  "N_bfi",
  "N_precip",
  "L_temp",
  "L_fl_slope",
  "L_maxelev",
  "NB_nlcd11_41_43", #*
  "N_nlcd11_11",
  "N_nlcd11_90_95",
  "N_nlcd11_21_24", #*
  "N_nlcd11_81",
  "N_nlcd11_82",
  "N_pop11den",
  "N_allepa_den",
  "N_allmine_den",
  "N_total_p_yield",
  "N_totww",
  "UDOR", #TODO: look for these.
  "UNDR",
  "DMD",
  "DM2D",
  "N_rx_stlen_den")
#Set a seed number so that the process can be repeatable
set.seed(10)

predictors_fluvial <- mutate(predictors_fluvial,
                             NB_nlcd11b_41_43 = NB_nlcd11b_41+ NB_nlcd11b_42+ NB_nlcd11b_43,
                             N_nlcd11_90_95 = NB_nlcd11b_90+NB_nlcd11b_95,
                             N_nlcd11_21_24 = NB_nlcd11b_21+ NB_nlcd11b_22+ NB_nlcd11b_23+ NB_nlcd11b_24,
                             N_nlcd11_11c = N_nlcd11_11,
                             totww_mgalc = N_totww,
                             DM2D_Fishtail = DM2D
)

include <- c(
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
  "N_rx_stlen_den")

# 2. For each stream reach, run model with incremental changes in input variables



# function that takes in single catchment, runs predictions with fitted model using variables as they are
# Also changes value of forest buffer variable by user-supplied amount, and re-runs predictions
# returns both sets of predictions
# CW: try adapting this for sensitivity analysis
# predict_catchment <- function(id, value){
#   input <- filter(predictors_fluvial, comid == id)%>%
#     # select(i_variables)
#     select(all_of(include))
#   print(paste('current value', input$NB_nlcd11_41_43))
#   preds_b <- predict(BR, input, n.trees=BR$gbm.call$best.trees,type="response")
#   input$NB_nlcd11_41_43 <- value
#   print(paste('new value', input$NB_nlcd11_41_43))
#   preds_a <- predict(BR, input, n.trees=BR$gbm.call$best.trees,type="response")
#   return(c(preds_b, preds_a))
# }

predictor_species  = predictors_fluvial %>% 
  select(include) #%>%
  #slice_head(n=10)

# starting from here - wrap in a function that inputs idx, outputs dataframe
# use lapply to apply function to huc8

#bind rows of output dfs
# do multiple covars at a time, have a column that has the covar name and one that has the covar value in addition to the other columns for comid, probabiliyt


# pick out one stream to test sensitivity analysis on
idx = 50
#var = "N_nlcd11_11"
test_reach = predictor_species[idx,]
current_cond = test_reach$NB_nlcd11b_41_43

# duplicate this stream 10x, varying NB_nlcd11_41_43 by 10% each time
test_reach = predictor_species[rep(idx,101),]
test_comid = predictors_fluvial$comid[idx]
test_reach$NB_nlcd11b_41_43 = seq(from = 0, to = 100, by = 1)
test_reach$RunID = sprintf("%d_%06d", test_comid, 1:nrow(test_reach))
test_reach

test_reach_input = test_reach %>% select(include)


predict_native_region<-predict(BR,test_reach_input,n.trees=BR$gbm.call$best.trees,type="response" )
#CW: compare outputs of this line to AGAP's outputs. They should match.

#predict_native_region_prob<-as.data.frame(cbind(predictors_fluvial["comid"],predict_native_region))
predict_native_region_prob <- as.data.frame(cbind(test_reach[c("RunID","NB_nlcd11b_41_43")],predict_native_region))
names(predict_native_region_prob)[3]<-c("predict_prob")
head(predict_native_region_prob)

# save probability predictions
#write.csv(predict_native_region_prob, sprintf("K:/GIS/AFWA_BrookTrout/Data/Analysis/Sensitivity_analysis/AGAP_model/%d_NB_nlcd11_41_43_SA.csv", test_comid))

# approximate loess function by fitting a 3rd degree polynomial equation
fit_fn = lm(predict_prob ~ poly(NB_nlcd11b_41_43, 3), data = predict_native_region_prob)
coeff = coef(fit_fn)
print(fit_fn)

# create a knockoff partial dependence plot (I think that's what I'm doing at least?)
#
ggplot(predict_native_region_prob, aes(x=NB_nlcd11b_41_43, y=predict_prob)) +
  geom_point() +
  labs(
    title = sprintf("Stream %d", test_comid),
    x = "Network buffer forest land cover (%)",
    y = "Probability of brook trout occurrence"
  ) +
  geom_smooth(method="lm", formula = y~poly(x,3), col="blue") +
  geom_smooth(method="loess", col="purple") + 
  geom_vline(xintercept = current_cond, col="orange")



