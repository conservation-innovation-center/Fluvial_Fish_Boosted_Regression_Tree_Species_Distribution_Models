#R code to develop species distribution models for fluvial fishes based on their native ranges using Boosted Regression Trees (BRTs)
#Developed in support of the USGS Aquatic GAP project
#Code developed by Hao Yu & Arthur Cooper, Research Associates, Department of Fisheries and Wildlife, Michigan State University

#Set working directory and remove all objects from the workspace
setwd("K:/GIS/AFWA_BrookTrout/Data/Raw_Data")
rm(list = ls(all = TRUE))

#Load and attach necessary libraries
library(dismo) #'gbm.step' function to generate BRT models
library(labdsv) #'matrify' function to flip species data table orientation
library(dplyr)
library(stringr)
library(gbm)
library(arrow)
library(nhdplusTools)
####################################
#Import fish data and predictor data
####################################
# fish <- read.csv("Fish_data.csv", header = T) #input fish data table
# fish <- read.csv('fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_data.csv')
fish <- read.csv("USGS_data/fish_list.csv",header=T)
# fish_name_itis<-fish[,c(6:8)] #subset table that includes fish species ITIS (Integrated Taxonomic Information System) code, common name, and scientific name
fish_name_itis <- select(fish, itis_tsn, common_name, scientific_name)
fish_name_itis <- read.csv('USGS_data/species_list_v2_0.csv')
#Separate unrestricted data (no sharing restriction/used in BRT model development) from restricted data (cannot be publicly shared/not used in BRT model development)
fish_unrestricted <- fish
# fish_unrestricted<-fish[is.na(fish$restricted),] #unrestricted fish data
# fish_restricted<-fish[!is.na(fish$restricted),] #restricted fish data

#Input predictor variables table for fluvial stream reaches 
# predictors_fluvial<-read.csv("Predictors.csv",header=T) 
predictors_fluvial<-read_parquet("nhdplusv2_agap_landscape_characterstics.parquet.gzip")
# reachcodes <- get_nhdplus(comid = predictors_fluvial$comid, skip_geometry = TRUE, properties = 'reachcode')%>%
#   mutate(HUC8 = str_sub(reachcode, 1, 8))

#######################################
#Import USGS NAS native HUC8 range data
#######################################


# HUC8<-read.csv("HUC8_ranges.csv",header=T) #input species HUC8 range data table
# HUC8<-read.csv("fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_nas_ranges.csv")
# HUC8<-HUC8[,-4] #remove unneeded field
#Restrict HUC8 range data to 'native' HUC8s only (HUC8s with 'introduced' status representing non-native range are removed)
# HUC8_native<-HUC8[HUC8$origin_status=="Native",]
# rm(HUC8)

#HUC8 range data obtained from USGS Non-indigenous Aquatic Species (NAS) program
HUC8_native <- read.csv("fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_nas_ranges.csv")%>%
  mutate(HUC8_code = as.character(HUC8))%>%
  mutate(
    HUC8_code = replace_when(HUC8_code, nchar(HUC8_code) < 8 ~ paste0("0",HUC8_code))
  )  

##########################################################
#Import HUC8s for NHDPlusV2 and Fish List
#Combine HUC locations and predictors
##########################################################

# Need to join HUC8 info to NHDPlus stream reaches
#---- 
HUC8_reaches <- readRDS('HUC8_streams.rds')
predictors_fluvial_HUC8<-merge(predictors_fluvial,HUC8_reaches,by.x="comid",by.y="COMID")
#----

predictors_fluvial_HUC8 <- predictors_fluvial
#Create fish list and associated HUC8 native ranges

fish_list <- read.csv("USGS_data/species_list_v2_0.csv")
# fish_list<-read.csv("USGS_data/fish_list.csv",header=T)
# HUC8_native<-HUC8_native[HUC8_native$scientific_name %in% fish_list$scientific_name,]
# 
# HUC8_native<-HUC8_native[,1:2]

fish_name<-unique(fish_name_itis[fish_name_itis$scientific_name %in% HUC8_native$scientific_name,])
fish_name<-unique(fish_name_itis[fish_name_itis$scientific_name %in% c('Salvelinus fontinalis'),])

############################################
#Develop Boosted Regression Tree (BRT) model
############################################

#Remove NAs if present
# SpeciesData<-na.omit(subset(fish_unrestricted,select=c(comid_v2,itis_tsn,sp_count)))

#Changes input species data table orientation from 'stacked' (species records in rows) to species data in individual columns
# Species_Matrix<-matrify(SpeciesData)
SpeciesMatrix <- read.csv("USGS_data/agap_fish_dataset_v2_0.csv")

#Create copy of species matrix that will later be converted from abundances to binary (0/1) presence/absence
# Species01<-Species_Matrix
# nfhp_Species01<-cbind(rownames(Species01),Species01)
# names(nfhp_Species01)[1]<-c("comid")
colnames(SpeciesMatrix) <- str_replace_all(colnames(SpeciesMatrix), pattern = "X", "")
# nfhp_Species01 <- SpeciesMatrix

#Create list of species and native HUC8s
fish_listrange<-merge(fish_name,HUC8_native, by.x="scientific_name",by.y="scientific_name",all=FALSE)
fish_listrange<-fish_listrange[order(fish_listrange$itis_tsn.x),]

#Begin BRT loop
for(i in unique(fish_listrange$itis_tsn.x)){
  # Extract scientific and common name
  currentfish<-as.matrix(unique(subset(fish_listrange,itis_tsn.x==i,select=c(scientific_name,common_name))))
  scientific<-currentfish[1,1]
  scientific<-gsub(" ", "_", scientific)
  common<-currentfish[1,2]
  common<-gsub(" ", "_", common)

#Pull HUC8 range data for target species
HUC_rangeNATIVE<-subset(fish_listrange,itis_tsn.x==i,select=c(common_name,HUC8))
# HUC_rangeNATIVE<-filter(eel_brook_ranges, itis_tsn == i)%>%select(comid)
  
#Use HUC8 range to create range-wide predictor variable table 
# HUC_predictorsNATIVE<-merge(HUC_rangeNATIVE,predictors_fluvial_HUC8,by="HUC8",all=FALSE)
HUC_predictorsNATIVE <- merge(HUC_rangeNATIVE, predictors_fluvial, by = "comid", all = FALSE)

fish_species01<-subset(SpeciesMatrix,select=c("comid",i))

#Use HUC8 range to create predictor variable table for presence-absence locations
SpeciesVariableNATIVE<-merge(fish_species01,HUC_predictorsNATIVE,by="comid",all=FALSE)%>%
  mutate(PA = '162003')%>%
  mutate(NB_nlcd11_41_43 = NB_nlcd11b_41+ NB_nlcd11b_42+ NB_nlcd11b_43,
         N_nlcd11_90_95 = NB_nlcd11b_90+NB_nlcd11b_95,
         N_nlcd11_21_24 = NB_nlcd11b_21+ NB_nlcd11b_22+ NB_nlcd11b_23+ NB_nlcd11b_24
         )
#   names(SpeciesVariableNATIVE)[2]<-c("Abundance")
# SpeciesVariableNATIVE$PA[SpeciesVariableNATIVE$Abundance>0]<-1
# SpeciesVariableNATIVE$PA[SpeciesVariableNATIVE$Abundance<1]<-0

stat.species<-c()

print(unique(SpeciesVariableNATIVE[4]))
i_variables<-setdiff(c(1:ncol(SpeciesVariableNATIVE)),c(1:4,27))
exclude <- c("comid", i, "L_areasqkm","N_areasqkm", "PA")
i_variables<-setdiff(colnames(SpeciesVariableNATIVE),exclude) 
include <- c(
  "N_areasqkm",
  "N_bfi",
  "N_precip",
  "L_temp",
  "L_fl_slope",
  "L_maxelev",
  "NB_nlcd11_41_43",
  "N_nlcd11_11",
  "N_nlcd11_90_95",
  "N_nlcd11_21_24",
  "N_nlcd11_81",
  "N_nlcd11_82",
  "N_pop11den",
  "N_allepa_den",
  "N_allmine_den",
  "N_total_p_yield",
  "N_totww",
  # "UDOR",
  # "UNDR",
  # "DMD",
  # "DM2D",
  "N_rx_stlen_den")
#Set a seed number so that the process can be repeatable
set.seed(10)
#Check the number of species presences to determine the starting lr ("learning") rate in initial BRT model
np<-length(which(SpeciesVariableNATIVE[,"PA"]==1))#number of presences
nab<-length(which(SpeciesVariableNATIVE[,"PA"]==0))#number of absences
percentp<-np/nrow(SpeciesVariableNATIVE)
if (np < 100){lr<-0.01}else{lr<-0.05}

print(paste("no=",np,",percentp=",percentp,",lr=",lr,sep=" "))
BR<-NULL

#Build BRT model
#Reduce learning rate by half if the best tree model has > 1,000 trees
#Maximum number of trees is set at 10,000
count<-0
while(is.null(BR)){
BR<-gbm.step(data=SpeciesVariableNATIVE, gbm.x=include, gbm.y="PA",
             family = "bernoulli", tree.complexity = 5, learning.rate = lr, max.trees = 10000,
             plot.main=FALSE, keep.fold.models=TRUE, keep.fold.vector=TRUE, keep.fold.fit=TRUE)
       
      
if(!is.null(BR)){
BR_stat<-as.data.frame(cbind(BR$gbm.call$tree.complexity,BR$gbm.call$learning.rate,BR$gbm.call$best.trees))
        names(BR_stat)<-c("tree.complexity","learning.rate","best.n.trees")
print(BR_stat)
if(BR$gbm.call$best.trees<1000) {
BR<-NULL
}
}
      
count<-count+1 # to avoid endless if attempting with too few species presences
if(count>=10){
BR$gbm.call$tree.complexity<-9999
BR$gbm.call$learning.rate<-lr
BR$gbm.call$best.trees<-np
}
lr<-lr*0.5
}

saveRDS(BR, sprintf("results/%s.rds",i))
    
dev_exp<-1-(BR$cv.statistics$deviance.mean/BR$self.statistics$mean.null)#model deviance explained

#Gather BRT variable contributions output
varcont<-as.data.frame(BR$contributions)
varcont<-varcont[order(varcont$var),]
varcont<-as.data.frame(t(varcont))
varcont<-varcont[-1,]
varcont["Name"]<-scientific
varcont<-varcont[c(23,1:22)]
    
#Write results of variable contributions to a .csv table
if(i==min(unique(fish_listrange$itis_tsn.x))){
write.table(varcont,paste("results/BRT_VarContributions.csv",sep=""),sep=",",row.names=FALSE)
}else{
write.table(varcont,paste("results/BRT_VarContributions.csv",sep=""),sep=",",row.names=FALSE,
            col.names=FALSE,append=TRUE)
}

#Compile BRT model statistics
BR_stat<-as.data.frame(cbind(np,nab,percentp,BR$gbm.call$tree.complexity,BR$gbm.call$learning.rate,BR$gbm.call$best.trees,dev_exp,BR$self.statistics$discrimination,BR$cv.statistics$discrimination.mean))
names(BR_stat)<-c("presences", "absences", "prevalence", "tree.complexity","learning.rate","best.n.trees","dev_exp","train.auc","cv.auc")
stat<-cbind(BR_stat)
row.names(stat)<-scientific
stat.species<-rbind(stat.species,stat)

#Write BRT statistics to a new table
if(i==min(unique(fish_listrange$itis_tsn))){
write.table(stat.species,paste("results/BRT_Stats.csv",sep=""),sep=",",row.names=TRUE,col.names=NA)
}else{
write.table(stat.species,paste("results/BRT_Stats.csv",sep=""),sep=",",row.names=TRUE,
                  col.names=F,append=TRUE)
}

#Develop partial dependence plot of the top 12 predictors
gbm.plot(BR, n.plots=12, write.title = FALSE, common.scale = FALSE,plot.layout=c(3, 4))
  mtext(paste(gsub("_", " ", scientific)," (",gsub("_", " ", common),")",sep=""), outer = TRUE, line=-2, cex = 1.5)
  
savePlot(filename=paste("results/BRT_",scientific,"_",common,"_plots.pdf",sep=""),type=c("pdf"), device=dev.cur())

##############################
#Cross validation of BRT model
##############################    

pred_sum_BRT<-c()
predict_all<-c()
n.fold<-10
k<-0
for(k in 1:n.fold){
selector<-BR$fold.vector
i_fold_t<-which(selector!=k)
i_fold_v<-which(selector==k)
k.cv.tdata<-SpeciesVariableNATIVE[i_fold_t,]
k.cv.vdata<-SpeciesVariableNATIVE[i_fold_v,]

k.cv.fit<-BR$fold.models[[k]]

r_TEST<-predict(k.cv.fit, newdata=k.cv.vdata,n.trees=k.cv.fit$n.trees, type="response")
d_TEST <- as.data.frame(cbind(k.cv.vdata$PA, r_TEST))
dd_TEST<-as.data.frame(cbind(k.cv.vdata, r_TEST))
pred_sum_BRT<-rbind(pred_sum_BRT,d_TEST)
predict_all<-rbind(predict_all,dd_TEST)
}

Deviance_TEST_BRT<-calc.deviance(obs=pred_sum_BRT[,1], pred=pred_sum_BRT[,2], family="bernoulli",calc.mean=TRUE)
pres_TEST_BRT<-pred_sum_BRT[pred_sum_BRT[,1]==1, 2]
abs_TEST_BRT<-pred_sum_BRT[pred_sum_BRT[,1]==0, 2]

e_TEST_BRT <- evaluate(p=pres_TEST_BRT, a=abs_TEST_BRT)

#Plot and save AUC results
plot(e_TEST_BRT, 'ROC')

savePlot(filename=paste("results/BRT_",scientific,"_",common,"_AUC.tiff",sep=""),type=c("tiff"), device=dev.cur())

########################################################
#Develop presence/absence cutoffs and diagnostic metrics
########################################################

t_TEST_BRT <- threshold(e_TEST_BRT)
cutoff_TEST_BRT<-np/(np+nab)#the prevalence from "threshold" is modeled prevalence

colnames(pred_sum_BRT)<-c("PA","Predict")
pred_sum_BRT<-as.data.frame(pred_sum_BRT)

pred_sum_BRT$predict_PA<-ifelse(pred_sum_BRT$Predict>cutoff_TEST_BRT,1,0)
predict_all$predict_PA<-ifelse(predict_all$r_TEST>cutoff_TEST_BRT,1,0)
names(predict_all)[28]<-c("predict_prob")

conf_d_TEST_BRT<-table(pred_sum_BRT$predict_PA,pred_sum_BRT$PA)

sensitivity_TEST_BRT<-conf_d_TEST_BRT[2,2]/(conf_d_TEST_BRT[2,2]+conf_d_TEST_BRT[2,1])
specificity_TEST_BRT<-conf_d_TEST_BRT[1,1]/(conf_d_TEST_BRT[1,1]+conf_d_TEST_BRT[1,2])

TSS_TEST_BRT<-sensitivity_TEST_BRT+specificity_TEST_BRT-1

det_cv_fold_BRT<-as.data.frame(cbind(Deviance_TEST_BRT, np, nab,e_TEST_BRT@auc,e_TEST_BRT@cor,cutoff_TEST_BRT,sensitivity_TEST_BRT,specificity_TEST_BRT,TSS_TEST_BRT))

colnames(det_cv_fold_BRT)<-c("Deviance","np","na","auc","cor","threshold","sensitivity","specificity","TSS")
rownames(det_cv_fold_BRT)<-scientific
if(i==min(unique(fish_listrange$itis_tsn))){
write.table(det_cv_fold_BRT,paste("results/BRT_CV.csv",sep=""),sep=",",row.names=TRUE,col.names=NA)
}else{
write.table(det_cv_fold_BRT,paste("results/BRT_CV.csv",sep=""),sep=",",row.names=TRUE,
                  col.names=F,append=TRUE)
}

#Write output cross validation tables
write.csv(pred_sum_BRT,paste("results/BRT_CV_predict_",scientific,"_",common,".csv"))
write.csv(predict_all,paste("results/BRT_CV_predict_all_",scientific,"_",common,".csv"))


########################################################################
#Project model results to all fluvial stream reaches within native range
########################################################################
BR <- readRDS("results/162003.rds")
predictors_fluvial <- mutate(predictors_fluvial,
                             NB_nlcd11_41_43 = NB_nlcd11b_41+ NB_nlcd11b_42+ NB_nlcd11b_43,
       N_nlcd11_90_95 = NB_nlcd11b_90+NB_nlcd11b_95,
       N_nlcd11_21_24 = NB_nlcd11b_21+ NB_nlcd11b_22+ NB_nlcd11b_23+ NB_nlcd11b_24
)
predict_catchment <- function(id, value){
  input <- filter(predictors_fluvial, comid == id)%>%
    # select(i_variables)
    select(all_of(include))
  print(paste('current value', input$NB_nlcd11_41_43))
  preds_b <- predict(BR, input, n.trees=BR$gbm.call$best.trees,type="response")
  input$NB_nlcd11_41_43 <- value
  print(paste('new value', input$NB_nlcd11_41_43))
  preds_a <- predict(BR, input, n.trees=BR$gbm.call$best.trees,type="response")
  return(c(preds_b, preds_a))
}
predictor_species<-HUC_predictorsNATIVE[,c(4:25)]
predict_native_region<-predict(BR,predictor_species,n.trees=BR$gbm.call$best.trees,type="response" )

predict_native_region_prob<-as.data.frame(cbind(HUC_predictorsNATIVE[,c(1:3)],predict_native_region))
names(predict_native_region_prob)[4]<-c("predict_prob")

predict_native_region_prob$predict_PA<-ifelse(predict_native_region_prob$predict_prob>cutoff_TEST_BRT,1,0)

write.csv(predict_native_region_prob,paste("results/BRT_",scientific,"_",common,"_native_prediction.csv",sep=""))

save(BR,file=paste("results/BRT_",scientific,"_",common,"_.RData",sep=""))# Save to R data file
}#BRT model loop

#END 
################

