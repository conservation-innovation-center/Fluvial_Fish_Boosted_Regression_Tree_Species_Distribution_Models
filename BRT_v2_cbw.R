#Code to develop Boosted Regression Tree (BRT) models for fluvial species within their native ranges
#Developed in support of the USGS Aquatic GAP project
#Code developed by Hao Yu, Research Associate, Department of Fisheries and Wildlife, Michigan State University

# [hopefully very minor] edits made by C. Weinstein to work with what we already have
# Fingers crossed nothing weird happens.
# August 2026.

#Set working directory and remove previous files if present
#setwd('D:/Backup_20211012/Lab/Data_new/AGAP_2/')
setwd("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/AGAP_downloads_Jul2026")
rm(list = ls(all = TRUE))

#Import necessary libraries
library(labdsv)#'matrify' function to flip species data table orientation
library(dplyr)# 'distinct' function to remove replicates 
library(dismo)#'gbm.step' function to generate BRT models
library(PresenceAbsence)#'presence.absence.accuracy' function to calculate model evaluation metrics
library(arrow) # CW: used for reading in parquet file

####################################
#Import fish data and predictir data
####################################
#fish<-read.csv("fish_all_clean_0820.csv", header = T)#input fish data table
# CW: the '_clean' part of this filename has me concerned. guess we'll just use the fish data we have though.
fish <- read.csv("BRT/fish_list_v2_0.csv",header=T)

#Separate unrestricted data (no sharing restriction/used in BRT model development) from restricted data (cannot be publicly shared/not used in BRT model development)
#CW: fish data that we're using is downloaded from ScienceBase, presumably only includes unrestricted data. skip the filtering stuff.
#fish_unrestricted<-fish[is.na(fish$restricted),]#unrestricted fish data
#fish_restricted<-fish[!is.na(fish$restricted),]#restricted fish data
#fish_name_itis<-fish_unrestricted[,c(6:8)]#subset table that inclues fish species ITIS (Integrated Taxonomic Information System) code, common name, and scientific name
fish_name_itis <- fish[,c(1:3)] # CW: chose these columns based on the comment above. column ordering in our file is apparently different from what Hao's using here

#Input predictor variables table for fluvial stream reaches 
#predictors_fluvial<-read.csv("predictors_all_fluvial_0524_2022.csv",header=T)#predictor variables

# CW: replaced the above with loading in landscape & dam metrics, joining tables together.
predictors_landscape<-read_parquet("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/nhdplusv2_agap_landscape_characterstics.parquet.gzip")
#dam_metrics <- read_csv_arrow("dam_fragmentation_metrics_nhdplusv21.csv")
dam_metrics <- read_csv_arrow("Dam_metrics/dam_fragmentation_metrics_nhdplusv21.csv")
predictors_fluvial <- predictors_landscape %>%
  merge(dam_metrics, by.x = "comid", by.y = "COMID")

#########################################
#link HUC8 ranges and predictor variables
#########################################
#spatial_HUC8<-read.csv("nhdplusv2_comid_huc8_2022_spatial_filter_no_name.csv",header=T)#input species HUC8 range data table

#CW: hopefully this is the equivalent of whatever file Hao's reading in above
spatial_HUC8 <- read.csv("BRT/fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_nas_ranges.csv")%>%
  mutate(HUC8_code = as.character(HUC8))%>%
  mutate(
    HUC8_code = replace_when(HUC8_code, nchar(HUC8_code) < 8 ~ paste0("0",HUC8_code))
  )
#CW: crap, this does not have COMID which is needed for the merge below. sigh. It does have HUC8 though. I guess this is where we need to relate HUC to reach to comid?
# yes ok this is where we read in the table that Mike made. Need to redo this part I guess(?), except why is this step even needed?

#CW: start by reading in the NHD flowlines I guess. We just need a table that has both reachcode and comid in it. then we can get huc8 from reachcode
NHD_flowlines <- read.csv("C:/Users/cweinstein/Documents/Projects/AFWA_2026/AGAP_downloads_Jul2026/BRT/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07/NHDFlowline_Network_ExportTable.csv")

#CW: based on the filename Hao uses for spatial_HUC8, I'm going to guess that this is supposed to have been filtered already. Let's try doing a rough filtering based on what Jared sent us
# original number of rows: 2691339
# Filtering steps: remove FTYPE of "Coastline" and "Pipeline"

#TODO: actually, skip the filtering for now. if it turns out to be necessary, add it here

# isolate just the COMID and REACHCODE. Use REACHCODE to identify HUC8, relating COMID to HUC8
# when REACHCODE is < 14 digits long, add in leading zero that likely got truncated because HUC naming conventions are stupid
HUC8_comid_df <- NHD_flowlines %>%
  dplyr::select(COMID, REACHCODE) %>%
  mutate(REACHCODE_char = as.character(REACHCODE)) %>%
  mutate(
    REACHCODE_char = replace_when(REACHCODE_char, nchar(REACHCODE_char) < 14 ~ paste0("0",REACHCODE_char))
  ) %>%
  mutate(HUC8 = substr(REACHCODE_char, 1, 8)) %>%
  dplyr::select(COMID, HUC8) %>%
  left_join(spatial_HUC8, by = c("HUC8" = "HUC8_code"))


#predictors_fluvial_HUC8<-merge(predictors_fluvial,spatial_HUC8,by.x="comid",by.y="COMID")#merge HUC8 range and predictor variables
predictors_fluvial_HUC8 <- merge(predictors_fluvial, HUC8_comid_df, by.x="comid", by.y="COMID") #CW: I think this is what we're trying to achieve, but I have questions about why we're doing it this way

#yeah ok the above line failed because the resulting file is too big. That makes sense. Try subsetting to just safo first since that's the only species we care about?

names(predictors_fluvial_HUC8)[24]<-c("HUC8")



##########################
#Import native HUC8 ranges
##########################
new_species<-read.csv("round_2_NAS_species_BRT.csv", header = T)

new_species_HUC<-read.csv("round_2_NAS_HUC_BRT.csv", header = T)

new_species_HUC_new<-new_species_HUC[new_species_HUC$Scientific_name %in% new_species$Scientific_name,]

new_species_HUC_ITIS<-merge(new_species_HUC_new,new_species,by.x="Scientific_name",by.y="Scientific_name", all.x=T)

HUC8_1_native<-new_species_HUC_ITIS

fish_name_1<-unique(fish_name_itis[fish_name_itis$itis_species %in% HUC8_1_native$itis_species,])

############################################
#Develop Boosted Regression Tree (BRT) model
############################################
#Remove NAs if present
SpeciesData<-na.omit(subset(fish_unrestricted,select=c(comid_v2,itis_species,sp_count)))

#Changes input species data table orientation from 'stacked' (species records in rows) to species data in individual columns
Species_Matrix<-matrify(SpeciesData) 

#Create copy of species matrix that will later be converted from abundances to binary (0/1) presence/absence
Species01<-Species_Matrix
nfhp_Species01<-cbind(rownames(Species01),Species01)#Add comid
names(nfhp_Species01)[1]<-c("comid")

#Create list of species and native HUC8s
fish_listrange_1<-HUC8_1_native
fish_listrange_1<-fish_listrange_1[order(fish_listrange_1$itis_species),]
names(fish_listrange_1)<-c("scientific_name","HUC8","itis_species","common_name")

#Begin BRT loop
for(i in unique(fish_listrange_1$itis_species)){
# Extract scientific and common name
  currentfish<-as.matrix(unique(subset(fish_listrange_1,itis_species==i,select=c(scientific_name,common_name))))
  scientific<-currentfish[1,1]
  scientific<-gsub(" ", "_", scientific)
  common<-currentfish[1,2]
  common<-gsub(" ", "_", common)

#Pull HUC8 range data for target species
HUC_rangeNATIVE<-subset(fish_listrange_1,itis_species==i,select=c(common_name,HUC8))
  
#Use HUC8 range to create range-wide predictor variable table  
HUC_predictorsNATIVE<-merge(HUC_rangeNATIVE,predictors_fluvial_HUC8,by="HUC8",all=FALSE)
fish_species01<-subset(nfhp_Species01,select=c("comid",i))

#Use HUC8 range to create predictor variable table for presence-absence locations
SpeciesVariableNATIVE<-merge(fish_species01,HUC_predictorsNATIVE,by="comid",all=FALSE)
names(SpeciesVariableNATIVE)[2]<-c("Abundance")
SpeciesVariableNATIVE$PA[SpeciesVariableNATIVE$Abundance>0]<-1
SpeciesVariableNATIVE$PA[SpeciesVariableNATIVE$Abundance<1]<-0

print(unique(SpeciesVariableNATIVE[4]))
i_variables<-setdiff(c(1:ncol(SpeciesVariableNATIVE)),c(1:4,27))

#checking percentage of presence of species to determine the starting lr rate.
np<-length(which(SpeciesVariableNATIVE[,"PA"]==1))#number of presences
nab<-length(which(SpeciesVariableNATIVE[,"PA"]==0))#number of absences
percentp<-np/nrow(SpeciesVariableNATIVE)
if (np < 100){lr<-0.01}else{lr<-0.05}

print(paste("no=",np,",percentp=",percentp,",lr=",lr,sep=" "))
BR<-NULL

#Build BRT model
#Reduce learning rate by half if the best tree model has < 1000 trees
count<-0
while(is.null(BR)){
set.seed(10)#set random seed to 10
BR<-gbm.step(data=SpeciesVariableNATIVE, gbm.x=i_variables, gbm.y="PA",
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
      
count<-count+1 # to avoid endless while loop when only few presences of species
if(count>=10){
BR$gbm.call$tree.complexity<-9999
BR$gbm.call$learning.rate<-lr
BR$gbm.call$best.trees<-np
}
lr<-lr*0.5
}

dev_exp<-1-(BR$cv.statistics$deviance.mean/BR$self.statistics$mean.null)#model deviance explained

#Gather BRT variable contributions output
varcont<-as.data.frame(BR$contributions)
varcont<-varcont[order(varcont$var),]
varcont<-as.data.frame(t(varcont))
varcont<-varcont[-1,]
varcont["Name"]<-scientific
varcont<-varcont[c(23,1:22)]

#Write results of variable contributions to a .csv table    
if(i==min(unique(fish_listrange_1$itis_species))){
write.table(varcont,paste("BRT_coarse_range/BRT_VarContributions.csv",sep=""),sep=",",row.names=FALSE)
}else{
write.table(varcont,paste("BRT_coarse_range/BRT_VarContributions.csv",sep=""),sep=",",row.names=FALSE,
            col.names=FALSE,append=TRUE)
}

#Compile BRT model statistics
BR_stat<-as.data.frame(cbind(np,nab,percentp,BR$gbm.call$tree.complexity,BR$gbm.call$learning.rate,BR$gbm.call$best.trees,dev_exp,
                      BR$self.statistics$discrimination,BR$cv.statistics$discrimination.mean,
                      BR$cv.statistics$correlation.mean,BR$cv.statistics$correlation.se))
names(BR_stat)<-c("presences", "absences", "prevalence", "tree.complexity","learning.rate","best.n.trees","dev_exp","train.auc","cv.auc","cv.corr","cv.corr.se")
stat<-BR_stat
row.names(stat)<-scientific

#Write BRT statistics to a new table  
if(i==min(unique(fish_listrange_1$itis_species))){
write.table(stat,paste("BRT_non_coarse_range/BRT_Stats.csv",sep=""),sep=",",row.names=TRUE,col.names=NA)
}else{
write.table(stat,paste("BRT_non_coarse_range/BRT_Stats.csv",sep=""),sep=",",row.names=TRUE,
                  col.names=F,append=TRUE)
}

#Develop partial dependence plot of the top 12 predictors
gbm.plot(BR, n.plots=12, write.title = FALSE, common.scale = FALSE,plot.layout=c(3, 4))
mtext(paste(gsub("_", " ", scientific)), outer = TRUE, line=-2, cex = 1.5,font=3)

#Export the plot to the assigned folder
savePlot(filename=paste("BRT_non_coarse_range/BRT_",scientific,"_",common,"_plots.pdf",sep=""),type=c("pdf"), device=dev.cur())

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

#Deviance_TEST_BRT<-calc.deviance(obs=pred_sum_BRT[,1], pred=pred_sum_BRT[,2], family="bernoulli",calc.mean=TRUE)
pres_TEST_BRT<-pred_sum_BRT[pred_sum_BRT[,1]==1, 2]
abs_TEST_BRT<-pred_sum_BRT[pred_sum_BRT[,1]==0, 2]

e_TEST_BRT <- evaluate(p=pres_TEST_BRT, a=abs_TEST_BRT)

#Plot and save AUC results
plot(e_TEST_BRT, 'ROC')
savePlot(filename=paste("BRT_non_coarse_range/BRT_",scientific,"_",common,"_AUC.tiff",sep=""),type=c("tiff"), device=dev.cur())

########################################################
#Develop presence/absence cutoffs and diagnostic metrics
########################################################

t_TEST_BRT <- threshold(e_TEST_BRT)
cutoff_TEST_BRT<-t_TEST_BRT$spec_sens

colnames(pred_sum_BRT)<-c("PA","Predict")
pred_sum_BRT<-as.data.frame(pred_sum_BRT)

pred_sum_BRT$predict_PA<-ifelse(pred_sum_BRT$Predict>cutoff_TEST_BRT,1,0)
predict_all$predict_PA<-ifelse(predict_all$r_TEST>cutoff_TEST_BRT,1,0)
names(predict_all)[28]<-c("predict_prob")

DATA<-cbind(row.names(pred_sum_BRT),pred_sum_BRT[,c(1,2)])
names(DATA)[1]<-c("ID")

#Calculate sensitivity, specificity,TSS, Kappa, and PCC
PA_package<-presence.absence.accuracy(DATA,threshold=cutoff_TEST_BRT,st.dev=FALSE)
conf_d_TEST_BRT<-table(pred_sum_BRT$predict_PA,pred_sum_BRT$PA)
sensitivity_TEST_BRT<-PA_package$sensitivity
specificity_TEST_BRT<-PA_package$specificity
TSS_TEST_BRT<-sensitivity_TEST_BRT+specificity_TEST_BRT-1
Kappa_TEST<-PA_package$Kappa
PCC_TEST<-PA_package$PCC

det_cv_fold_BRT<-as.data.frame(cbind(np, nab,e_TEST_BRT@auc,e_TEST_BRT@cor,cutoff_TEST_BRT,sensitivity_TEST_BRT,specificity_TEST_BRT,TSS_TEST_BRT,Kappa_TEST,PCC_TEST))

colnames(det_cv_fold_BRT)<-c("np","na","auc","cor","threshold","sensitivity","specificity","TSS","Kappa","PCC")
rownames(det_cv_fold_BRT)<-scientific

#Write the diagnostic metrics in a CSV file
if(i==min(unique(fish_listrange_1$itis_species))){
write.table(det_cv_fold_BRT,paste("BRT_non_coarse_range/BRT_CV.csv",sep=""),sep=",",row.names=TRUE,col.names=NA)
}else{
write.table(det_cv_fold_BRT,paste("BRT_non_coarse_range/BRT_CV.csv",sep=""),sep=",",row.names=TRUE,
                  col.names=F,append=TRUE)
}

#Write output crossvalidation tables
write.csv(pred_sum_BRT,paste("BRT_non_coarse_range/BRT_CV_predict_",scientific,"_",common,".csv"),row.names=F)
write.csv(predict_all,paste("BRT_non_coarse_range/BRT_CV_predict_all_",scientific,"_",common,".csv"),row.names=F)

########################################################################
#Project model results to all fluvial stream reaches within native range
########################################################################
predictor_species<-HUC_predictorsNATIVE[,c(4:25)]
predict_native_region<-predict(BR,predictor_species,n.trees=BR$gbm.call$best.trees,type="response" )

predict_native_region_prob<-as.data.frame(cbind(HUC_predictorsNATIVE[,c(1:3)],predict_native_region))
names(predict_native_region_prob)[4]<-c("predict_prob")

predict_native_region_prob$predict_PA<-ifelse(predict_native_region_prob$predict_prob>cutoff_TEST_BRT,1,0)

write.csv(predict_native_region_prob,paste("BRT_non_coarse_range/BRT_",scientific,"_",common,"_native_prediction.csv",sep=""),row.names=F)

save.image(paste("BRT_non_coarse_range/BRT_",scientific,"_",common,".RData",sep=""))# Save to R data file
}#BRT model loop

#####
#END#
#####







