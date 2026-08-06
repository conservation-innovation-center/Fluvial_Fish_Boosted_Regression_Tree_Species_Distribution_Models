#R code to develop species distribution models for fluvial fishes based on their native ranges using Boosted Regression Trees (BRTs)
#Developed in support of the USGS Aquatic GAP project
#Code developed by Hao Yu & Arthur Cooper, Research Associates, Department of Fisheries and Wildlife, Michigan State University

#Set working directory and remove all objects from the workspace
#setwd("K:/GIS/AFWA_BrookTrout/Data/Raw_Data")
setwd("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/AGAP_downloads_Jul2026")

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
#fish <- read.csv("USGS_data/fish_list.csv",header=T)
fish <- read.csv("BRT/fish_list_v2_0.csv",header=T)
# fish_name_itis<-fish[,c(6:8)] #subset table that includes fish species ITIS (Integrated Taxonomic Information System) code, common name, and scientific name
#fish_name_itis <- select(fish, itis_tsn, common_name, scientific_name)
#fish_name_itis <- read.csv('USGS_data/species_list_v2_0.csv')
fish_name_itis <- read.csv('PA/species_list_v2_0.csv')#CW note: why are we overwriting fish_name_itis?
#Separate unrestricted data (no sharing restriction/used in BRT model development) from restricted data (cannot be publicly shared/not used in BRT model development)
#CW note: where is our restricted data from? Is all the restricted & unrestricted data from AGAP?
fish_unrestricted <- fish #CW note: wait, was unrestricted data already isolated in this csv? What's the point of this?
# fish_unrestricted<-fish[is.na(fish$restricted),] #unrestricted fish data
# fish_restricted<-fish[!is.na(fish$restricted),] #restricted fish data

#Input predictor variables table for fluvial stream reaches 
# predictors_fluvial<-read.csv("Predictors.csv",header=T) 
#predictors_landscape<-read_parquet("nhdplusv2_agap_landscape_characterstics.parquet.gzip")
predictors_landscape<-read_parquet("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/nhdplusv2_agap_landscape_characterstics.parquet.gzip")
#dam_metrics <- read_csv_arrow("dam_fragmentation_metrics_nhdplusv21.csv")
dam_metrics <- read_csv_arrow("Dam_metrics/dam_fragmentation_metrics_nhdplusv21.csv")
predictors_fluvial <- predictors_landscape %>%
  merge(dam_metrics, by.x = "comid", by.y = "COMID")
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
#HUC8_native <- read.csv("fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_nas_ranges.csv")%>%
HUC8_native <- read.csv("BRT/fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_nas_ranges.csv")%>%
  mutate(HUC8_code = as.character(HUC8))%>%
  mutate(
    HUC8_code = replace_when(HUC8_code, nchar(HUC8_code) < 8 ~ paste0("0",HUC8_code)) # CW note: got an error here, even after specifically calling from dplyr. Probably because I'm not using Mike's environment
  )  

#CW: does this need to be subset to brook trout? Already appears to have been subset to native status

#CW note: need to update dplyr to at least 1.2.0 for the above code to work. Sigh.
# might need to fully update R to get update of dplyr to stick?
# wait... when I open R studio, it shows an updated version of dplyr (1.2.1). But then when I run the code, it reverts to 1.1.4? What?
# oh, it's because rlang is a dependency of dplyr and was out of date. Weird that it didn't automatically update when I updated dplyr.

##########################################################
#Import HUC8s for NHDPlusV2 and Fish List
#Combine HUC locations and predictors
##########################################################

# Need to join HUC8 info to NHDPlus stream reaches
#---- 
#HUC8_reaches <- readRDS('HUC8_streams.rds') 
HUC8_reaches <- readRDS("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/HUC8_streams.rds")
#HUC8_reaches2 = select(read.csv("BRT/fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_data.csv"), c("comid", "huc8"))
predictors_fluvial_HUC8<-merge(predictors_fluvial,HUC8_reaches,by.x="comid",by.y="comid") #CW: got an error here about needing a uniquely valid column. changing by.y from "COMID" to "comid" to match column name in HUC8 layer fixed it

#CW: testing whether first 8 digits of reach code is equal to HUC8 code
HUC8_reaches$reach8dig = substr(HUC8_reaches$reachcode, 1, 8)
HUC8_reaches$reachHUC_match <- HUC8_reaches$HUC8 == HUC8_reaches$reach8dig
unique(HUC8_reaches$reachHUC_match) # ok all the rows evaluated to TRUE, so maybe we're back to the drawing board?

#CW: let's try using NHD flowlines table instead
flowlines_path = "C:/Users/cweinstein/Documents/Projects/AFS_workshop/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07/NHDPlusNationalData/NHDPlusV21_National_Seamless_Flattened_Lower48.gdb/NHDSnapshot/NHDFlowline_Network"

#----

#predictors_fluvial_HUC8 <- predictors_fluvial
#Create fish list and associated HUC8 native ranges

#fish_list <- read.csv("USGS_data/species_list_v2_0.csv") # CW: why are we reading this again?
# fish_list<-read.csv("USGS_data/fish_list.csv",header=T)
# HUC8_native<-HUC8_native[HUC8_native$scientific_name %in% fish_list$scientific_name,]
# 
# HUC8_native<-HUC8_native[,1:2]

fish_name<-unique(fish_name_itis[fish_name_itis$scientific_name %in% HUC8_native$scientific_name,]) # CW: what is this?
fish_name<-unique(fish_name_itis[fish_name_itis$scientific_name %in% c('Salvelinus fontinalis'),])

############################################
#Develop Boosted Regression Tree (BRT) model
############################################

#Remove NAs if present
# SpeciesData<-na.omit(subset(fish_unrestricted,select=c(comid_v2,itis_tsn,sp_count)))

#Changes input species data table orientation from 'stacked' (species records in rows) to species data in individual columns
# Species_Matrix<-matrify(SpeciesData)
#SpeciesMatrix <- read.csv("USGS_data/agap_fish_dataset_v2_0.csv") # ***THESE ARE THE OBSERVATION POINTS
SpeciesMatrix <- read.csv("PA/agap_fish_dataset_v2_0.csv") # ***THESE ARE THE OBSERVATION POINTS
# TODO: look into filtering these
# TODO: double-check list of predictor variables is the same as AGAP's
# CW: SpeciesMatrix is already binary presence/absence. sum of brook trout column is 5090, which doesn't seem right?

#Create copy of species matrix that will later be converted from abundances to binary (0/1) presence/absence
# Species01<-Species_Matrix
# nfhp_Species01<-cbind(rownames(Species01),Species01)
# names(nfhp_Species01)[1]<-c("comid")
colnames(SpeciesMatrix) <- str_replace_all(colnames(SpeciesMatrix), pattern = "X", "")
# nfhp_Species01 <- SpeciesMatrix

#Create list of species and native HUC8s # CW: this part filters native huc8s to just brook trout
# TODO: try changing this so it filters to just stream reaches in AGAP's prediction outputs
AGAP_pred = read.csv("BRT/fluvial_fish_brt_predictions_v2_0/fluvial_fish_brt_predictions_v2_0.csv")
agap_pred_safo = read.csv("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/USGS_data/fluvial_fish_brt_predictions_v2_0_SAFO_ONLY.csv")
fish_listrange<-merge(fish_name,HUC8_native, by.x="scientific_name",by.y="scientific_name",all=FALSE)
fish_listrange<-fish_listrange[order(fish_listrange$itis_tsn.x),]

#save to separate csv for debugging
#write.csv(fish_listrange, "safo_nativerange.csv")
alicia_huc8s = c('01010006','01010007','01010008','01010009','01010010','01010011','01050004','01070006','01100005','01100006','02020001','02020002','02020003','02020004','02020005','02020006','02020007','02020008','02030101','02030102','02030103','02030201','02030202','02040101','02040102','02040104','02050101','02050102','02050103','02050104','02050105','04010101','04010102','04030103','04030104','04030108','04030109','04120101','04120102','04120103','04130001','04130002','04130003','04140101','04140102','04140201','04140202','04140203','04140301','04140302','04270101','04290001','04290002','04290003','04290004','04290005','04290006','04290007','04290008','04300101','04300102','04300103','04300104','04300105','04300106','04300107','04300108','04300109','04300201','04300202','04310001','05010001','05010002','05010004','07040002','09030002','01010001','01010002','01010003','01010004','01010005','01020001','01020002','01020003','01020004','01020005','01030001','01030002','01030003','01040001','01040002','01050001','01050002','01050003','01060001','01060002','01060003','01070001','01070002','01070003','01070004','01070005','01080101','01080102','01080103','01080104','01080105','01080106','01080107','01080201','01080202','01080203','01080204','01080205','01080206','01080207','01090001','01090002','01090003','01090004','01090005','01100001','01100002','01100003','01100004','01100005','01100006','01110000','02010001','02010002','02010003','02010004','02010005','02010006','02010007','02020001','02020002','02020003','02020004','02020005','02020006','02020007','02020008','02030101','02030102','02030103','02030104','02030105','02030201','02030202','02040101','02040102','02040103','02040104','02040105','02040106','02040203','02050101','02050102','02050103','02050104','02050105','02050106','02050107','02050201','02050202','02050203','02050204','02050205','02050206','02050301','02050302','02050303','02050304','02050305','02060003','02070001','02070002','02070003','02070004','02070005','02070006','02070007','02070008','02070009','02080103','02080201','02080202','02080203','02080204','03010101','03010103','03040101','03050101','03050105','03060101','03060102','04010201','04010202','04010301','04010302','04020101','04020102','04020103','04020104','04020105','04020201','04020202','04020203','04020300','04030101','04030102','04030104','04030105','04030106','04030107','04030108','04030110','04030111','04030112','04030201','04030202','04030203','04030204','04040002','04040003','04060103','04060104','04060105','04060106','04060107','04060200','04070001','04070002','04070003','04070004','04070005','04070006','04070007','04080300','04110003','04120101','04120102','04120103','04120104','04130001','04130002','04130003','04140101','04140102','04140201','04140202','04140203','04150101','04150102','04150200','04150301','04150302','04150303','04150304','04150305','04150306','04150307','05010001','05010002','05010003','05010005','05010006','05010007','05020001','05020004','05020006','05050001','05050002','05050003','05050005','05050007','06010101','06010102','06010103','06010105','06010106','06010107','06010108','06010201','06010202','06010203','06010204','06020002','06020003','07010206','07030001','07030002','07030003','07030005','07040001','07040003','07040004','07040005','07040006','07040007','07040008','07050001','07050002','07050003','07050004','07050005','07050006','07050007','07060001','07060002','07060003','07060004','07060005','07060006','07070001','07070002','07070003','07070004','07070005','07070006','07090001','07090003','07090004','07120006','09030001')
agap_huc8s = fish_listrange$HUC8_code
length(setdiff(alicia_huc8s, fish_listrange$HUC8_code))
length(setdiff(fish_listrange$HUC8_code, alicia_huc8s))

alicia_huc8s_df = data.frame(HUC8_code = unique(alicia_huc8s), CIC_NAS = 1)
agap_huc8s_df = data.frame(HUC8_code = fish_listrange$HUC8_code, AGAP_NR = 1)
AGAP_AFWA_NR_HUC8_comparison = merge(alicia_huc8s_df, agap_huc8s_df, by = "HUC8_code", all=TRUE)
#write.csv(AGAP_AFWA_NR_HUC8_comparison, "safo_nativerange_AGAP_CIC_comparison.csv")

#reading in HUC8s exported from joined NHD WBD layer
NHD_HUC8_NR = read.csv("K:/GIS/AFWA_BrookTrout/Data/Raw_Data/HUC8_NR_comparison_join_export.csv")
length(unique(NHD_HUC8_NR$HUC_8))
setdiff(agap_huc8s_df$HUC8_code, NHD_HUC8_NR$HUC8_code)
length(unique(setdiff(agap_huc8s_df$HUC8_code, NHD_HUC8_NR$HUC8_code)))

#exported WBD hU8 selected by AGAP HUCs:
WBD_HUC8_AGAP_exp = read.csv("K:/GIS/AFWA_BrookTrout/Data/Analysis/BRT_analysis/WBDHU8_AGAP_HUC8_select.csv") %>%
  mutate(HUC8_code = as.character(HUC8))%>%
  mutate(
    HUC8_code = replace_when(HUC8_code, nchar(HUC8_code) < 8 ~ paste0("0",HUC8_code)) # CW note: got an error here, even after specifically calling from dplyr. Probably because I'm not using Mike's environment
  )
head(WBD_HUC8_AGAP_exp)
WBD_HUC8_AGAP_list = WBD_HUC8_AGAP_exp$HUC8_code

setdiff(agap_huc8s, WBD_HUC8_AGAP_list)



#Begin BRT loop
# CW: this loop fails. Sigh.
# cW: this doesn't seem like it even needs to be a loop. We already subset this to just Brook Trout.
# CW: getting this out of a loop would certainly make it easier to figure out which line is throwing an error
i = 162003
#for(i in unique(fish_listrange$itis_tsn.x)){
  # Extract scientific and common name
  currentfish<-as.matrix(unique(subset(fish_listrange,itis_tsn.x==i,select=c(scientific_name,common_name))))
  scientific<-currentfish[1,1]
  scientific<-gsub(" ", "_", scientific) # CW: spaces already get replaced with underscores here. Make sure this line gets run before saving anything so we don't end up with spaces in filenames
  common<-currentfish[1,2]
  common<-gsub(" ", "_", common)

#Pull HUC8 range data for target species
HUC_rangeNATIVE<-subset(fish_listrange,itis_tsn.x==i,select=c(common_name,HUC8_code))
# HUC_rangeNATIVE<-filter(eel_brook_ranges, itis_tsn == i)%>%select(comid)

#write.csv(select(HUC_predictorsNATIVE,"HUC8_code", "reachcode"), "HUC_predictorsNATIVE_sub.csv")
  
#Use HUC8 range to create range-wide predictor variable table 
HUC_predictorsNATIVE<-merge(HUC_rangeNATIVE,predictors_fluvial_HUC8,by.x="HUC8_code", by.y="HUC8",all=FALSE)
#CW: there are only 265 unique HUC8s now. That doesn't explain why we ended up with more observations though?
#CW: also there are duplicate rows where one comid shows up multiple times - only 271231 unique comids, even though this table has 335377 rows. Weird.

#CW: side quest to double-check that first 8 digits of each reach code is still the same as the HUC8
HUC_predictorsNATIVE_sub = select(HUC_predictorsNATIVE, "HUC8_code", "comid", "reachcode")
HUC_predictorsNATIVE_sub$reach8dig = substr(HUC_predictorsNATIVE_sub$reachcode, 1, 8)
HUC_predictorsNATIVE_sub$reachHUC_match <- HUC_predictorsNATIVE_sub$HUC8 == HUC_predictorsNATIVE_sub$reach8dig
unique(HUC_predictorsNATIVE_sub$reachHUC_match) #ok these all still evaluate to TRUE

comid_text = paste(unique(HUC_predictorsNATIVE$comid), collapse=",")
writeLines(comid_text, "NATIVE_comid_txt.txt")


#HUC_predictorsNATIVE <- merge(HUC_rangeNATIVE, predictors_fluvial, by = "comid", all = FALSE) #CW: aha! this is the line that fails!
#Cw: ah. HUC_rangeNATIVE doesn't have comid column. so we need to relate HUC8 to comid in order to merge these 2 tables. Okie dokie.
#CW: But i thought this script makes that relate table somewhere? Check farther up and/or check what I did in my script.
#CW: ok I think I got HUCs from the REACHCODE column in flowlines. Why didn't I just use the HUC & comid columns in catchments?
#CW wait, but also comid & huc8 are already related to each other in agap_fish_dataset_v2_0.csv, which is already loaded as SpeciesMatrix. Am I missing something?

fish_species01<-subset(SpeciesMatrix,select=c("comid",i))

#Use HUC8 range to create predictor variable table for presence-absence locations
# CW: number of HUC8s decreased from 270 in HUC_rangeNATIVE to 254 in Species_VariableNATIVE. Why?
SpeciesVariableNATIVE<-merge(fish_species01,HUC_predictorsNATIVE,by="comid",all=FALSE)%>%
  #mutate(PA = fish_species01$162003)%>%
  mutate(NB_nlcd11_41_43 = NB_nlcd11b_41+ NB_nlcd11b_42+ NB_nlcd11b_43,
         N_nlcd11_90_95 = NB_nlcd11b_90+NB_nlcd11b_95, # TODO: double-check whether we should be using N_ or NB_ nlcd values
         N_nlcd11_21_24 = NB_nlcd11b_21+ NB_nlcd11b_22+ NB_nlcd11b_23+ NB_nlcd11b_24
         )
SpeciesVariableNATIVE$PA = SpeciesVariableNATIVE$`162003` #CW: what is the purpose of doing this? it gets removed 9 lines down
#   names(SpeciesVariableNATIVE)[2]<-c("Abundance")
# SpeciesVariableNATIVE$PA[SpeciesVariableNATIVE$Abundance>0]<-1
# SpeciesVariableNATIVE$PA[SpeciesVariableNATIVE$Abundance<1]<-0

stat.species<-c()

print(unique(SpeciesVariableNATIVE[4]))
#i_variables<-setdiff(c(1:ncol(SpeciesVariableNATIVE)),c(1:4,27))
exclude <- c("comid", i, "L_areasqkm","N_areasqkm", "PA")
i_variables<-setdiff(colnames(SpeciesVariableNATIVE),exclude) 
# CW: from the paper - 9 natural and 13 anthropogenic inputs. So 22 total, which matches the length we have here. Huh.
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
  "UDOR", #TODO: look for these.
  "UNDR",
  "DMD",
  "DM2D",
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
    
#CW: skip next several lines of code
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
#CW: pick back up here
BR <- readRDS("results/162003.rds")
predictors_fluvial <- mutate(predictors_fluvial, # CW: why are we mutating these again?
                             NB_nlcd11_41_43 = NB_nlcd11b_41+ NB_nlcd11b_42+ NB_nlcd11b_43,
       N_nlcd11_90_95 = NB_nlcd11b_90+NB_nlcd11b_95,
       N_nlcd11_21_24 = NB_nlcd11b_21+ NB_nlcd11b_22+ NB_nlcd11b_23+ NB_nlcd11b_24
)
# function that takes in single catchment, runs predictions with fitted model using variables as they are
# Also changes value of forest buffer variable by user-supplied amount, and re-runs predictions
# returns both sets of predictions
# CW: try adapting this for sensitivity analysis
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
predictor_species  = predictors_fluvial %>% select(all_of(include)) 
predict_native_region<-predict(BR,predictor_species,n.trees=BR$gbm.call$best.trees,type="response" )
#CW: compare outputs of this line to AGAP's outputs. They should match.

predict_native_region_prob<-as.data.frame(cbind(predictors_fluvial["comid"],predict_native_region))
names(predict_native_region_prob)[2]<-c("predict_prob")

#predict_native_region_prob$predict_PA<-ifelse(predict_native_region_prob$predict_prob>cutoff_TEST_BRT,1,0)

write.csv(predict_native_region_prob,paste("results/BRT_",scientific,"_",common,"_native_prediction.csv",sep=""))

save(BR,file=paste("results/BRT_",scientific,"_",common,"_.RData",sep=""))# Save to R data file
#BRT model loop

#END 
################

