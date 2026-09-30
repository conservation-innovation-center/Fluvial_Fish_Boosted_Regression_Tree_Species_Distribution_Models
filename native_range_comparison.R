# comparing the too-many versions of native range files that we have


# read in each version
brt_coarse_path = "K:/GIS/AFWA_BrookTrout/Data/Raw_Data/AGAP_downloads_Jul2026/BRT/fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_coarse_ranges.csv"
brt_nas_path = "K:/GIS/AFWA_BrookTrout/Data/Raw_Data/AGAP_downloads_Jul2026/BRT/fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_nas_ranges.csv"
agap_coarse_path = "K:/GIS/AFWA_BrookTrout/Data/Raw_Data/AGAP_downloads_Jul2026/Coarse_range_maps/agap_coarse_fish_ranges_v2_0.csv"
nas_range_path = "K:/GIS/AFWA_BrookTrout/Data/Raw_Data/AGAP_downloads_Jul2026/NAS/Data/Salvelinus_fontinalis_range_20171013.csv"

brt_coarse <- read.csv(brt_coarse_path)
brt_nas <- read.csv(brt_nas_path)
agap_coarse <- read.csv(agap_coarse_path)
nas_range <- read.csv(nas_range_path)

head(brt_coarse)
head(brt_nas)
head(agap_coarse)
head(nas_range)

#ok well just looking at the column names here it seems like brt_nas and/or nas_range are probably what we want to use, but let's double-check

# subset to just brook trout (SAFO)

brt_nas_safo = subset(brt_nas, itis_tsn == "162003")

head(brt_nas_safo)

unique(brt_nas_safo$origin_status) #only Native origin status, so all HUC8s listed should be native range
length(unique(brt_nas_safo$HUC8)) #270...good lord that's a lot

# nas_range has already been subset to just safo
unique(nas_range$ORIGIN) # both native and introduced are in this one. subset to just native
nas_range_native = subset(nas_range, ORIGIN == "Native")
head(nas_range_native)
length(unique(nas_range_native$HUC8))
# ok this also has 270 HUCs. Promising that they match.

# maybe let's do a set difference just to see whether they're exactly the same. If so I think we can just use one of these.
setdiff(brt_nas_safo$HUC8, nas_range_native$HUC8)
# no difference! Hooray!

# add a leading zero to the HUC8 column, then save brt_nas_safo to a separate csv for later use
brt_nas_safo_toSave <- brt_nas_safo %>%
  mutate(HUC8_code = as.character(HUC8))%>%
  mutate(
    HUC8_code = replace_when(HUC8_code, nchar(HUC8_code) < 8 ~ paste0("0",HUC8_code))
  )

write.csv(brt_nas_safo_toSave, "K:/GIS/AFWA_BrookTrout/Data/Raw_Data/AGAP_downloads_Jul2026/BRT/fluvial_fish_brt_model_artifacts_v2_0/brt_model_inputs/brt_fish_nas_ranges_safo.csv")
