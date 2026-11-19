# Upstream network buffer area workflow

# 1. Generate 90m buffer for NHD flowlines
# 2. Calculate area of 90m buffer that occurs in catchment of each stream reach
# 3. Use NHD csvs to identify reaches that are upstream of each given stream reach
# 4. Sum buffer areas of each upstream stream reach.


import arcpy
import os
import time
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
arcpy.env.overwriteOutput=True
arcpy.env.outputCoordinateSystem = r"K:\GIS\AFWA_BrookTrout\Data\SA\All_NR_HUC8.shp"
arcpy.env.parallelProcessingFactor = "75%"
arcpy.env.compression = "LZ77"

# Easy function to make it more straightforward to print timestamps
def get_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S %m/%d/%Y", t)
    return current_time

print(get_time(),'Getting started!')
start_time = time.time()

generate_HUC8 = False
generate_90m_buff = False
generate_network_buffer_areas = False
summarize_network_forest_buffer = True

root_dir = r"C:\Users\cweinstein\Documents\Projects\AFWA_2026"
buffer_gdb = os.path.join(root_dir, "Sensitivity_analysis", "buffers_90m", "NHD_buffers_90m_v2.gdb")

NR_HUC8 = r"K:\GIS\AFWA_BrookTrout\Data\Raw_Data\AGAP_downloads_Jul2026\BRT\fluvial_fish_brt_model_artifacts_v2_0\brt_model_inputs\brt_fish_nas_ranges_safo.csv"

NHD_flowlines = os.path.join(root_dir, "AGAP_downloads_Jul2026", "BRT", "NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07", "NHDPlusNationalData", "NHDPlusV21_National_Seamless_Flattened_Lower48.gdb", "NHDSnapshot", "NHDFlowline_Network")
NHD_areas = os.path.join(root_dir, "AGAP_downloads_Jul2026", "BRT", "NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07", "NHDPlusNationalData", "NHDPlusV21_National_Seamless_Flattened_Lower48.gdb", "NHDSnapshot", "NHDArea")
NHD_waterbodies = os.path.join(root_dir, "AGAP_downloads_Jul2026", "BRT", "NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07", "NHDPlusNationalData", "NHDPlusV21_National_Seamless_Flattened_Lower48.gdb", "NHDSnapshot", "NHDWaterbody")

NHD_catchments = os.path.join(root_dir, "AGAP_downloads_Jul2026", "BRT", "NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07", "NHDPlusNationalData", "NHDPlusV21_National_Seamless_Flattened_Lower48.gdb", "NHDPlusCatchment", "Catchment_working")

NLCD_tif = os.path.join(root_dir, "Data", "Annual_NLCD_LndCov_2011_CU_C1V1", "Annual_NLCD_LndCov_2011_CU_C1V1.tif") # 2011 NLCD (year used for AGAP)
#NLCD_tif = os.path.join(root_dir, "Data", "Annual_NLCD_LndCov_2025_CU_C1V2", "Annual_NLCD_LndCov_2025_CU_C1V2.tif") # 2025 NLCD (most up to date layer at time of this analysis)

# initialize flowlines
NHD_lyr = arcpy.MakeFeatureLayer_management(NHD_flowlines, "NHD_lyr")
NHD_areas_lyr = arcpy.MakeFeatureLayer_management(NHD_areas, "NHD_areas_lyr")
NHD_waterbodies_lyr = arcpy.MakeFeatureLayer_management(NHD_waterbodies, "NHD_waterbodies_lyr")
Catchments_lyr = arcpy.MakeFeatureLayer_management(NHD_catchments, "Catchments_lyr")


# read in HUC8s, making sure to read HUC8_code as text so the leading zero doesn't get stripped off
# would love to know if there's a reason why so many HUC codes start with a zero. seems vastly inconvenient.
NR_HUC8_list = pd.read_csv(NR_HUC8, dtype={'HUC8_code': str})['HUC8_code']
for HUC8 in NR_HUC8_list:
    NR_flowlines = arcpy.SelectLayerByAttribute_management(NHD_lyr,"ADD_TO_SELECTION", f"HUC8 = '{HUC8}'")

def catchment_summary(Catchments_lyr, buff_out, out_fc, field_name, field_alias):
        print(get_time(), "summarizing area within each catchment")

        # select catchments using flowlines, since catchments table doesn't include HUC8
        NR_catchments = arcpy.management.SelectLayerByLocation(Catchments_lyr, "INTERSECT", NR_flowlines)

        # summarize 90m buffer area within each catchment

        arcpy.analysis.SummarizeWithin(in_polygons=NR_catchments, in_sum_features=buff_out, out_feature_class=out_fc, shape_unit = "SQUAREKILOMETERS")
        
        # change field name to something that makes more sense
        arcpy.management.AlterField(out_fc, field = "sum_Area_SQUAREKILOMETERS", new_field_name = field_name, new_field_alias = field_alias)

        print(get_time(), f"Finished summarizing area in catchments")

if generate_HUC8:
    print(get_time(), "generating HUC8 codes in NHD table")
    # create HUC8 from reachcode
    arcpy.management.CalculateField(in_table = NHD_lyr, field = "HUC8", expression = "!REACHCODE![0:8]", expression_type = "PYTHON3", field_type = "TEXT")

if generate_90m_buff:
    # Let's try just generating buffer for all our HUC8s in a single file. What could possibly go wrong?


    print(get_time(), f"generating 90m buffer")
    # set output filepath
    temp_buff_out = os.path.join(buffer_gdb, f"NR_safo_01mBuff")
    merge_out = os.path.join(buffer_gdb, f"NR_safo_footprints")
    buff_out = os.path.join(buffer_gdb, f"NR_safo_90mBuff")

    

    # clear flowlines of any selections, just in case
    arcpy.SelectLayerByAttribute_management(NHD_lyr, "CLEAR_SELECTION")
    
    # select flowlines that are in native range HUC8s
    for HUC8 in NR_HUC8_list:
        NR_flowlines = arcpy.SelectLayerByAttribute_management(NHD_lyr,"ADD_TO_SELECTION", f"HUC8 = '{HUC8}'")
        NR_areas = arcpy.management.SelectLayerByLocation(NHD_areas_lyr, "INTERSECT", NR_flowlines)
        NR_waterbodies = arcpy.management.SelectLayerByLocation(NHD_waterbodies_lyr, "INTERSECT", NR_flowlines)

        arcpy.analysis.PairwiseBuffer(NR_flowlines, temp_buff_out, ".01 Meters", dissolve_option = "ALL")

        arcpy.management.Merge(inputs = [NR_areas, NR_waterbodies, temp_buff_out], output = merge_out)

    # generate 90m buffer for native range flowlines
    arcpy.analysis.PairwiseBuffer(merge_out, buff_out, "90 Meters", dissolve_option = "ALL")

    # summarize 90m buffer area within each catchment
    out_fc = os.path.join(buffer_gdb, f"NR_safo_90mBuff_sumwithin")
    catchment_summary(Catchments_lyr, buff_out, out_fc, "buffer90m_km2", "90m buffer area, sq km")



    """
    print(get_time(), f"generating 90m buffer for HUC {HUC8}")
    # set output filepath
    buff_out = os.path.join(buffer_gdb, f"HUC_{HUC8}_90mBuff")

    # select flowlines for HUC8
    HUC_flowlines = arcpy.SelectLayerByAttribute_management(NHD_lyr,"NEW_SELECTION", f"HUC8 = '{HUC8}'")

    # select catchments using flowlines, since catchments table doesn't include HUC8
    HUC_catchments = arcpy.management.SelectLayerByLocation(Catchments_lyr, "INTERSECT", HUC_flowlines)

    # generate 90m buffer for HUC8 flowlines
    arcpy.analysis.PairwiseBuffer(HUC_flowlines, buff_out, "90 Meters", dissolve_option = "ALL")

    print(get_time(), "summarizing buffer area within each catchment")

    # summarize 90m buffer area within each catchment
    out_fc = os.path.join(buffer_gdb, f"HUC_{HUC8}_90mBuff_sumwithin")
    arcpy.analysis.SummarizeWithin(in_polygons=HUC_catchments, in_sum_features=buff_out, out_feature_class=out_fc, shape_unit = "SQUAREKILOMETERS")

    # change field name to something that makes more sense
    arcpy.management.AlterField(out_fc, field = "sum_Area_SQUAREKILOMETERS", new_field_name = "buffer90m_km2", new_field_alias = "90m buffer area, sq km")

    print(get_time(), f"Finished summarizing buffer area in catchments")
"""


if generate_network_buffer_areas:

    # generate list of comids in HUC8 to use later
    # identify HUC catchments layer
    buffered_flowlines = os.path.join(buffer_gdb, f"NR_safo_90mBuff")
    NR_catchments_sumwithin = os.path.join(buffer_gdb, f"NR_safo_90mBuff_sumwithin")
    NR_catchments_lyr = arcpy.MakeFeatureLayer_management(NR_catchments_sumwithin, "Catchments_lyr")
    #HUC_catchments = os.path.join(buffer_gdb, f"HUC_{HUC8}_90mBuff_sumwithin")
    #HUC_catchments_lyr = arcpy.MakeFeatureLayer_management(HUC_catchments, "Catchments_lyr")

    # make list of all FEATUREID (COMID) values in HUC catchments layer
    comids = []
    with arcpy.da.SearchCursor(NR_catchments_lyr, ["FEATUREID"]) as cursor:
        for row in cursor:
            id = row[0]
            comids.append(id)

    print(get_time(), f"Native range has {len(comids)} stream reaches")

    # make empty list, to be filled with one dataframe for each comid
    df_list = []
    missing_comid_list = []

    with arcpy.da.SearchCursor(NR_catchments_lyr, ['FEATUREID', 'buffer90m_km2']) as cursor:
        reach_buffer_lookup = {row[0]:row[1] for row in cursor}

    for comid in tqdm(comids):
        #print(get_time(), f"Identifying reaches upstream of comid {comid}")
        upstream_comids_csv = rf"K:\GIS\AFWA_BrookTrout\Data\NHD_Tools_CSVs\{comid}.csv"
        df_path = os.path.join(root_dir, "Sensitivity_analysis", "buffers_90m", "NB_sqkm_tables", f"comid_{comid}_NB_area.csv")
        # note any comids that don't have a table of upstream reaches
        if not os.path.exists(upstream_comids_csv):
            missing_comid_list.append(comid)
            upstream_comids = [comid]

            
        else:
            
            

            # check if this reach has already been run. if so, skip calculation and just read in result. if not, do full calculations
            if os.path.exists(df_path): 
                df = pd.read_csv(df_path)
                df_list.append(df)
                continue

            # read in upstream comids, convert to numpy array for reasons that I've forgotten
            upstream_comids = pd.read_csv(upstream_comids_csv)['x'].to_numpy()

        upstream_areas_list = []
        for id in upstream_comids:
            # NOTE: stream reaches without a catchment cause an error due to mismatch between flowlines and catchment layers
            # Mismatch causes a key error that's handled below. 
            try:
                upstream_areas_list.append(reach_buffer_lookup[id])
            except KeyError:
                pass
                #print(f"lookup failed for id {id}")

                # TODO: generate a list of problematic comids, look at length etc
        upstream_buffer_area = sum(area for area in upstream_areas_list if area is not None)
        df = pd.DataFrame({'COMID': [comid], 'NB_areasqkm': [upstream_buffer_area], 'upstream_comids': [upstream_comids_csv]})
        df_list.append(df)
        # save as separate csv to figure out why mega-csv looks weird
        df.to_csv(df_path)

    # list of dfs becomes single df
    print(get_time(), f"making single table for native range")
    mega_df = pd.concat(df_list, ignore_index=True)

    # save mega df as separate csv
    mega_df_path = os.path.join(root_dir, "Sensitivity_analysis", "buffers_90m", "NB_sqkm_tables", f"NR_NB_area.csv")
    
    mega_df.to_csv(mega_df_path) 

    # save csv of missing comid values to investigate later
    missing_comid_data = {
        'COMID': missing_comid_list
    }
    missing_comid_df = pd.DataFrame(missing_comid_data)         
    missing_comid_CSV = os.path.join(root_dir, "Sensitivity_analysis", "buffers_90m", "NB_sqkm_tables", f"NR_missing_comids.csv")     
    missing_comid_df.to_csv(missing_comid_CSV)


# TODO: figure out how to run this for one HUC8 at a time
# TODO: then figure out how to iterate through multiple HUC8s


    # Output: data frame with 2 columns: comid & total network buffer area
    # append to list of dfs 
    # after for loop: smush together with pandas concat


if summarize_network_forest_buffer:
    # rasterize 90m buffers
    NR_safo_90mBuff_raster = os.path.join(buffer_gdb, f"NR_safo_90mBuff_rast")
    NLCD_90mBuff = os.path.join(buffer_gdb, f"NLCD_90mBuff")

    if not os.path.exists(NR_safo_90mBuff_raster):
        print("rasterizing 90m buffer polygons")
        arcpy.conversion.PolygonToRaster(
            in_features=os.path.join(buffer_gdb, f"NR_safo_90mBuff"),
            value_field = "OBJECTID",
            out_rasterdataset = NR_safo_90mBuff_raster,
            cellsize=NLCD_tif
        )

    # use rasterized 90m buffers as a mask on NLCD
    
    if not arcpy.Exists(NLCD_90mBuff):
        print(get_time(), "clipping NLCD raster to 90m buffers")

        # TODO: try rasterizing buffer polygons & using extract by mask instead of clip raster
        NLCD_90mBuff_extr = arcpy.sa.ExtractByMask(
            in_raster = NLCD_tif,
            in_mask_data = NR_safo_90mBuff_raster,
            extraction_area = "INSIDE"
        )

        NLCD_90mBuff_extr.save(NLCD_90mBuff)
        print(get_time(), "finished extracting NLCD raster for 90m buffer area")


    # reclass clipped NLCD raster to binary layer of just plantable classes
    # NLCD plantable classes: 21 (developed, open space); 22 (developed, low intensity); 81 (pasture/hay); 82 (cultivated crops)
    print(get_time(), "reclassifying NLCD buffer raster")
    NLCD_pa_remap = [[21, 1], [22, 1], [81, 1], [82, 1]]
    NLCD_pa_90mbuff_raster = arcpy.sa.Reclassify(
        in_raster = NLCD_90mBuff, 
        reclass_field = "Value", 
        remap = arcpy.sa.RemapValue(NLCD_pa_remap), 
        missing_values = "NODATA"
        )

    #summarize plantable areas by catchment
    out_fc = os.path.join(buffer_gdb, f"NR_safo_90mBuff_plantable_areas")
    field_name = "buffer90m_pa_km2"
    field_alias = "90m buffer plantable area, sq km"
    print(get_time(), "summarizing NLCD buffer by catchment")
    catchment_summary(Catchments_lyr, NLCD_pa_90mbuff_raster, out_fc, field_name, field_alias)
    print(get_time(), "finishing summarizing NLCD buffer by catchment")

    # TODO: turn NB area section into a function so we can just re-use that code here instead of copying everything



    # summarize plantable area by catchment similar to what I did for NB area calculation

    # generate network buffer plantable areas similar to what I did for NB area calculation
    



print(get_time(), "All done!")

# print how long it took script to run
end_time = time.time()
duration = end_time - start_time
print(f"Script ran for {duration/60} minutes")