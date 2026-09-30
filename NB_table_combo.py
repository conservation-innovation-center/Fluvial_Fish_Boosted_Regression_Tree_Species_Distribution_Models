# combine existing NB tables while figuring out how to get the rest to work

import os
import time
import sys
import pandas as pd
import numpy as np
import glob

def get_time():
    t = time.localtime()
    current_time = time.strftime("%H:%M:%S %m/%d/%Y", t)
    return current_time

print(get_time(),'Getting started!')
start_time = time.time()


folder_path = r"C:\Users\cweinstein\Documents\Projects\AFWA_2026\Sensitivity_analysis\buffers_90m\NB_sqkm_tables"

# Get all items, then filter out directories
# files = [f for f in os.listdir(folder_path)]
# print(files)

file_list = glob.iglob(r"C:\Users\cweinstein\Documents\Projects\AFWA_2026\Sensitivity_analysis\buffers_90m\NB_sqkm_tables\comid_*.csv", recursive=True)
#print(get_time(), f"trying to combine {len(file_list)} files")
df_list = []
# Yields files one by one instead of creating a massive list
for file_path in file_list:
    #print(file_path)
    df = pd.read_csv(file_path)
    df_list.append(df)

print(get_time(), "made list of dfs to combine")
# Merge them vertically
combined_df = pd.concat(df_list, ignore_index=True)

print(get_time(), "combined dfs!")
outpath = r"C:\Users\cweinstein\Documents\Projects\AFWA_2026\Sensitivity_analysis\buffers_90m\NB_sqkm_tables\NR_NB_area_test.csv"

combined_df.to_csv(outpath)

#print(combined_df)

print(get_time(), "All done!")

# print how long it took script to run
end_time = time.time()
duration = end_time - start_time
print(f"Script ran for {duration/60} minutes")